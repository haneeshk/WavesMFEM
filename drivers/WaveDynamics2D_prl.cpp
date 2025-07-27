//                        Dynamics of a 2D surface.
//                        Central difference scheme will be used

#include "mfem.hpp"
#include "omp.h"
#include <filesystem> // C++17
#include <nlohmann/json.hpp>
#include <fstream>
#include <iostream>
#include <string>
#include <iomanip> // at the top of your file if not already included
#include <cmath>   // Include cmath for pow function
#include <map>
#include <initializer_list>
#include <stdexcept>
#include "save.hpp"
#include "quess.hpp"
#include "mfemplus.hpp"

namespace fs = std::filesystem;
using namespace mfem;
using namespace std;
using json = nlohmann::json;

constexpr double π = M_PI;

struct matProps
{
	double λ;
	double μ;
	double ρ;
	double E;
	double ν;
};

struct simProps
{
	double t_final; // final t4ime
	double Δt;		// time step
	double dt_save;
	int N; // Number of time steps
	int n_save;
};

std::tuple<int, int, int> GetNodeInfo(const mfem::GridFunction *nodes_ptr);

int createNodeGridFunction(mfem::GridFunction &sel_nodes, const mfem::Array<int> &ess_tdof_list);

double hat_u_x_zero(const mfem::Vector &pt);
double hat_u_y_zero(const mfem::Vector &pt);
int determineDirichletDof(const mfem::ParMesh &pmesh, const mfem::ParFiniteElementSpace &pfespace, mfem::Array<int> &ess_tdof_listbcx, mfem::Array<int> &ess_tdof_listbcy);
int createDirichletVals(const GridFunction *nodes_ptr, mfem::Array<int> &ess_tdof_listbcx, mfem::Array<int> &ess_tdof_listbcy, std::function<double(const mfem::Vector &)> hat_u_x_func, std::function<double(const mfem::Vector &)> hat_u_y_func, mfem::Vector &hat_u_x, mfem::Vector &hat_u_y);
std::error_code logInWaveDynamics2D(json &inputParameters, matProps &brainMatProps, simProps &brainSimProps, mfem::Mesh *Ω);
int saveResults(json &inputParameters, ParaViewDataCollection &pvdc, std::ofstream &pointDisplacement, const int order, mfem::ParGridFunction &u, mfem::ParGridFunction &ϵ, mfem::ParGridFunction &σ, GridFunction &nodes);
FunctionCoefficient initialize_velocity(const json &inputParameters);
void InitializeOldSolution(const HypreParVector &u_hypre, const HypreParVector &v_hypre, const HypreParMatrix &M, const HypreParMatrix &K, const HypreParVector &F, double dt, HypreParVector &u_old_hypre);
void InitializeOldSolution(const ParGridFunction &u, const ParGridFunction &v, const HypreParMatrix &M, const HypreParMatrix &K, double dt, ParGridFunction &u_old);
bool readInputParameters(int argc, char *argv[], json &inputParameters);
double dirichletamp(const json &inputParameters, const double &t);
double tractionamp(const json &inputParameters, const double &t);

int main(int argc, char *argv[])
{
	// Initialize MPI and HYPRE.
	Mpi::Init(argc, argv);
	int num_procs = Mpi::WorldSize();
	int myid = Mpi::WorldRank();
	Hypre::Init();
	// Print OpenMP threads
	omp_set_num_threads(1);

	json inputParameters;
	matProps brainMatProps;
	simProps brainSimProps;
	readInputParameters(argc, argv, inputParameters);
	bool computeStressStrain = inputParameters.value("computeStressStrain", true);

	std::string mesh_file_str = inputParameters["Mesh Parameters"]["meshFileName"];
	const char *mesh_file = mesh_file_str.c_str();
	int order = inputParameters["Mesh Parameters"]["order"];
	bool static_cond = inputParameters["Mesh Parameters"]["static_cond"];
	bool visualization = inputParameters["Mesh Parameters"]["visualization"];
	int ref_levels = inputParameters["Mesh Parameters"]["ref_levels"];
	const char *device_config = "cpu";

	Device device(device_config);
	if (myid == 0)
		device.Print();

	// Read the mesh from the given mesh file on all processors.
	Mesh *mesh = new Mesh(mesh_file, 1, 1);
	int dim = mesh->Dimension(); // dim should equal 2.
	mesh->EnsureNodes();
	mesh->SetCurvature(order, false);	 // Enable high-order geometry
	for (int l = 0; l < ref_levels; l++) // Refine mesh
		mesh->UniformRefinement();

	// Define a parallel mesh on all MPI tasks by partitioning the serial mesh. Delete serial mesh.
	ParMesh *pmesh = new ParMesh(MPI_COMM_WORLD, *mesh);

	logInWaveDynamics2D(inputParameters, brainMatProps, brainSimProps, mesh);
	ConstantCoefficient lambda_coeff(brainMatProps.λ), mu_coeff(brainMatProps.μ);
	ConstantCoefficient E_coeff(brainMatProps.E), NU_coeff(brainMatProps.ν);
	ConstantCoefficient density(brainMatProps.ρ);

	// Data logged so delete mesh now.
	delete mesh;
	//------------------------------------------------------------------------------------------------------------------
	// Define parallel finite element spaces.

	FiniteElementCollection *fec = new H1_FECollection(order, dim); // Define a H1 finite element space for nodal displacements.
	ParFiniteElementSpace *pfespace = new ParFiniteElementSpace(pmesh, fec, dim);
	pmesh->SetNodalFESpace(pfespace);
	// fespacebc = new FiniteElementSpace(mesh, fec, 1); // Define a finite element space to project bcs in the loop.

	// All of these need to be checked, since pfespace and fespace work differently.
	cout << "Number of finite element unknowns in MPI_Rank: " << Mpi::WorldRank() << " is " << pfespace->GetTrueVSize() << endl;
	GridFunction *nodes = pmesh->GetNodes();
	MFEM_VERIFY(nodes, "Mesh must be created with generate_nodes = true");
	const int vdim = nodes->VectorDim();		 // Should be 2 for 2D
	const int num_nodes = pfespace->GetNVDofs(); // Num of nodes in each MPI rank.
	cout << "vdim" << vdim << endl;
	cout << "node->Size()" << nodes->Size() << endl; // Num of nodes in each MPI rank.
	cout << "num_nodes" << num_nodes << endl;		 // Num of nodes in each MPI rank.

	// Define an L2 finite element space for element strain and stress.
	// Dimension of the L2 finite element space is the number of strain and stress components.
	// For 2D, we have 3 components.

	FiniteElementCollection *l2fec;
	ParFiniteElementSpace *pl2fespace;

	l2fec = new L2_FECollection(0, dim);
	pl2fespace = new ParFiniteElementSpace(pmesh, l2fec, 3); // 3 since 3 strain and stress components
	const int num_els = pl2fespace->GetNDofs();

	Array<int> ess_bdr_x(pmesh->bdr_attributes.Max()), ess_bdr_y(pmesh->bdr_attributes.Max());

	ess_bdr_x = 0;
	ess_bdr_x[2] = 1; // right edge

	ess_bdr_y = 0;
	ess_bdr_y[2] = 1; // right edge

	Array<int> ess_tdof_listbcx, ess_tdof_listbcy;

	// These need to be checked too for parallel running.
	determineDirichletDof(*pmesh, *pfespace, ess_tdof_listbcx, ess_tdof_listbcy);
	Vector hat_u_x(ess_tdof_listbcx.Size()), hat_u_y(ess_tdof_listbcy.Size());
	createDirichletVals(nodes, ess_tdof_listbcx, ess_tdof_listbcy, hat_u_x_zero, hat_u_y_zero, hat_u_x, hat_u_y);

	// VectorFunctionCoefficient sel_nodes_coeff(2, [](const Vector &x, Vector &v)
	// 										  {
	// 											  v.SetSize(x.Size());
	// 											  v = 0.0;
	// 											  mfem::Vector c{2.5, 2.5}, pt(x);
	// 											  // 	for (int i = 0; i < ndofs; ++i)
	// 											  // 	{
	// 											  // 		std::cout << i;
	// 											  //
	// 											  // 		pt[0] = (*nodes)(i);
	// 											  // 		pt[1] = (*nodes)(i + ndofs);
	// 											  // 		pt -= c;
	// 											  //

	// 											  // Example: assign based on position x = (x[0], x[1])
	// 											  if (pt.Norml2() > 2.5){
	// 											  v[0] = 1.0; // x-component
	// 											  v[1] = 1.0; // y-component
	// 											  } });

	ParGridFunction sel_nodes(pfespace);
	sel_nodes = 0.0;

	for (int i = 0; i < ess_tdof_listbcx.Size(); i++)
	{
		int idx = ess_tdof_listbcx[i];
		std::cout << "idx :" << idx << "\n";
		sel_nodes(idx) = 1.0;
	}

	// Force vector if non-zero tractions are applied.
	VectorArrayCoefficient traction(dim);
	Vector tr(pmesh->bdr_attributes.Size());
	tr = 0.0;
	cout << tractionamp(inputParameters, 0.0) << endl;
	tr(0) = tractionamp(inputParameters, 0.0); // boundary attribute.

	traction.Set(0, new PWConstCoefficient(tr));   // x component
	traction.Set(1, new ConstantCoefficient(0.0)); // y component

	ParLinearForm *f;
	f = new ParLinearForm(pfespace);
	f->AddBoundaryIntegrator(new VectorBoundaryLFIntegrator(traction));
	f->Assemble();
	HypreParVector F = *(f->ParallelAssemble());
	delete f;

	// Mass and stiffness matrices
	ParBilinearForm m(pfespace);
	m.AddDomainIntegrator(new VectorMassIntegrator(density));
	m.Assemble();
	m.Finalize();
	HypreParMatrix *M = (m.ParallelAssemble());

	HypreBoomerAMG *amg = new HypreBoomerAMG();
	amg->SetSystemsOptions(dim);
	amg->SetPrintLevel(0);
	amg->SetOperator(*M);

	HyprePCG *pcg = new HyprePCG(MPI_COMM_WORLD);
	pcg->SetOperator(*M);
	pcg->SetTol(1e-15);
	pcg->SetMaxIter(500);
	pcg->SetPrintLevel(0);
	pcg->SetPreconditioner(*amg);

	ParBilinearForm k(pfespace);
	k.AddDomainIntegrator(new mfemplus::IsotropicElasticityIntegrator(E_coeff, NU_coeff));
	k.Assemble();
	k.Finalize();
	HypreParMatrix *K = (k.ParallelAssemble());

	ParGridFunction u(pfespace), u_old(pfespace), u_new(pfespace), v(pfespace), a(pfespace);
	u = u_old = u_new = v = a = 0.0;
	HypreParVector u_hypre(pfespace), u_old_hypre(pfespace), u_new_hypre(pfespace), v_hypre(pfespace), a_hypre(pfespace), rhs(pfespace);
	u_hypre = u_old_hypre = u_new_hypre = v_hypre = a_hypre = rhs = 0.0;
	// Need to check these parts. What will rhs llok like if it is initialized as it currently is?
	// InitializeOldSolution(u, v, M, K, dt, u_old); // comment if non-zero tractions are prescribed.
	InitializeOldSolution(u_hypre, v_hypre, *M, *K, F, brainSimProps.Δt, u_old_hypre); // un-comment if non-zero tractions are prescribed.

	// for (int i = 0; i < u_hypre.Size(); i++)
	// {
	// 	u_hypre(i) = u(i);
	// 	u_old_hypre(i) = u_old(i);
	// 	u_new_hypre(i) = u_new(i);
	// 	v_hypre(i) = v(i);
	// 	a_hypre(i) = a(i);
	// }

	// Initialize GlobalStressStrain object. Initialize stress and strain grid functions.
	mfemplus::GlobalStressStrain WaveDynamics(pmesh, pfespace);
	ParGridFunction eps(pl2fespace), sig(pl2fespace);
	if (computeStressStrain)
	{
		eps = sig = 0.0;
		WaveDynamics.GlobalStrain(u_new, eps);
		WaveDynamics.GlobalStress(eps, E_coeff, NU_coeff, sig);
	}

	// Inititialize Paraview object
	ParaViewDataCollection pvdc("Waves2D", pmesh);
	ofstream pointDisplacement;

	saveResults(inputParameters, pvdc, pointDisplacement, order, u, eps, sig, sel_nodes);

	const int selectNode = 2;
	int cycle = 0;
	int t = 0;
	pvdc.SetCycle(cycle); // Record time step number
	pvdc.SetTime(t);	  // Record simulation time
	pvdc.Save();

	cycle += 1;

	// Time-stepping loop
	cout << "Time stepping loop beginning" << endl;

	for (double t = brainSimProps.Δt; t <= brainSimProps.t_final; t += brainSimProps.Δt)
	{

		tr = 0.0;
		tr(0) = tractionamp(inputParameters, t);
		traction.Set(0, new PWConstCoefficient(tr));   // x component
		traction.Set(1, new ConstantCoefficient(0.0)); // y component

		ParLinearForm f(pfespace);
		f.AddBoundaryIntegrator(new VectorBoundaryLFIntegrator(traction));
		f.Assemble();
		HypreParVector F = *(f.ParallelAssemble());
		// rhs = 0.0;
		K->Mult(u_hypre, rhs); // rhs = F - K u
		// rhs *= -1; // comment this if prescribing non-zero tractions.
		add(F, -1.0, rhs, rhs); // un-comment this if prescribing non-zero tractions.

		pcg->Mult(rhs, a_hypre); // solve M a = rhs
		// Central difference: u_{n+1} = 2u_n - u_{n-1} + brainSimProps.Δt^2 * a
		HypreParVector temp(pfespace);
		temp = 0.0;
		add(u_hypre, u_hypre, temp);										  // 2u_n
		add(temp, -1.0, u_old_hypre, temp);									  // 2u_n - u_{n-1}
		add(temp, brainSimProps.Δt * brainSimProps.Δt, a_hypre, u_new_hypre); // u_{n+1}

		{
			int i = 0;
			for (auto nd : ess_tdof_listbcx)
			{
				u_new_hypre(nd) = hat_u_x(i);
				i++;
			}

			i = 0;
			for (auto nd : ess_tdof_listbcy)
			{
				u_new_hypre(nd + num_nodes) = hat_u_y(i);
				i++;
			}
		}

		u_old_hypre = u_hypre;
		u_hypre = u_new_hypre;

		u_old = u_old_hypre;
		u = u_hypre;
		u_new = u_new_hypre;

		if (computeStressStrain)
		{
			WaveDynamics.GlobalStrain(u_new, eps);
			WaveDynamics.GlobalStress(eps, E_coeff, NU_coeff, sig);

			// Convert engineering shear strain to true shear strain.
			for (int i = (2 * num_els) - 1; i < eps.Size(); i++)
			{
				eps(i) = eps(i) / 2;
			}
		}

		if ((cycle % brainSimProps.n_save) == 0)
		{

			pointDisplacement << std::fixed << std::setprecision(10) << t << "\t" << u(selectNode) << "\n";
			pvdc.SetCycle(cycle); // Record time step number
			pvdc.SetTime(t);	  // Record simulation time
			pvdc.Save();
			cout << "cycle\t:" << cycle << "/" << brainSimProps.N << endl;
		}
		cycle++;
	}

	pointDisplacement.close();
	ofstream finaldisp("../results/" + inputParameters["testName"].get<std::string>() + "/finaldisp_prl.dat");
	u.Save(finaldisp);
	finaldisp.close();

	cout << "Output saved. " << endl;

	delete pfespace;
	delete fec;
	delete pmesh;

	return 0;
}

FunctionCoefficient initialize_velocity(const json &inputParameters)
{
	// Extract required parameters
	double A = inputParameters["Boundary Conditions"]["A"];
	double epsilon = inputParameters["Initial Conditions"]["epsilon"];
	double omega = inputParameters["Physical Parameters"]["omega"];

	auto epsilon_mollifier = [epsilon](double x) -> double
	{
		if (std::abs(x) >= epsilon)
			return 0.0;
		double r = x / epsilon;
		double peak_inv = std::exp(1.0);
		return std::exp(-1.0 / (1.0 - r * r)) * peak_inv;
	};

	// Define and return the FunctionCoefficient
	return FunctionCoefficient([A, omega, epsilon_mollifier](const Vector &x) -> double
							   { return A * omega * epsilon_mollifier(x[0]); });
}
// Initializes u_old = u - dt * v + 0.5 + dt² * (M⁻¹ (f - K u))
void InitializeOldSolution(const HypreParVector &u_hypre, const HypreParVector &v_hypre, const HypreParMatrix &M, const HypreParMatrix &K, const HypreParVector &F, double dt, HypreParVector &u_old_hypre)
{
	// MFEM_ASSERT(u.FESpace() == u_old.FESpace(), "u and u_old must use same FESpace");

	HypreParVector Ku(F), rhs(F), a0(F);
	Ku = rhs = a0 = 0.0;

	K.Mult(u_hypre, Ku);
	add(F, -1.0, Ku, rhs);
	// Ku *= -1;

	HypreBoomerAMG *amg = new HypreBoomerAMG(M);
	// amg->SetSystemsOptions(u.FESpace()->GetVDim());
	amg->SetPrintLevel(0);

	HyprePCG *pcg = new HyprePCG(MPI_COMM_WORLD);
	pcg->SetOperator(M);
	pcg->SetTol(1e-12);
	pcg->SetMaxIter(500);
	pcg->SetPrintLevel(0);
	pcg->SetPreconditioner(*amg);
	pcg->Mult(rhs, a0); // a0 = M^{-1} (F - K u)

	u_old_hypre = u_hypre;
	u_old_hypre.Add(-dt, v_hypre);
	u_old_hypre.Add(0.5 * dt * dt, a0); // ✅ correct sign
}
// Initializes u_old = u - dt * v + 0.5 - dt² * (M⁻¹K u)
void InitializeOldSolution(const ParGridFunction &u, const ParGridFunction &v, const HypreParMatrix &M, const HypreParMatrix &K, double dt, ParGridFunction &u_old)
{
	MFEM_ASSERT(u.FESpace() == u_old.FESpace(), "u and u_old must use same FESpace");

	HypreParVector Ku(u.ParFESpace()), a0(u.ParFESpace());
	K.Mult(u, Ku);
	a0 = 0.0;

	HypreBoomerAMG amg(M);
	amg.SetSystemsOptions(u.FESpace()->GetVDim());
	amg.SetPrintLevel(0);

	HyprePCG pcg(MPI_COMM_WORLD);
	pcg.SetOperator(M);
	pcg.SetTol(1e-12);
	pcg.SetMaxIter(500);
	pcg.SetPrintLevel(0);
	pcg.SetPreconditioner(amg);
	pcg.Mult(Ku, a0); // a0 = M^{-1} (- K u)
	u_old = u;
	u_old.Add(-dt, v);
	u_old.Add(-0.5 * dt * dt, a0); // ✅ correct sign
}
// This function needs to be modified to prescribe appropriate Dirichlet BCs. Make sure parameters are in JSON file.
double dirichletamp(const json &inputParameters, const double &t)
{
	double A = inputParameters["Boundary Conditions"]["A"];
	double omega = inputParameters["Physical Parameters"]["omega"];

	return A * sin(omega * t);
}

double tractionamp(const json &inputParameters, const double &t)
{
	double T = inputParameters["Boundary Conditions"]["T"];
	double t0 = inputParameters["Boundary Conditions"]["t0"];
	double zero = 0.0;

	if (t < t0)
		return T;
	else
		return 0;
}

bool readInputParameters(int argc, char *argv[], json &inputParameters)
{
	std::cout << "Inside the readInputParameters file" << std::endl;
	//  (argc > 1) ? argv[1] :
	std::string json_path = "../input/InputParameters/WaveDynamics2D.json"; // Change JSON file path accordingly.

	std::ifstream infile(json_path);
	if (!infile.is_open())
	{
		std::cerr << "Error: Could not open " << json_path << std::endl;
		return 1;
	}

	infile >> inputParameters;

	return 0;
}

int saveResults(json &inputParameters, ParaViewDataCollection &pvdc, std::ofstream &pointDisplacement, const int order, mfem::ParGridFunction &u, mfem::ParGridFunction &ϵ, mfem::ParGridFunction &σ, mfem::GridFunction &nodes)
{
	std::string resultsFolder = "../results/" + inputParameters["testName"].get<std::string>();
	std::string pointDataFile = resultsFolder + "/" + "pointDisplacement.dat";

	pointDisplacement.open(pointDataFile);
	if (!pointDisplacement.is_open())
	{
		std::cerr << "Failed to open: " << pointDataFile << std::endl;
		return 1; // or some error code
	}

	pvdc.SetPrefixPath(resultsFolder); // Directory to save data
	pvdc.SetLevelsOfDetail(order);	   // Optional: for visualization
	pvdc.SetHighOrderOutput(true);	   // Keep high-order info
	pvdc.RegisterField("u", &u);	   // Associate displacement field with data collection

	pvdc.RegisterField("nodes", &nodes);

	bool computeStressStrain = inputParameters.value("computeStressStrain", true);
	if (computeStressStrain)
	{
		pvdc.RegisterField("strain", &ϵ);
		pvdc.RegisterField("stress", &σ);
	}

	return 0;
}

std::error_code logInWaveDynamics2D(json &inputParameters, matProps &brainMatProps, simProps &brainSimProps, mfem::Mesh *Ω)
{

	std::string resultsFolder = "../results/" + inputParameters["testName"].get<std::string>();
	if (fs::exists(resultsFolder))
	{
		std::error_code ec;
		fs::remove_all(resultsFolder, ec); // remove all contents recursively
		if (ec)
		{
			std::cerr << "Error removing folder: " << ec.message() << std::endl;
			return ec;
		}
	}

	// Recreate a fresh empty directory
	fs::create_directories(resultsFolder);

	brainMatProps.λ = inputParameters["Physical Parameters"]["lambda"];
	brainMatProps.μ = inputParameters["Physical Parameters"]["mu"];

	brainMatProps.E = inputParameters["Physical Parameters"]["youngmod"];
	brainMatProps.ν = inputParameters["Physical Parameters"]["poissonratio"];
	brainMatProps.ρ = inputParameters["Physical Parameters"]["Density"];
	std::cout << "the density is" << brainMatProps.ρ << "\n";

	cout << "Physical parameters read" << endl;

	brainSimProps.t_final = inputParameters["Simulation Parameters"]["Total Simulation Time"];
	brainSimProps.Δt = inputParameters["Simulation Parameters"]["Simulation Time Step"];
	brainSimProps.N = static_cast<int>(std::round(brainSimProps.t_final / brainSimProps.Δt));
	brainSimProps.dt_save = inputParameters["Simulation Parameters"]["Data Storage Time Interval"];
	brainSimProps.n_save = static_cast<int>(std::round(brainSimProps.dt_save / brainSimProps.Δt));
	cout << "Simulation parameters read" << endl;

	// std::cout << "dt_save / dt" << dt_save / dt << std::endl;
	// std::cout << "n_save\t:" << n_save << std::endl;

	ofstream u_data(resultsFolder + "/" + "SimulationDetails.txt");
	{

		double λ, μ, ρ;
		λ = brainMatProps.λ;
		μ = brainMatProps.μ;
		ρ = brainMatProps.ρ;

		double h = Ω->GetElementSize(1, /*type=*/1);
		double c_p = std::sqrt((λ + 2 * μ) / ρ);
		double c_s = std::sqrt(μ / ρ);
		double CFL_dt = h / c_p;
		u_data << "saved output to " + resultsFolder << endl;
		u_data << "t_final\t" << brainSimProps.t_final << endl;
		u_data << "dt\t" << brainSimProps.Δt << endl;
		u_data << "c_p\t" << c_p << endl;
		u_data << "c_s\t" << c_s << endl;
		u_data << "CFL time\t" << CFL_dt << endl;
		u_data << "Time Steps \t" << brainSimProps.N << endl;
		u_data << "Data storage time interval \t" << brainSimProps.dt_save << endl;
		u_data << "\n\n";
		u_data << "Mesh Details";
		u_data << "\n\n";
		u_data << "order: " << inputParameters["Mesh Parameters"]["order"] << "\n";
		u_data << "ref_levels: " << inputParameters["Mesh Parameters"]["ref_levels"] << "\n";
		u_data << "Number of elements: " << Ω->GetNE() << "\n"; // This will only return number of elements in a specific task.
	}

	u_data.close();
	return {};
};

int determineDirichletDof(const mfem::ParMesh &pmesh, const mfem::ParFiniteElementSpace &pfespace,
						  mfem::Array<int> &ess_tdof_listbcx, mfem::Array<int> &ess_tdof_listbcy)
{

	mfem::Array<int> all_bdr_attr_selector(pmesh.bdr_attributes.Max());
	all_bdr_attr_selector = 1;

	mfem::Array<int> bdr_nodes_list;
	pfespace.GetEssentialTrueDofs(all_bdr_attr_selector, bdr_nodes_list, 0);
	const mfem::GridFunction *nodes_ptr = pmesh.GetNodes();
	// pmesh.GetDof

	auto [num_dofs, vdim, num_nodes] = GetNodeInfo(nodes_ptr);

	// MFEM_VERIFY(nodes, "Mesh has no nodes. Did you call SetCurvature or SetNodalFESpace?");
	mfem::Vector c{2.5, 2.5};
	cout << "nodes_ptr.Size(): " << nodes_ptr->Size() << "\n";
	double ϵ = 0.0001;
	for (auto nd : bdr_nodes_list)
	{
		mfem::Vector pt(2);

		pt[0] = (*nodes_ptr)(nd);
		pt[1] = (*nodes_ptr)(nd + num_nodes);

		if (std::abs(pt[0] - 5.0) < ϵ /*&& std::abs(pt[1] - 2.5) < 1.0*/)
		{
			std::cout << "nd is: " << nd << "{" << pt[0] << "," << pt[1] << "}," << std::endl;
			ess_tdof_listbcx.Append(nd);
		}

		if (std::abs(pt[0] - 5.0) < ϵ /*&& std::abs(pt[1] - 2.5) < 1.0*/)
		{
			std::cout << "nd is: " << nd << "{" << pt[0] << "," << pt[1] << "}," << std::endl;

			ess_tdof_listbcy.Append(nd);
		}
	}
	return 0;
}

int createNodeGridFunction(mfem::GridFunction &sel_nodes, const mfem::Array<int> &ess_tdof_list)
{

	for (auto i : ess_tdof_list)
	{
		sel_nodes(i) = 1.0;
	}

	return 0;
}

double hat_u_x_zero(const mfem::Vector &pt)
{
	double x1 = pt[0];
	double x2 = pt[1];
	return 0.0;
}

double hat_u_y_zero(const mfem::Vector &pt)
{
	double x1 = pt[0];
	double x2 = pt[1];
	return 0.0;
}

int createDirichletVals(const GridFunction *nodes_ptr, mfem::Array<int> &ess_tdof_listbcx, mfem::Array<int> &ess_tdof_listbcy, std::function<double(const mfem::Vector &)> hat_u_x_func, std::function<double(const mfem::Vector &)> hat_u_y_func, mfem::Vector &hat_u_x, mfem::Vector &hat_u_y)
{
	auto [num_dofs, vdim, num_nodes] = GetNodeInfo(nodes_ptr);
	mfem::Vector pt{0.0, 0.0};
	int i = 0;
	for (auto nd : ess_tdof_listbcx)
	{

		pt[0] = (*nodes_ptr)(nd);
		pt[1] = (*nodes_ptr)(nd + num_nodes);
		hat_u_x(i) = hat_u_x_func(pt);
		i++;
	}
	i = 0;
	for (auto nd : ess_tdof_listbcy)
	{

		pt[0] = (*nodes_ptr)(nd);
		pt[1] = (*nodes_ptr)(nd + num_nodes);
		hat_u_y(i) = hat_u_y_func(pt);
		i++;
	}

	return 0;
}

std::tuple<int, int, int> GetNodeInfo(const mfem::GridFunction *nodes_ptr)
{
	if (!nodes_ptr)
	{
		throw std::invalid_argument("Error: nodes_ptr is null.");
	}

	int num_dofs = nodes_ptr->Size();  // total number of entries
	int vdim = nodes_ptr->VectorDim(); // dimension of coordinates
	int num_nodes = num_dofs / vdim;   // number of actual nodes

	return std::make_tuple(num_dofs, vdim, num_nodes);
}

// DONE: clean yo "createNodeGridFunction". There is a lot of commented stuff in it.

// TODO: Remove the u_y displacement boundary condition. However, now the problem has no constraings in the y direction.
// TODO: need to be able to start a simulation from the results of a previous simulation.
// TODO: Change geometry to CDisk.
// TODO: Change geometry to Disk with two holes.
// TODO: Add spatially varing Elastic Properties.
