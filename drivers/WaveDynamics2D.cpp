//                        Dynamics of a 2D surface.
//                        Central difference scheme will be used

#include "mfem.hpp"
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

void JacobiEigenvalues(const mfem::DenseMatrix &A, mfem::Vector &eigenvalues, int max_iter = 100, double tol = 1e-12)
{
	const int n = A.Height();
	MFEM_VERIFY(A.Width() == n, "Matrix must be square");

	mfem::DenseMatrix M(A); // make a copy to work on
	eigenvalues.SetSize(n);

	for (int iter = 0; iter < max_iter; ++iter)
	{
		// Find largest off-diagonal element
		int p = 0, q = 1;
		double max_val = std::abs(M(0, 1));
		for (int i = 0; i < n; ++i)
		{
			for (int j = i + 1; j < n; ++j)
			{
				double val = std::abs(M(i, j));
				if (val > max_val)
				{
					max_val = val;
					p = i;
					q = j;
				}
			}
		}

		// Converged?
		if (max_val < tol)
			break;

		double θ = 0.5 * std::atan2(2.0 * M(p, q), M(q, q) - M(p, p));
		double c = std::cos(θ);
		double s = std::sin(θ);

		// Rotate
		for (int i = 0; i < n; ++i)
		{
			if (i != p && i != q)
			{
				double Mip = M(i, p);
				double Miq = M(i, q);
				M(i, p) = M(p, i) = c * Mip - s * Miq;
				M(i, q) = M(q, i) = s * Mip + c * Miq;
			}
		}

		double Mpp = M(p, p);
		double Mqq = M(q, q);
		double Mpq = M(p, q);
		M(p, p) = c * c * Mpp - 2 * s * c * Mpq + s * s * Mqq;
		M(q, q) = s * s * Mpp + 2 * s * c * Mpq + c * c * Mqq;
		M(p, q) = M(q, p) = 0.0;
	}

	// Fill in diagonal as eigenvalues
	for (int i = 0; i < n; ++i)
		eigenvalues[i] = M(i, i);
}
int createNodeGridFunction(mfem::GridFunction &sel_nodes, const mfem::Array<int> &ess_tdof_list);

int determineDirichletDof(const mfem::Mesh &mesh,
						  const mfem::FiniteElementSpace &fespace,
						  mfem::Array<int> &ess_tdof_listx,
						  mfem::Array<int> &ess_tdof_listy);
std::error_code logInWaveDynamics2D(json &inputParameters, matProps &brainMatProps, simProps &brainSimProps, mfem::Mesh *Ω);
int saveResults(json &inputParameters, ParaViewDataCollection &pvdc, std::ofstream &pointDisplacement, const int order, mfem::GridFunction &u, mfem::GridFunction &ϵ, mfem::GridFunction &σ, GridFunction &nodes);
FunctionCoefficient initialize_velocity(const json &inputParameters);
void InitializeOldSolution(const GridFunction &u, const GridFunction &v, const SparseMatrix &M, const SparseMatrix &K, const Vector &F, double dt, GridFunction &u_old);
void InitializeOldSolution(const GridFunction &u, const GridFunction &v, const SparseMatrix &M, const SparseMatrix &K, double dt, GridFunction &u_old);
bool readInputParameters(int argc, char *argv[], json &inputParameters);
double dirichletamp(const json &inputParameters, const double &t);
double tractionamp(const json &inputParameters, const double &t);

int main(int argc, char *argv[])
{
	json inputParameters;
	matProps brainMatProps;
	simProps brainSimProps;
	readInputParameters(argc, argv, inputParameters);
	bool computeStressStrain = inputParameters.value("computeStressStrain", true);
	//   Parse command-line options.
	// const char *mesh_file = "./input/meshes/Rectangle-quad.mesh"; // Change this line, include 2D mesh file.
	std::string mesh_file_str = inputParameters["Mesh Parameters"]["meshFileName"];
	const char *mesh_file = mesh_file_str.c_str();
	int order = inputParameters["Mesh Parameters"]["order"];
	bool static_cond = inputParameters["Mesh Parameters"]["static_cond"];
	bool visualization = inputParameters["Mesh Parameters"]["visualization"];
	int ref_levels = inputParameters["Mesh Parameters"]["ref_levels"];

	// Read the mesh from the given mesh file.
	Mesh *mesh = new Mesh(mesh_file, 1, 1);
	int dim = mesh->Dimension();	  // dim should equal 2.
	mesh->SetCurvature(order, false); // Enable high-order geometry
	for (int l = 0; l < ref_levels; l++)
		mesh->UniformRefinement();

	logInWaveDynamics2D(inputParameters, brainMatProps, brainSimProps, mesh);
	ConstantCoefficient lambda_coeff(brainMatProps.λ), mu_coeff(brainMatProps.μ);
	ConstantCoefficient E_coeff(brainMatProps.E), NU_coeff(brainMatProps.ν);
	ConstantCoefficient density(brainMatProps.ρ);

	//------------------------------------------------------------------------------------------------------------------
	// Define finite element spaces.

	FiniteElementCollection *fec;
	FiniteElementSpace *fespace, *fespacebc;

	fec = new H1_FECollection(order, dim); // Define a H1 finite element space for nodal displacements.
	fespace = new FiniteElementSpace(mesh, fec, dim);
	mesh->SetNodalFESpace(fespace);
	cout << "Mesh created" << endl;

	fespacebc = new FiniteElementSpace(mesh, fec, 1); // Define a finite element space to project bcs in the loop.

	cout << "Number of finite element unknowns: " << fespace->GetTrueVSize() << endl;
	GridFunction *nodes = mesh->GetNodes();
	MFEM_VERIFY(nodes, "Mesh must be created with generate_nodes = true");
	const int vdim = nodes->VectorDim(); // Should be 2 for 2D
	const int ndofs = nodes->FESpace()->GetNDofs();
	cout << "vdim" << vdim << endl;
	cout << "node->Size()" << nodes->Size() << endl;
	cout << "ndofs" << ndofs << endl;

	// Define an L2 finite element space for element strain and stress.
	// Dimension of the L2 finite element space is the number of strain and stress components.
	// For 2D, we have 3 components.

	FiniteElementCollection *l2fec;
	FiniteElementSpace *l2fespace;

	l2fec = new L2_FECollection(0, dim);
	l2fespace = new FiniteElementSpace(mesh, l2fec, 3);
	const int num_els = l2fespace->GetNDofs();

	//------------------------------------------------------------------------------------------------------------------

	//------------------------------------------------------------------------------------------------------------------
	//    Determine the list of essential boundary dofs for each vector dimension.

	Array<int> ess_bdr_x(mesh->bdr_attributes.Max()), ess_bdr_y(mesh->bdr_attributes.Max());

	ess_bdr_x = 0;
	ess_bdr_x[2] = 1; // right edge

	ess_bdr_y = 0;
	ess_bdr_y[2] = 1; // right edge

	// fespace->GetEssentialTrueDofs(ess_bdr_x, ess_tdof_listx, 0);
	// fespace->GetEssentialTrueDofs(ess_bdr_y, ess_tdof_listy, 1);

	// To project BCs in the loop, define scalar tdof lists for each component.
	Array<int> ess_tdof_listbcx, ess_tdof_listbcy;
	fespacebc->GetEssentialTrueDofs(ess_bdr_x, ess_tdof_listbcx, 0);
	fespacebc->GetEssentialTrueDofs(ess_bdr_y, ess_tdof_listbcy, 0);
	// determineDircheletDof(mesh, ess_tdof_listbcx, ess_tdof_listbcy);
	// determineDirichletDof(*mesh, *fespace, ess_tdof_listbcx, ess_tdof_listbcy);
	// {

	// 	if (!nodes)
	// 	{
	// 		std::cerr << "Mesh has no nodes! Did you forget SetCurvature()?" << std::endl;
	// 		std::cout << "Press Enter to continue...";
	// 		std::cin.get(); // Waits until Enter is pressed
	// 	}
	// 	else
	// 	{
	// 		std::cout << "Nodes exists";
	// 		std::cout << "Press Enter to continue...";
	// 		std::cin.get(); // Waits until Enter is pressed
	// 	}
	// }

	// Define vectors for each displacement to apply BCs in the loop.
	Vector dispbcx(mesh->bdr_attributes.Size()), dispbcy(mesh->bdr_attributes.Size());
	//------------------------------------------------------------------------------------------------------------------

	//------------------------------------------------------------------------------------------------------------------
	// Create and initialize displacement, velocity, and acceleration grid functions.
	GridFunction u(fespace), v(fespace), a(fespace);

	// Initial conditions

	FunctionCoefficient u0_coeff([](const Vector &x)
								 { return 0.0; });
	FunctionCoefficient v0_coeff([](const Vector &x)
								 { return 0.0; });
	// u.ProjectCoefficient(u0_coeff);
	v.ProjectCoefficient(v0_coeff);
	u = 0.0;
	v = 0.0;
	a = 0.0;
	{
		cout << "the norms of u, v, and a are" << u.Norml2() << ", " << v.Norml2() << ", " << a.Norml2() << "\n"; // L2 (Euclidean) norm
	}

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

	GridFunction sel_nodes(fespacebc);
	sel_nodes = 0.0;
	// for (int n = 0; n < ndofs; n++)
	for (int i = 0; i < ess_tdof_listbcx.Size(); i++)
	{
		int idx = ess_tdof_listbcx[i];
		// std::cout << "idx :" << idx << "\n";
		sel_nodes(idx) = 1.0;
	}

	// for (auto i : ess_tdof_listbcx)
	// {
	// 	std::cout << i << "\n";
	// 	sel_nodes(i) = 1.0;
	// 	std::cout << sel_nodes(i) << "\n";
	// }
	// sel_nodes.ProjectCoefficient(sel_nodes_coeff);
	// createNodeGridFunction(sel_nodes, ess_tdof_list);

	// Force vector if non-zero tractions are applied.
	VectorArrayCoefficient traction(dim);
	Vector tr(mesh->bdr_attributes.Size());
	tr = 0.0;
	cout << tractionamp(inputParameters, 0.0) << endl;
	tr(0) = tractionamp(inputParameters, 0.0); // boundary attribute.

	traction.Set(0, new PWConstCoefficient(tr));   // x component
	traction.Set(1, new ConstantCoefficient(0.0)); // y component

	LinearForm *f;
	f = new LinearForm(fespace);
	f->AddBoundaryIntegrator(new VectorBoundaryLFIntegrator(traction));
	f->Assemble();
	Vector F(*f);
	delete f;

	// Mass and stiffness matrices
	BilinearForm m(fespace);
	m.AddDomainIntegrator(new VectorMassIntegrator(density));
	m.Assemble();
	m.Finalize();
	SparseMatrix M(m.SpMat());

	CGSolver cg;
	GSSmoother prec(M);
	cg.SetPreconditioner(prec);
	cg.SetOperator(M);
	cg.SetRelTol(1e-12);
	cg.SetAbsTol(1e-15);
	cg.SetMaxIter(200);
	cg.SetPrintLevel(0);

	BilinearForm k(fespace);
	k.AddDomainIntegrator(new mfemplus::IsotropicElasticityIntegrator(E_coeff, NU_coeff));
	// k.AddDomainIntegrator(new ElasticityIntegrator(lambda_coeff, mu_coeff));
	k.Assemble();
	k.Finalize();
	SparseMatrix K(k.SpMat());

	GridFunction u_old(fespace);
	// InitializeOldSolution(u, v, M, K, dt, u_old); // comment if non-zero tractions are prescribed.
	InitializeOldSolution(u, v, M, K, F, brainSimProps.Δt, u_old); // un-comment if non-zero tractions are prescribed.

	GridFunction u_new(u);
	Vector rhs(fespace->GetVSize()); // Full-sized RHS
	rhs = 0.0;
	ConstantCoefficient zero(0.0);

	// Initialize GlobalStressStrain object. Initialize stress and strain grid functions.
	mfemplus::GlobalStressStrain WaveDynamics(mesh, fespace);
	GridFunction eps(l2fespace), sig(l2fespace);
	if (computeStressStrain)
	{
		eps = sig = 0.0;
		WaveDynamics.GlobalStrain(u_new, eps);
		WaveDynamics.GlobalStress(eps, E_coeff, NU_coeff, sig);
	}

	// Inititialize Paraview object
	ParaViewDataCollection pvdc("Waves2D", mesh);
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
		// std::cout << t << std::endl;

		// Specify and assemble traction on rhs.
		tr = 0.0;
		tr(0) = tractionamp(inputParameters, t);

		traction.Set(0, new PWConstCoefficient(tr));   // x component
		traction.Set(1, new ConstantCoefficient(0.0)); // y component

		LinearForm *f;
		f = new LinearForm(fespace);
		f->AddBoundaryIntegrator(new VectorBoundaryLFIntegrator(traction));
		f->Assemble();
		Vector F(*f);
		delete f;

		// rhs = F - K u
		K.Mult(u, rhs);
		// rhs *= -1; // comment this if prescribing non-zero tractions.
		add(F, -1.0, rhs, rhs); // un-comment this if prescribing non-zero tractions.

		cg.Mult(rhs, a); // solve M a = rhs
		// Central difference: u_{n+1} = 2u_n - u_{n-1} + brainSimProps.Δt^2 * a
		Vector temp(u.Size());
		add(u, u, temp);										  // 2u_n
		add(temp, -1.0, u_old, temp);							  // 2u_n - u_{n-1}
		add(temp, brainSimProps.Δt * brainSimProps.Δt, a, u_new); // u_{n+1}

		// Project Dirichlet BCs back onto GridFunction. Each vector component on a bdr_attribute can have different BC.
		// A rather convoluted approach is taken here since I could not find a way with MFEM's data structures.
		// Perhaps a neater way can be found in the future.

		dispbcx = dispbcy = 0.0; // Vectors to prescribe Dirichlet BCs by bdr_attribute
		// dispbcx(1) = dirichletamp(inputParameters, t); // Constant Non-zero Dirichlet BCs only on bdr_attributes in the x/z component.

		// Lambda functions to prescribe spatially dependent Dirichlet BCs.

		// FunctionCoefficient dircx([t, &mu, &rho, &l, &B1, &eps](const Vector &xcoords)
		//                           {
		//     if (xcoords(2) < 0.00001)
		//         return - B1 * xcoords(1) * cos(eps + (π * t * sqrt(mu/rho))/l);
		//     else return 0.0; });

		// FunctionCoefficient dircy([t, &mu, &rho, &l, &B1, &eps](const Vector &xcoords)
		//                           {
		//     if (xcoords(2) < 0.00001)
		//         return B1 * xcoords(0) * cos(eps + (π * t * sqrt(mu/rho))/l);
		//     else return 0.0; });

		// PWConstCoefficient for Dirichlet BCs by boundary attribute.
		PWConstCoefficient dircx(dispbcx), dircy(dispbcy);
		GridFunction ux(fespacebc), uy(fespacebc);
		ux = uy = 0.0;

		ux.ProjectBdrCoefficient(dircx, ess_bdr_x);
		uy.ProjectBdrCoefficient(dircy, ess_bdr_y);

		// For loops to replace vector grid function degrees of freedom with presribed Dirichlet values
		for (int i = 0; i < ess_tdof_listbcx.Size(); i++)
		{
			int vdof = fespace->DofToVDof(ess_tdof_listbcx[i], 0);
			u_new(vdof) = ux(ess_tdof_listbcx[i]);
		}
		for (int i = 0; i < ess_tdof_listbcy.Size(); i++)
		{
			int vdof = fespace->DofToVDof(ess_tdof_listbcy[i], 1);
			u_new(vdof) = uy(ess_tdof_listbcy[i]);
		}

		u_old = u;
		u = u_new;

		// Compute  strain and stress

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
			// std::cout << "saving data" << std::endl;
			pointDisplacement << std::fixed << std::setprecision(10) << t << "\t" << u(selectNode) << "\n";
			pvdc.SetCycle(cycle); // Record time step number
			pvdc.SetTime(t);	  // Record simulation time
			pvdc.Save();
			cout << "cycle\t:" << cycle << "/" << brainSimProps.N << endl;
		}
		cycle++;
	}

	pointDisplacement.close();

	//   Free the used memory.
	cout << "Output saved. " << endl;

	delete fespace;
	delete fespacebc;
	delete fec;
	delete mesh;

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
void InitializeOldSolution(const GridFunction &u, const GridFunction &v, const SparseMatrix &M, const SparseMatrix &K, const Vector &F, double dt, GridFunction &u_old)
{
	MFEM_ASSERT(u.FESpace() == u_old.FESpace(), "u and u_old must use same FESpace");

	Vector Ku(u.Size()), rhs(u.Size());
	K.Mult(u, Ku);
	add(F, -1.0, Ku, rhs);
	// Ku *= -1;

	Vector a0(u.Size());
	a0 = 0.0;
	CGSolver cg;
	cg.SetOperator(M);
	cg.SetRelTol(1e-12);
	cg.SetMaxIter(200);
	cg.SetPrintLevel(0);
	GSSmoother M_prec(M);
	cg.SetPreconditioner(M_prec);
	cg.Mult(rhs, a0); // a0 = M^{-1} (F - K u)

	u_old = u;
	u_old.Add(-dt, v);
	u_old.Add(0.5 * dt * dt, a0); // ✅ correct sign
}
// Initializes u_old = u - dt * v + 0.5 - dt² * (M⁻¹K u)
void InitializeOldSolution(const GridFunction &u, const GridFunction &v, const SparseMatrix &M, const SparseMatrix &K, double dt, GridFunction &u_old)
{
	MFEM_ASSERT(u.FESpace() == u_old.FESpace(), "u and u_old must use same FESpace");

	Vector Ku(u.Size());
	K.Mult(u, Ku);

	Vector a0(u.Size());
	a0 = 0.0;
	CGSolver cg;
	cg.SetOperator(M);
	cg.SetRelTol(1e-12);
	cg.SetMaxIter(200);
	cg.SetPrintLevel(0);
	GSSmoother M_prec(M);
	cg.SetPreconditioner(M_prec);
	cg.Mult(Ku, a0); // a0 = M^{-1}(-K u)

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
	std::string json_path = "./input/InputParameters/WaveDynamics2D.json"; // Change JSON file path accordingly.

	std::ifstream infile(json_path);
	if (!infile.is_open())
	{
		std::cerr << "Error: Could not open " << json_path << std::endl;
		return 1;
	}

	infile >> inputParameters;

	return 0;
}

int saveResults(json &inputParameters, ParaViewDataCollection &pvdc, std::ofstream &pointDisplacement, const int order, mfem::GridFunction &u, mfem::GridFunction &ϵ, mfem::GridFunction &σ, mfem::GridFunction &nodes)
{
	std::string resultsFolder = "./results/" + inputParameters["testName"].get<std::string>();
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

	std::string resultsFolder = "./results/" + inputParameters["testName"].get<std::string>();
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
		u_data << "Number of elements: " << Ω->GetNE() << "\n";
	}

	u_data.close();
	return {};
};

int determineDirichletDof(const mfem::Mesh &mesh,
						  const mfem::FiniteElementSpace &fespace,
						  mfem::Array<int> &ess_tdof_listx,
						  mfem::Array<int> &ess_tdof_listy)
{

	mfem::Array<int> all_bdr_attr_selector(mesh.bdr_attributes.Max());
	all_bdr_attr_selector = 1;

	mfem::Array<int> bdr_nodes_list;
	fespace.GetEssentialTrueDofs(all_bdr_attr_selector, bdr_nodes_list);
	// const mfem::GridFunction *nodes = mesh.GetNodes();
	// MFEM_VERIFY(nodes, "Mesh has no nodes. Did you call SetCurvature or SetNodalFESpace?");
	mfem::Vector c{2.5, 2.5};
	// for (int i = 0; i < ndofs; ++i)
	// {
	// 	std::cout << i << "\n";
	// 	mfem::Vector pt(2);
	// 	pt[0] = (*nodes)(i);
	// 	pt[1] = (*nodes)(i + ndofs);
	// 	// pt -= c;

	// 	if (pt[0] > 4.5 && pt[0] < 5.5 && pt[1] > 1.0)
	// 	{
	// 		std::cout << pt[0] << " " << pt[1] << std::endl;
	// 		ess_tdof_listx.Append(i);
	// 		ess_tdof_listy.Append(i);
	// 	}
	// }
	return 0;
}

int createNodeGridFunction(mfem::GridFunction &sel_nodes, const mfem::Array<int> &ess_tdof_list)
{
	// const GridFunction *nodes_ptr = mesh->GetNodes();

	// if (!nodes_ptr)
	// {
	// 	std::cerr << "Error: mesh does not have nodes. Did you call mesh.SetCurvature(...)?" << std::endl;
	// 	return 1;
	// }

	// int num_nodes = nodes_ptr->Size(); // total # of entries
	// int vdim = nodes_ptr->VectorDim(); // dimension of coordinates
	// int ndofs = num_nodes / vdim;	   // number of actual nodes

	for (auto i : ess_tdof_list)
	{
		sel_nodes(i) = 1.0;
	}
	// for (int i = 0; i < ndofs; ++i)
	// {
	// 	out << i;
	// 	for (int d = 0; d < vdim; ++d)
	// 	{
	// 		nodes(i + d * ndofs) = (*nodes_ptr)(i + d * ndofs);
	// 	}
	// }

	return 0;
}

#include <cmath>
#include <limits>

// Jacobi rotation to compute eigenvalues of a symmetric matrix
