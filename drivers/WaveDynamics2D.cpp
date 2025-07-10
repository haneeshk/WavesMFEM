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

FunctionCoefficient initialize_velocity(const json &inputParameters);
void InitializeOldSolution(const GridFunction &u, const GridFunction &v, const SparseMatrix &M, const SparseMatrix &K, const Vector &F, double dt, GridFunction &u_old);
void InitializeOldSolution(const GridFunction &u, const GridFunction &v, const SparseMatrix &M, const SparseMatrix &K, double dt, GridFunction &u_old);
bool readInputParameters(int argc, char *argv[], json &inputParameters);
double dirichletamp(const json &inputParameters, const double &t);
double tractionamp(const json &inputParameters, const double &t);

int main(int argc, char *argv[])
{
	//   Parse command-line options.
	const char *mesh_file = "../input/meshes/Rectangle-quad.mesh"; // Change this line, include 2D mesh file.
	int order = 2;
	bool static_cond = false;
	bool visualization = 0;
	int ref_levels = 3;

	json inputParameters;
	readInputParameters(argc, argv, inputParameters);

	OptionsParser args(argc, argv);
	args.AddOption(&order, "-o", "--order",
				   "Finite element order (polynomial degree).");
	args.AddOption(&visualization, "-vis", "--visualization", "-no-vis",
				   "--no-visualization",
				   "Enable or disable GLVis visualization.");
	args.AddOption(&ref_levels, "-r", "--ref_levels",
				   "Number of uniform mesh refinements.");
	args.Parse();

	// Read the mesh from the given mesh file.
	Mesh *mesh = new Mesh(mesh_file, 1, 1);
	mesh->SetCurvature(order, false); // Enable high-order geometry

	cout << "Mesh created" << endl;

	int dim = mesh->Dimension(); // dim should equal 2.
	cout << dim << endl;

	for (int l = 0; l < ref_levels; l++)
		mesh->UniformRefinement();

	//------------------------------------------------------------------------------------------------------------------
	// Define finite element spaces.
	// Define a H1 finite element space for nodal displacements.
	FiniteElementCollection *fec;
	FiniteElementSpace *fespace, *fespacebc;

	fec = new H1_FECollection(order, dim);
	fespace = new FiniteElementSpace(mesh, fec, dim);

	// Define a finite element space to project bcs in the loop.
	fespacebc = new FiniteElementSpace(mesh, fec);

	cout << "Number of finite element unknowns: " << fespace->GetTrueVSize() << endl;

	// Define an L2 finite element space for element strain and stress.
	// Dimension of the L2 finite element space is the number of strain and stress components.
	// For 2D, we have 3 components.
	FiniteElementCollection *l2fec;
	FiniteElementSpace *l2fespace;

	l2fec = new L2_FECollection(0, dim);
	l2fespace = new FiniteElementSpace(mesh, l2fec, 3);
	//------------------------------------------------------------------------------------------------------------------

	//------------------------------------------------------------------------------------------------------------------
	//    Determine the list of essential boundary dofs for each vector dimension.

	Array<int> ess_tdof_listx, ess_tdof_listy,
		ess_bdr_x(mesh->bdr_attributes.Max()), ess_bdr_y(mesh->bdr_attributes.Max());

	ess_bdr_x = 0;
	ess_bdr_x[2] = 1; // right edge

	ess_bdr_y = 0;
	ess_bdr_y[2] = 1; // right edge

	fespace->GetEssentialTrueDofs(ess_bdr_x, ess_tdof_listx, 0);
	fespace->GetEssentialTrueDofs(ess_bdr_y, ess_tdof_listy, 1);

	Array<int> ess_tdof_list; // Arrange all true dofs in one array.

	ess_tdof_list.Append(ess_tdof_listx);
	ess_tdof_list.Append(ess_tdof_listy);

	// To project BCs in the loop, define scalar tdof lists for each component.
	Array<int> ess_tdof_listbcx, ess_tdof_listbcy;
	fespacebc->GetEssentialTrueDofs(ess_bdr_x, ess_tdof_listbcx, 0);
	fespacebc->GetEssentialTrueDofs(ess_bdr_y, ess_tdof_listbcy, 0);

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
	// initialize_velocity(inputParameters);
	u.ProjectCoefficient(u0_coeff);
	v.ProjectCoefficient(v0_coeff);
	//------------------------------------------------------------------------------------------------------------------

	//------------------------------------------------------------------------------------------------------------------
	// Read simulation parameters.
	const float lambda = inputParameters["Physical Parameters"]["lambda"];
	const float mu = inputParameters["Physical Parameters"]["mu"];
	ConstantCoefficient lambda_coeff(lambda), mu_coeff(mu);

	const float youngmod = inputParameters["Physical Parameters"]["youngmod"];
	const float poissonratio = inputParameters["Physical Parameters"]["poissonratio"];
	ConstantCoefficient E_coeff(youngmod), NU_coeff(poissonratio);

	const float rho = inputParameters["Physical Parameters"]["Density"];
	ConstantCoefficient density(rho);

	cout << "Physical parameters read" << endl;

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

	// Make sure all of the below parameters (or the ones that are necessary) are included in JSON file.
	const double t_final = inputParameters["Simulation Parameters"]["Total Simulation Time"];
	const double TimeSteps = inputParameters["Simulation Parameters"]["Time Steps"];
	const double dt = t_final / TimeSteps;

	GridFunction u_old(fespace);
	// InitializeOldSolution(u, v, M, K, dt, u_old); // comment if non-zero tractions are prescribed.
	InitializeOldSolution(u, v, M, K, F, dt, u_old); // un-comment if non-zero tractions are prescribed.

	GridFunction u_new(u);
	Vector rhs(fespace->GetVSize()); // Full-sized RHS
	ConstantCoefficient zero(0.0);

	// Initialize GlobalStressStrain object. Initialize stress and strain grid functions.
	mfemplus::GlobalStressStrain WaveDynamics(mesh, fespace);
	GridFunction eps(l2fespace), sig(l2fespace);
	eps = sig = 0.0;

	GridFunction eps_temp, sig_temp;
	WaveDynamics.GlobalStrain(u_new, eps_temp);
	WaveDynamics.GlobalStress(eps_temp, E_coeff, NU_coeff, sig_temp);

	eps = eps_temp;
	sig = sig_temp;
	// Inititialize Paraview object

	std::string resultsFolder = "../results/" + inputParameters["testName"].get<std::string>();
	fs::create_directories(resultsFolder);
	ParaViewDataCollection pvdc("Waves2D", mesh);

	pvdc.SetPrefixPath(resultsFolder);	// Directory to save data
	pvdc.SetLevelsOfDetail(order);		// Optional: for visualization
	pvdc.SetHighOrderOutput(true);		// Keep high-order info
	pvdc.RegisterField("u", &u);		// Associate displacement field with data collection
	pvdc.RegisterField("strain", &eps); // Associate strain field with data collection
	pvdc.RegisterField("stress", &sig); // Associate stress field with data collection

	ofstream pointDisplacement(resultsFolder + "/" + "pointDisplacement.dat");
	const int selectNode = 2;

	int cycle = 0;
	int t = 0;
	pvdc.SetCycle(cycle); // Record time step number
	pvdc.SetTime(t);	  // Record simulation time
	pvdc.Save();

	cycle += 1;

	// Time-stepping loop
	cout << "Time stepping loop beginning" << endl;

	for (double t = dt; t <= t_final; t += dt)
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
		// Central difference: u_{n+1} = 2u_n - u_{n-1} + dt^2 * a
		Vector temp(u.Size());
		add(u, u, temp);			  // 2u_n
		add(temp, -1.0, u_old, temp); // 2u_n - u_{n-1}
		add(temp, dt * dt, a, u_new); // u_{n+1}

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

		// Compute strain and stress
		GridFunction eps_temp, sig_temp;
		WaveDynamics.GlobalStrain(u_new, eps_temp);
		WaveDynamics.GlobalStress(eps_temp, E_coeff, NU_coeff, sig_temp);

		eps = eps_temp;
		sig = sig_temp;

		if ((cycle % 50) == 0)
		{
			pointDisplacement << std::fixed << std::setprecision(10) << t << "\t" << u(selectNode) << "\n";
			pvdc.SetCycle(cycle); // Record time step number
			pvdc.SetTime(t);	  // Record simulation time
			pvdc.Save();
			// cout << "cycle\t:" << cycle << "/" << TimeSteps << endl;
		}
		cycle++;
	}

	cout << "Loop completed." << endl;

	pointDisplacement.close();

	ofstream u_data(resultsFolder + "/" + "SimulationDetails.txt");
	u_data << "saved output to " + resultsFolder << endl;
	u_data << "The  time period T is \t" << inputParameters["Physical Parameters"]["Time Period"] << endl;
	u_data << "The  simulation time is t_final is \t" << t_final << endl;
	u_data.close();

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
