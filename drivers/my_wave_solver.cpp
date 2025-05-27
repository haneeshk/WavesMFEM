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

namespace fs = std::filesystem;
using namespace mfem;
using namespace std;
using json = nlohmann::json;

constexpr double π = M_PI;

// Standing wave
// double omega_l= (nl+1./2)*M_PI *(c_o/ell);

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

// Initializes u_old = u - dt * v + 0.5 * dt² * (M⁻¹ K u)
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
	cg.Mult(Ku, a0); // a0 = M^{-1} K u

	u_old = u;
	u_old.Add(-dt, v);
	u_old.Add(-0.5 * dt * dt, a0); // ✅ correct sign
}

std::function<FunctionCoefficient(double)> CreateLeftBC_t(const json &inputParameters)
{
	double A = inputParameters["Boundary Conditions"]["A"];
	double epsilon = inputParameters["Initial Conditions"]["epsilon"];
	double omega = inputParameters["Physical Parameters"]["omega"];

	std::function<FunctionCoefficient(double)> f =
		[A, omega](double t)
	{
		return FunctionCoefficient([A, omega, t](const Vector &x) -> double
								   { return A * sin(omega * t); });
	};

	return f;
}

void initializeQuess(json &inputParameters)
{

	double l = inputParameters["Physical Parameters"]["Bar Length"];
	double n = inputParameters["Boundary Conditions"]["n"];
	double c = inputParameters["Physical Parameters"]["Bar Wave Speed"];

	double lambda = (4 * l) / (2 * n + 1);
	double gamma = (2.0 * π) / lambda;

	inputParameters["Physical Parameters"]["lambda"] = lambda;
	inputParameters["Physical Parameters"]["gamma"] = gamma;
	inputParameters["Physical Parameters"]["omega"] = gamma * c;
	double omega = inputParameters["Physical Parameters"]["omega"];
	inputParameters["Physical Parameters"]["Time Frequency"] = omega / (2 * π);
	inputParameters["Physical Parameters"]["Time Period"] = (2 * π) / omega;

	return;
}

bool readInputParameters(int argc, char *argv[], json &inputParameters)
{
	std::string json_path = (argc > 1) ? argv[1] : "input.json";

	std::ifstream infile(json_path);
	if (!infile.is_open())
	{
		std::cerr << "Error: Could not open " << json_path << std::endl;
		return 1;
	}

	infile >> inputParameters;
	initializeQuess(inputParameters);
	return 0;
}

void CheckSolution(const double t, const GridFunction &u, const json inputParameters)
{ //   ElementAndIntegrationPoint ixi= ComputeElementAndIntegrationPoint(u,0.0);
	double A = inputParameters["Boundary Conditions"]["A"];
	double omega = inputParameters["Physical Parameters"]["omega"];
	double u_leftEnd = u(0);
	double expected = A * sin(omega * t);

	std::cout << "u(0, t) = " << u_leftEnd << ", expected = " << expected
			  << ", error = " << std::abs(u_leftEnd - expected) << std::endl;

	return;
}

int main(int argc, char *argv[])
{

	int order2 = 3;

	json inputParameters;
	readInputParameters(argc, argv, inputParameters);
	cout
		<< inputParameters.dump(2) << std::endl;

	Mesh mesh2 = Mesh::MakeCartesian1D(inputParameters["Simulation Parameters"]["Mesh Parameters"]["Mesh Size N"], 1.0);
	mesh2.SetCurvature(order2, false); // Enable high-order geometry

	H1_FECollection fec2(order2, mesh2.Dimension());
	FiniteElementSpace fespace2(&mesh2, &fec2);

	// Output mesh for verification
	ofstream mesh_ofs2("my_mesh2.mesh");
	mesh2.Print(mesh_ofs2);

	const GridFunction &nodes2 = *mesh2.GetNodes();

	if (mesh2.GetNodes() == nullptr)
	{
		cout << "mesh2.GetNodes() = nullptr\n";
	}

	cout << "mesh2 order: " << order2 << endl;
	cout << "mesh2 dimension: " << mesh2.Dimension() << endl;
	cout << "fespace2.GetTrueVSize(): " << fespace2.GetTrueVSize() << endl;
	cout << "mesh2.GetNV(): " << mesh2.GetNV() << endl;

	// Solution and velocity
	GridFunction u(&fespace2); // displacement
	GridFunction v(&fespace2); // velocity
	GridFunction a(&fespace2); // acceleration

	// Initial conditions

	FunctionCoefficient u0_coeff([](const Vector &x)
								 { return 0.0; });
	FunctionCoefficient v0_coeff = initialize_velocity(inputParameters);

	u.ProjectCoefficient(u0_coeff);
	v.ProjectCoefficient(v0_coeff);

	// Define essential boundary conditions
	Array<int> left_bdr(mesh2.bdr_attributes.Max());
	Array<int> right_bdr(mesh2.bdr_attributes.Max());
	left_bdr = 0;
	right_bdr = 0;
	left_bdr[0] = 1;  // attribute 1 (left)
	right_bdr[1] = 1; // attribute 2 (right)

	// Mass and stiffness matrices

	BilinearForm m(&fespace2);
	m.AddDomainIntegrator(new MassIntegrator);
	m.Assemble();
	m.Finalize();
	SparseMatrix M(m.SpMat());

	CGSolver cg;
	cg.SetOperator(M);
	cg.SetRelTol(1e-12);
	cg.SetAbsTol(1e-15);
	cg.SetMaxIter(200);
	cg.SetPrintLevel(0);

	GSSmoother prec(M);
	cg.SetPreconditioner(prec);

	BilinearForm k(&fespace2);
	k.AddDomainIntegrator(new DiffusionIntegrator);
	k.Assemble();
	k.Finalize();
	SparseMatrix K(k.SpMat());

	const double t_final = inputParameters["Simulation Parameters"]["Total Simulation Time"];
	double dt;
	{
		double T = inputParameters["Physical Parameters"]["Time Period"];
		double N = inputParameters["Simulation Parameters"]["Time Step Factor"];
		dt = T / N;
	}
	GridFunction u_old(&fespace2);
	InitializeOldSolution(u, v, M, K, dt, u_old);

	GridFunction u_new(u);
	Vector rhs(fespace2.GetVSize()); // Full-sized RHS
	ConstantCoefficient zero(0.0);

	std::string resultsFolder = "results/" + inputParameters["testName"].get<std::string>();
	fs::create_directories(resultsFolder);
	ParaViewDataCollection pvdc("configFiles", &mesh2);
	pvdc.SetPrefixPath(resultsFolder); // Directory to save data
	pvdc.SetLevelsOfDetail(order2);	   // Optional: for visualization
	pvdc.SetHighOrderOutput(true);	   // Keep high-order info
	pvdc.RegisterField("u", &u);	   // Associate field with data collection

	ofstream pointDisplacement(resultsFolder + "/" + "pointDisplacement.dat");
	const int selectNode = 282;

	// socketstream sol_sock("localhost", 19916, false);
	// sol_sock.precision(8);

	// Time-stepping loop
	int cycle = 1;
	double SimulationTimeSteps = t_final / dt;
	auto g_t_generator = CreateLeftBC_t(inputParameters);

	for (double t = dt; t <= t_final; t += dt)
	{

		// rhs = -K u
		K.Mult(u, rhs);
		rhs *= -1.0;

		cg.Mult(rhs, a); // solve M a = rhs
		// Central difference: u^{n+1} = 2u^n - u^{n-1} + dt^2 * a
		Vector temp(u.Size());
		add(u, u, temp);			  // 2u^n
		add(temp, -1.0, u_old, temp); // 2u^n - u^{n-1}
		add(temp, dt * dt, a, u_new); // u^{n+1}

		// Apply Dirichlet BCs
		FunctionCoefficient g_t = g_t_generator(t);
		u_new.ProjectBdrCoefficient(g_t, left_bdr);
		u_new.ProjectBdrCoefficient(zero, right_bdr);

		u_old = u;
		u = u_new;

		if ((cycle % 50) == 0)
		{
			// CheckSolution(t, u, inputParameters);
			// saveData(resultsFolder, cycle, u, nodes2);
			pointDisplacement << std::fixed << std::setprecision(10) << t << "\t" << u(selectNode) << "\n";
			pvdc.SetCycle(cycle); // Record time step number
			pvdc.SetTime(t);	  // Record simulation time
			pvdc.Save();
			cout << "cycle\t:" << cycle << "/" << SimulationTimeSteps << endl;
			// sol_sock << "solution\n"
			//  << mesh2 << u << flush;
		}

		cycle++;
	}

	pointDisplacement.close();

	// Output result
	ofstream u_data(resultsFolder + "/" + "SimulationDetails.txt");
	u_data << "saved output to " + resultsFolder << endl;
	u_data << "The  time period T is \t" << inputParameters["Physical Parameters"]["Bar Wave Speed"] << endl;
	u_data << "The  simulation time is t_final is \t" << t_final << endl;
	u_data.close();

	return 0;
}

// for (int i = 0; i < mesh.GetNV(); i++) // GetNV() gives the number of vertices
//        {
//            const double *coords = mesh.GetVertex(i);
//            cout << "Node " << i << ": (" << coords[0] << ", " << coords[1];
//            if (mesh.Dimension() == 3) // For 3D meshes
//                cout << ", " << coords[2];
//            cout << ")" << endl;
//        }

// Use f(x)
// double x1 = 0.1;
// std::cout << "f(" << x1 << ") = " << q.f(x1) << "\n";

// Use mollifier(x)
// double x2 = 0.05;
// std::cout << "mollifier(" << x2 << ") = " << q.mollifier(x2) << "\n";
// save_fx_data(q, "f_data.dat");