#include "mfem.hpp"
#include <filesystem>  // C++17
#include <fstream>
#include <string>
#include <iomanip>  // at the top of your file if not already included
#include <cmath>    // Include cmath for pow function

#include "save.hpp"

namespace fs = std::filesystem;
using namespace mfem;
using namespace std;


namespace standing {

	//   double omega_l= (nl+1./2)*M_PI *(c_o/ell);

}

namespace quess {  
	double l =1.0; //length of the string
	double n= 4;	
	 double c =1.0; //wave speed.
	 double lambda=(4*l)/(2*n+1);
	 double gamma=(2*M_PI)/lambda;
	 double A = 0.1; // Amplitude of the loading
	 double omega= gamma*c;
	 double epsilon=0.1;
 
 
	 auto f=[](double epsilon, double x){if (x <= epsilon) {
		 return 1.0 + (2.0 * x * x * x) / (epsilon * epsilon * epsilon)
					- (3.0 * x * x) / (epsilon * epsilon);
	 } else {
		 return 0.0;
	 }};
	 
 
	 auto mollifier = [](double epsilon, double x) -> double {
		 if (std::abs(x) >= epsilon) return 0.0;
		 double r = x / epsilon;
		 double peak_inv = std::exp(1.0);
		 return std::exp(-1.0 / (1.0 - r * r)) * peak_inv;
	 };
 
	 }
 


	 

FunctionCoefficient initialize_velocity(){
	
	FunctionCoefficient v0_coeff([](const Vector &x) { 
		using namespace quess;	
		return A * omega*mollifier(epsilon,x[0]); }
	);
	
	return v0_coeff;
}



// Initializes u_old = u - dt * v + 0.5 * dt² * (M⁻¹ K u)
void InitializeOldSolution(const GridFunction &u,
	const GridFunction &v,
	const SparseMatrix &M,
	const SparseMatrix &K,
	double dt,
	GridFunction &u_old)
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
cg.Mult(Ku, a0);  // a0 = M^{-1} K u

u_old = u;
u_old.Add(-dt, v);
u_old.Add(-0.5 * dt * dt, a0);  // ✅ correct sign
}


// struct ElementAndIntegrationPoint
// {
//     int element;
//     mfem::IntegrationPoint ip;

//     ElementAndIntegrationPoint() = default;

//     ElementAndIntegrationPoint(int elem, const mfem::IntegrationPoint &pt)
//         : element(elem), ip(pt) {}
// };



// ElementAndIntegrationPoint ComputeElementAndIntegrationPoint(mfem::GridFunction &u, double x_val=0.0){
// 	mfem::Mesh &mesh = *u.FESpace()->GetMesh();
//     const int NE = mesh.GetNE();

//     mfem::Vector x(1);
//     x(0) = x_val;

//     for (int i = 0; i < NE; i++)
//     {
//         mfem::ElementTransformation *T = mesh.GetElementTransformation(i);
//         mfem::IntegrationPoint ip;

//         if (T->TransformBack(x, ip))
//         {
//             return ElementAndIntegrationPoint(i,ip);  // ✅ Correct: returns double
//         }
//     }

//     throw std::runtime_error("x = 0.0 not found in any element");
// }

// double EvaluateAtX0(mfem::GridFunction &u, double x_val=0.0)
// {
//     mfem::Mesh &mesh = *u.FESpace()->GetMesh();
//     const int NE = mesh.GetNE();

//     mfem::Vector x(1);
//     x(0) = x_val;

//     for (int i = 0; i < NE; i++)
//     {
//         mfem::ElementTransformation *T = mesh.GetElementTransformation(i);
//         mfem::IntegrationPoint ip;
		

//         if (T->TransformBack(x, ip))
//         {
//             cout<<"i is \t :"<<i<<","<<endl;
// 			cout<<"x_val is \t :"<<x_val<<","<<endl;
// 			cout<<"ip.x is \t :"<<ip.x<<","<<endl;
// 			cout<<"u.GetValue(i,ip) is \t :"<<u.GetValue(i,ip)<<","<<endl;
// 			return u.GetValue(i,ip);  // ✅ Correct: returns double
//         }
//     }

//     throw std::runtime_error("x = 0.0 not found in any element");
// }

int main()
{
   int N = 100;
   int order2 = 3;
   Mesh mesh2 = Mesh::MakeCartesian1D(N, 1.0);
   mesh2.SetCurvature(order2, false); // Enable high-order geometry

   H1_FECollection fec2(order2, 1);
   FiniteElementSpace fespace2(&mesh2, &fec2);

   // Output mesh for verification
   ofstream mesh_ofs2("my_mesh2.mesh");
   mesh2.Print(mesh_ofs2);

   const GridFunction &nodes2 = *mesh2.GetNodes();
 	
   if (mesh2.GetNodes() == nullptr) { cout << "mesh2.GetNodes() = nullptr\n"; }

   cout << "mesh2 order: " << order2 << endl;
   cout << "mesh2 dimension: " << mesh2.Dimension() << endl;
   cout << "fespace2.GetTrueVSize(): " << fespace2.GetTrueVSize() << endl;
   cout << "mesh2.GetNV(): " << mesh2.GetNV() << endl;

   // Solution and velocity
   GridFunction u(&fespace2);      // displacement
   GridFunction v(&fespace2);      // velocity
   GridFunction a(&fespace2);      // acceleration

   // Initial conditions
 
 	
	
   FunctionCoefficient u0_coeff([](const Vector &x) { return 0.0; });
   FunctionCoefficient v0_coeff=initialize_velocity();




   

   u.ProjectCoefficient(u0_coeff);
   v.ProjectCoefficient(v0_coeff);

   


// Define essential boundary conditions
   Array<int> left_bdr(mesh2.bdr_attributes.Max());
   Array<int> right_bdr(mesh2.bdr_attributes.Max());
   left_bdr = 0;
   right_bdr = 0;
   left_bdr[0] = 1;   // attribute 1 (left)
   right_bdr[1] = 1;  // attribute 2 (right)

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

   // Time integration setup
   double T=2*M_PI/quess::omega;
   double t = 0.0, t_final = 5*(2*quess::l/quess::c), dt = T/10000;
   GridFunction u_old(&fespace2);
   InitializeOldSolution(u, v, M, K, dt, u_old);
 
 
   GridFunction u_new(u);
   Vector rhs(fespace2.GetVSize()); // Full-sized RHS
   ConstantCoefficient zero(0.0);

    

   std::string resultsFolder = "results/test9";
   fs::create_directories(resultsFolder); 
   ParaViewDataCollection pvdc("configFiles", &mesh2);
   pvdc.SetPrefixPath(resultsFolder);      // Directory to save data
   pvdc.SetLevelsOfDetail(order2);  // Optional: for visualization
   pvdc.SetHighOrderOutput(true);   // Keep high-order info
   pvdc.RegisterField("u", &u);    // Associate field with data collection
   
 

   ofstream pointDisplacement(resultsFolder+"/"+"pointDisplacement.dat");
   const int selectNode=282;
   


 // Time-stepping loop
 int cycle=1;
 
   
for (t = dt; t <= t_final; t += dt)
   {
      
		   
	
	
	// Apply Dirichlet BCs
      FunctionCoefficient g_t([t](const Vector &x) {
		using namespace quess;
         return  A* sin( omega*t); // g(t) at left, 0 at right
      });
    
      // rhs = -K u
      K.Mult(u, rhs);
      rhs *= -1.0;

      
	  cg.Mult(rhs, a);  // solve M a = rhs 
      // Central difference: u^{n+1} = 2u^n - u^{n-1} + dt^2 * a
      Vector temp(u.Size());
      add(u, u, temp);                 // 2u^n
      add(temp, -1.0, u_old, temp);   // 2u^n - u^{n-1}
      add(temp, dt*dt, a, u_new);     // u^{n+1}

      u_new.ProjectBdrCoefficient(g_t, left_bdr);
      u_new.ProjectBdrCoefficient(zero, right_bdr);

      u_old = u;
      u = u_new;
   

	//   ElementAndIntegrationPoint ixi= ComputeElementAndIntegrationPoint(u,0.0); 
	  double u_leftEnd = u(10);
	  double expected = quess::A * sin(quess::omega * t);
  
	  std::cout << "u(0, t) = " << u_leftEnd << ", expected = " << expected
			<< ", error = " << std::abs(u_leftEnd - expected) << std::endl;
	 


	  if ((cycle % 50) == 0)  {
		saveData(resultsFolder, cycle, u, nodes2);
		pointDisplacement << std::fixed << std::setprecision(10) << t << "\t" << u(selectNode) << "\n";
		pvdc.SetCycle(cycle);   // Record time step number
        pvdc.SetTime(t);       // Record simulation time
        pvdc.Save();   
	  }
	cout<<"cycle\t:"<<cycle<<endl;
	cycle++;
	
	}

	pointDisplacement.close();
	
	// Output result
   ofstream u_data(resultsFolder+"/"+"SimulationDetails.txt");
   u_data<<"saved output to "+resultsFolder<<endl;
   u_data<<"The  time period T is \t"<<T<<endl;
   u_data<<"The  simulation time is t_final is \t"<<t_final<<endl;
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



// for(auto i : boundary_dofs){
// 	cout<<i+1<<"\n";
// }



	// LinearForm b(&fespace2);
	// b.AddDomainIntegrator(new DomainLFIntegrator(quadratic_loading));
	// b.Assemble();

// // Step 3: Set up the linear form b(.) corresponding to the right-hand side
// ConstantCoefficient one(1.0);
// FunctionCoefficient quadratic_loading(
// 	[](const Vector &x){
// 		return (x[0]<0.75)? 0.0 : 1;
		
// 	});



// {double x_0_val; 
	
// 	x_0_val=nodes2(0);  // get coordinate of DOF 0
// 	x_0_val=nodes2(79);  // get coordinate of DOF 0
// 	 std::cout << "DOF 0 coordinate: " << x_0_val << std::endl;}
 