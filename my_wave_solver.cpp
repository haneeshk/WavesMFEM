#include "mfem.hpp"
#include <filesystem>  // C++17
#include <fstream>
#include <string>
#include <iomanip>  // at the top of your file if not already included

namespace fs = std::filesystem;
using namespace mfem;
using namespace std;



void saveData(std::string, const int, const Vector&, const GridFunction &);
int main()
{
   int N = 20;
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

	double nl=1.0;	
	double Am = 1.0;
	double ell =1.0;
	double c_o =1.0; 
	double omega_l= (nl+1./2)*M_PI *(c_o/ell);
   FunctionCoefficient u0_coeff([](const Vector &x) { return 0.0; });
   FunctionCoefficient v0_coeff([Am, omega_l](const Vector &x) { return Am * omega_l*cos(omega_l*x[0]); });
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

   BilinearForm k(&fespace2);
   k.AddDomainIntegrator(new DiffusionIntegrator);
   k.Assemble();
   k.Finalize();
   SparseMatrix K(k.SpMat());

   // Time integration setup
   double T=2*M_PI/omega_l;
   double t = 0.0, t_final = 5*T, dt = T/1000;
   GridFunction u_old(&fespace2);
   add(u,-dt,v,u_old);
   GridFunction u_new(u);
   Vector rhs(fespace2.GetVSize()); // Full-sized RHS
   ConstantCoefficient zero(0.0);

    

   std::string resultsFolder = "results/test3";
   fs::create_directories(resultsFolder); 
   ParaViewDataCollection pvdc("configFiles", &mesh2);
   pvdc.SetPrefixPath(resultsFolder);      // Directory to save data
   pvdc.SetLevelsOfDetail(order2);  // Optional: for visualization
   pvdc.SetHighOrderOutput(true);   // Keep high-order info
   pvdc.RegisterField("u", &u);    // Associate field with data collection
   
 
 // Time-stepping loop
 int cycle=1;
  
   for (t = 0.0; t <= t_final; t += dt)
   {
      // Apply Dirichlet BCs
      FunctionCoefficient g_t([omega_l,t](const Vector &x) {
         return  sin( omega_l*t); // g(t) at left, 0 at right
      });
      u.ProjectBdrCoefficient(g_t, left_bdr);
      u.ProjectBdrCoefficient(zero, right_bdr);

      // rhs = -K u
      K.Mult(u, rhs);
      rhs *= -1.0;

      // Solve M a = rhs
      GSSmoother M_inv(M);
      M_inv.Mult(rhs, a);

      // Central difference: u^{n+1} = 2u^n - u^{n-1} + dt^2 * a
      Vector temp(u.Size());
      add(u, u, temp);                 // 2u^n
      add(temp, -1.0, u_old, temp);   // 2u^n - u^{n-1}
      add(temp, dt*dt, a, u_new);     // u^{n+1}

      u_new.ProjectBdrCoefficient(g_t, left_bdr);
      u_new.ProjectBdrCoefficient(zero, right_bdr);

      u_old = u;
      u = u_new;
   

	


	  if ((cycle % 10) == 0)  {
		saveData(resultsFolder, cycle, u, nodes2);
		pvdc.SetCycle(cycle);   // Record time step number
        pvdc.SetTime(t);       // Record simulation time
        pvdc.Save();   
	  }
	cout<<"cycle\t:"<<cycle<<endl;
	cycle++;
	
	}

   // Output result
   ofstream u_data("u_data.txt");
   for (int i = 0; i < u.Size(); i++)
   {
      double x1 = nodes2(i);
      double val = u(i);
      u_data << x1 << " " << val << "\n";
   }
   u_data.close();
   cout<<"saved output to "+resultsFolder<<endl;
   cout<<"The  time period T is \t"<<T<<endl;
   cout<<"The  simulation time is t_final is \t"<<t_final<<endl;
   return 0;
}

void saveData(std::string folder, const int cycle, const Vector& u, const GridFunction &nodes2){
	
	std::ostringstream oss;
oss << std::setw(5) << 						std::setfill('0') << cycle;
	std::string padded = oss.str();  
	std::ofstream u_out(folder+"/u_data_t" + padded + ".dat");
		
		  for (int i = 0; i < u.Size(); i++)
		  {
			  double x1 = nodes2(i);
			  double val = u(i);
			  u_out << x1 << " " << val << "\n";
		  }
		  u_out.close();
	
	return;

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
