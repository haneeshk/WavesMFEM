#include "save.hpp" // for saveData(), PrintVectorFormatted(), etc.
#include "mfem.hpp"
#include <iostream>
#include <fstream>
#include <cstdlib> // For system()
#include <sstream>
#include <iomanip>

using namespace mfem;
using namespace std;

void SampleScalarGridFunction(Mesh *mesh,
							  const GridFunction &u,
							  const int samples_per_dir,
							  const std::string &filename)
{
	const int dim = mesh->Dimension();
	MFEM_VERIFY(dim == 2, "This code assumes a 2D mesh.");
	MFEM_VERIFY(mesh->GetNodes(), "Mesh must have nodes for this transformation.");
	MFEM_VERIFY(u.FESpace()->GetVDim() == 1, "GridFunction must be scalar-valued.");

	std::ofstream out(filename);
	out << std::scientific;

	for (int el = 0; el < mesh->GetNE(); el++)
	{
		ElementTransformation *T = mesh->GetElementTransformation(el);
		Vector phys_x;
		phys_x.SetSize(T->GetSpaceDim());

		for (int i = 0; i < samples_per_dir; i++)
		{
			double xi0 = (1.0 * i) / (samples_per_dir - 1.0);
			for (double j = 0; j < samples_per_dir; j++)
			{

				double xi1 = (1.0 * j) / (samples_per_dir - 1);

				IntegrationPoint ip;
				ip.Set2(xi0, xi1);
				T->SetIntPoint(&ip); // ← required for FEM evaluations

				T->Transform(ip, phys_x);

				double val = u.GetValue(el, ip);
				// std::cout << xi0 << "," << xi1 << endl;
				// std::cout << phys_x(0) << "," << phys_x(1) << " :" << val << " ....";
				// std::cout << endl;

				out << phys_x(0) << " " << phys_x(1) << " " << val << "\n";
			}
		}
	}

	out.close();
	std::cout << "Wrote sampled values to " << filename << std::endl;
}

int main(int argc, char *argv[])
{

	int nx = 2, ny = 2;

	Mesh *mesh = new Mesh(
		Mesh::MakeCartesian2D(
			nx, ny,
			Element::QUADRILATERAL, // element type
			true,					// generate_nodes
			1.0, 1.0				// domain size in x and y
			));

	mesh->PrintInfo();

	int dim = mesh->Dimension();
	int order = 2;
	FiniteElementCollection *fec;
	FiniteElementSpace *fespace, *scalar_fespace;
	fec = new H1_FECollection(order, dim);
	fespace = new FiniteElementSpace(mesh, fec, dim);
	scalar_fespace = new FiniteElementSpace(mesh, fec, 1);
	mesh->SetNodalFESpace(fespace);

	GridFunction *nodes = mesh->GetNodes();
	MFEM_VERIFY(nodes, "Mesh must be created with generate_nodes = true");
	const int vdim = nodes->VectorDim(); // Should be 2 for 2D
	const int ndofs = nodes->FESpace()->GetNDofs();
	cout << "vdim" << vdim << endl;
	cout << "node->Size()" << nodes->Size() << endl;
	cout << "ndofs" << ndofs << endl;

	ParaViewDataCollection pvdc("square", mesh);
	pvdc.SetPrefixPath("temp");
	pvdc.SetLevelsOfDetail(order);
	pvdc.SetHighOrderOutput(true);

	// Step 2: Save original mesh at time = 0
	pvdc.SetTime(0.0);
	pvdc.SetCycle(0);
	pvdc.Save();

	GridFunction displacement(fespace), scalar_basis_function(scalar_fespace);
	displacement = 0.0; // Initialize to zero
	scalar_basis_function = 0.0;

	// pvdc.RegisterField("displacement", &displacement); // Even if displacement is 0
	pvdc.RegisterField("basisFunction", &scalar_basis_function); // Even if
	pvdc.Save();
	displacement.Size();

	for (int n = 0; n < ndofs; n++)
	{
		// displacement = 0.0; // Initialize to zero
		scalar_basis_function = 0.0;
		displacement(n) = 0.1;
		scalar_basis_function(n) = 1;
		// cout << scalar_basis_function;
		cout << "ψ_" << n << endl;

		double t = (double)n;
		pvdc.SetTime(t);
		pvdc.SetCycle(n); // or any step index you like
		pvdc.Save();
		std::ostringstream filename;
		filename << "temp/basisFunction2DSecondOrderAt_" << std::setw(5) << std::setfill('0') << n << "_.dat";
		SampleScalarGridFunction(mesh, scalar_basis_function, 200, filename.str());

		// plot_fx_data(mesh, phi_n, "temp/phi" + std::to_string(n) + ".dat");
	}
	WriteDOFCoordinates(*mesh, "temp/nodes.dat");
	mesh->SetCurvature(order, true);
	std::ofstream mesh_out("temp/square.mesh");
	mesh->Print(mesh_out);
	mesh_out.close();

	return 0;
}
