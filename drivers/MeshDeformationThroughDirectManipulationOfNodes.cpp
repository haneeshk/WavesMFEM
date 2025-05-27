#include "mfem.hpp"
#include <iostream>
#include <fstream>
#include <cstdlib> // For system()

using namespace mfem;
using namespace std;

Vector CoordinatesOfDoF(const GridFunction &vfield, int i, std::string PreName = "X", std::string PostName = "Vec")
{
	int vdim = vfield.VectorDim();	  // Number of components (e.g., 2 for x/y)
	int ndofs = vfield.Size() / vdim; // Scalar DOFs per component
	MFEM_VERIFY(i >= 0 && i < ndofs, "Invalid DOF index");

	Vector coords(vdim); // To store [x, y] or [x, y, z], etc.
	for (int d = 0; d < vdim; ++d)
	{
		coords[d] = vfield(i + d * ndofs); // Component-wise layout
	}

	std::cout << PreName + std::to_string(i) + PostName + "={" << coords[0] << ", " << coords[1] << "}" << std::endl;
	// coords.Print(std::cout, 8); // Print nicely
	return coords;
}

Vector CoordinatesOfDoF(Mesh &mesh, const int i, std::string PreName = "X", std::string PostName = "Vec")
{
	GridFunction *nodes = mesh.GetNodes();
	int vdim = nodes->VectorDim();	  // Embedding dimension (e.g., 2 for x/y)
	int ndofs = nodes->Size() / vdim; // Number of FE nodes (scalar DOFs)
	MFEM_VERIFY(i >= 0 && i < ndofs, "Invalid DOF index");
	mfem::Vector coords(2); // Allocate and own size-2 array

	{

		double x = (*nodes)(i);			// x-coordinate
		double y = (*nodes)(i + ndofs); // y-coordinate (since coordinates are stored component-wise)
		coords[0] = x;
		coords[1] = y;

		std::cout << PreName + std::to_string(i) + PostName + "={" << x << ", " << y << "}" << std::endl;
	}
	return coords;
}

int main(int argc, char *argv[])
{

	int nx = 10, ny = 10;

	Mesh *mesh = new Mesh(
		Mesh::MakeCartesian2D(
			nx, ny,
			Element::QUADRILATERAL, // element type
			true,					// generate_nodes
			1.0, 1.0				// domain size in x and y
			));

	mesh->PrintInfo();

	int dim = mesh->Dimension();
	int order = 1;
	FiniteElementCollection *fec;
	FiniteElementSpace *fespace;
	fec = new H1_FECollection(order, dim);
	fespace = new FiniteElementSpace(mesh, fec, dim);

	mesh->SetNodalFESpace(fespace);

	GridFunction *nodes = mesh->GetNodes();
	MFEM_VERIFY(nodes, "Mesh must be created with generate_nodes = true");
	const int vdim = nodes->VectorDim(); // Should be 2 for 2D
	const int ndofs = nodes->FESpace()->GetNDofs();
	cout << "vdim" << vdim << endl;
	cout << "ndofs" << ndofs << endl;

	socketstream glvis("localhost", 19916);
	glvis << "mesh\n"
		  << *mesh << flush;
	ParaViewDataCollection pvdc("deformed_output", mesh);
	pvdc.SetPrefixPath("results");
	pvdc.SetLevelsOfDetail(order);
	pvdc.SetHighOrderOutput(true);

	// Step 2: Save original mesh at time = 0
	pvdc.SetTime(0.0);
	pvdc.SetCycle(0);

	GridFunction displacement(fespace);
	displacement = 0.0;								   // Initialize to zero
	pvdc.RegisterField("displacement", &displacement); // Even if displacement is 0
	pvdc.Save();

	double t = 1.0;
	pvdc.SetTime(t);
	pvdc.SetCycle(1); // or any step index you like

	Vector x(2); // assuming 2D

	for (int i = 0; i < ndofs; i++)
	{
		x[0] = (*nodes)(i);
		x[1] = (*nodes)(i + ndofs);
		(*nodes)(i) = x[0] + 0.3 * x[1];
		(*nodes)(i + ndofs) = x[1];
	}

	pvdc.Save();
	socketstream glvis2("localhost", 19916);
	glvis2 << "mesh\n"
		   << *mesh << flush;

	return 0;
}