#include "save.hpp" // for saveData(), PrintVectorFormatted(), etc.
#include "mfem.hpp"
#include <iostream>
#include <fstream>
#include <cstdlib> // For system()

using namespace mfem;
using namespace std;

int main(int argc, char *argv[])
{

	int nx = 3, ny = 3;

	Mesh *mesh = new Mesh(
		Mesh::MakeCartesian2D(
			nx, ny,
			Element::QUADRILATERAL, // element type
			true,					// generate_nodes
			1.0, 1.0				// domain size in x and y
			));

	mesh->PrintInfo();

	int dim = mesh->Dimension();
	int order = 3;
	FiniteElementCollection *fec;
	FiniteElementSpace *fespace;
	fec = new H1_FECollection(order, dim);
	fespace = new FiniteElementSpace(mesh, fec, dim);
	mesh->SetNodalFESpace(fespace);

	mesh_fespaceinfo(*mesh, *fespace, *fespace);

	GridFunction *nodes = mesh->GetNodes();
	MFEM_VERIFY(nodes, "Mesh must be created with generate_nodes = true");
	int vdim = nodes->VectorDim(); // Should be 2 for 2D
	const int ndofs = nodes->FESpace()->GetNDofs();
	cout << "vdim" << vdim << endl;
	cout << "ndofs" << ndofs << endl;

	ParaViewDataCollection pvdc("deformed_output", mesh);
	pvdc.SetPrefixPath("results");
	pvdc.SetLevelsOfDetail(order);
	pvdc.SetHighOrderOutput(true);

	// Step 2: Save original mesh at time = 0
	pvdc.SetTime(0.0);
	pvdc.SetCycle(0);

	GridFunction displacement(fespace), deform(fespace);
	deform = 0.0;
	displacement = 0.0;								   // Initialize to zero
	pvdc.RegisterField("displacement", &displacement); // Even if displacement is 0
	pvdc.Save();

	VectorFunctionCoefficient ripple(dim, [](const Vector &x, Vector &v)
									 {
										 double r = sqrt(x[0]*x[0] + x[1]*x[1]);
										 double amp = 0.05 * sin(6 * M_PI * r);
										 v[0] = 0.2* x[1];
										 v[1] =  0.2 * x[0]; });

	deform.ProjectCoefficient(ripple); // Project coefficient onto FE space

	const int SelectDOF = 10;
	CoordinatesOfDoF(*mesh, SelectDOF, "X", "Vec");
	*nodes += deform; // 6. Deform the mesh: mesh nodes += displacement
	GridFunctionValAtDOF(displacement, SelectDOF, "Displacement", "Vec");
	CoordinatesOfDoF(*mesh, SelectDOF, "x", "Vec");

	PrintVectorFormatted(displacement, 1, 8, 0, 1);

	double t = 1.0;
	pvdc.SetTime(t);
	pvdc.SetCycle(1); // or any step index you like

	pvdc.Save();

	return 0;
}