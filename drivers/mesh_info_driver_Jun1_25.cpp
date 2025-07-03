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
			double xi0 = i / (samples_per_dir - 1);
			for (int j = 0; j < samples_per_dir; j++)
			{

				double xi1 = j / (samples_per_dir - 1);

				IntegrationPoint ip;
				ip.Set2(xi0, xi1);
				T->SetIntPoint(&ip); // ← required for FEM evaluations

				T->Transform(ip, phys_x);

				double val = u.GetValue(el, ip);
				std::cout << xi0 << "," << xi1 << endl;
				std::cout << phys_x(0) << "," << phys_x(1) << " :" << val << " ....";
				std::cin.get();
				std::cout << endl;

				out << phys_x(0) << " " << phys_x(1) << " " << val << "\n";
			}
		}
	}

	out.close();
	std::cout << "Wrote sampled values to " << filename << std::endl;
}

int main(int argc, char *argv[])
{

	int nx = 1, ny = 1;

	Mesh *mesh = new Mesh(
		Mesh::MakeCartesian2D(
			nx, ny,
			Element::QUADRILATERAL, // element type
			true,					// generate_nodes
			1.2, 2.3				// domain size in x and y
			));
	Vector trans(2); // translation vector
	trans(0) = 1.0;	 // shift x by +1
	trans(1) = 2.0;	 // shift y by +2
	// mesh->MoveVertices(trans); // applies the translation to all vertices

	mesh->PrintInfo();

	int dim = mesh->Dimension();
	int order = 1;
	FiniteElementCollection *fec;
	FiniteElementSpace *fespace;
	fec = new H1_FECollection(order, dim);
	fespace = new FiniteElementSpace(mesh, fec, dim);
	// scalar_fespace = new FiniteElementSpace(mesh, fec, 1);
	mesh->SetNodalFESpace(fespace);

	mesh_fespaceinfo(*mesh, *fespace, *fespace);

	GridFunction *nodes = mesh->GetNodes();
	WriteDOFCoordinates(*mesh, "temp/nodes.txt");

	int el = 0; // pick the first element
	// MFEM_VERIFY(scalar_fespace->GetFE(0)->GetDof() ==
	// 				mesh->GetNodalFESpace()->GetFE(0)->GetDof(),
	// 			"Mismatch in FE definitions");
	const FiniteElement *fe = fespace->GetFE(el);
	std::cout << "Ref geometry: " << Geometry::Name[fe->GetGeomType()] << std::endl;
	std::cout << "fec->Name()" << fec->Name() << endl;

	std::cout << "DOFs on element " << el << ": " << fe->GetDof() << std::endl;

	IntegrationPoint ip;
	Vector shape(fe->GetDof());

	std::vector<std::pair<const mfem::real_t, const mfem::real_t>> test_points = {
		{-1.0, -1.0}, {1.0, -1.0}, {1.0, 1.0}, {-1.0, 1.0}, {0.0, 0.0}};

	for (auto [xi, eta] : test_points)
	{
		shape = 0.0;
		ip.Set2(xi, eta);
		fe->CalcShape(ip, shape);
		std::cout << "At (" << xi << ", " << eta << "): ";
		for (int i = 0; i < shape.Size(); i++)
			std::cout << shape(i) << " ";
		std::cout << "\n";
	}

	return 0;
}
