#include "mfem.hpp"
#include <iostream>
#include <fstream>
#include <cstdlib> // For system()

using namespace mfem;
using namespace std;

std::vector<double> linspace(double a, double b, int n)
{
	std::vector<double> result;
	if (n <= 0)
		return result;
	if (n == 1)
	{
		result.push_back(a);
		return result;
	}

	double step = (b - a) / (n - 1);
	for (int i = 0; i < n; ++i)
		result.push_back(a + i * step);

	return result;
}

void WriteNodeCoordinates(const Mesh &mesh, const std::string &filename)
{
	const GridFunction *nodes = mesh.GetNodes();

	if (!nodes)
	{
		std::cerr << "Error: mesh does not have nodes. Did you call mesh.SetCurvature(...)?" << std::endl;
		return;
	}

	int num_nodes = nodes->Size(); // total # of entries
	int vdim = nodes->VectorDim(); // dimension of coordinates
	int ndofs = num_nodes / vdim;  // number of actual nodes

	std::ofstream out(filename);
	if (!out)
	{
		std::cerr << "Error: could not open file '" << filename << "' for writing." << std::endl;
		return;
	}

	// out << "# NodeIndex";
	// for (int d = 0; d < vdim; d++)
	// {
	// 	out << " x" << d;
	// }
	// out << "\n";

	for (int i = 0; i < ndofs; ++i)
	{
		out << i;
		for (int d = 0; d < vdim; ++d)
		{
			out << " " << (*nodes)(i + d * ndofs);
		}
		out << "\n";
	}

	out.close();
	std::cout << "Wrote " << ndofs << " node coordinates to " << filename << std::endl;
}

bool plot_fx_data(Mesh &mesh, const GridFunction &u, string datafile = "temp/data.txt", bool plt = 0, double a = 0.0, double b = 1.0)
{

	ofstream out(datafile);
	out.precision(8);

	Vector x(1), ux(1);
	double x_min = 0.0, x_max = 1.0;

	Vector x_phys, u_val;
	for (int e = 0; e < mesh.GetNE(); e++)
	{
		ElementTransformation *T = mesh.GetElementTransformation(e);

		// Sample N points in the reference element [0,1]
		const int N = 20;
		for (int i = 0; i <= N; i++)
		{
			double xi = i / static_cast<double>(N); // reference coordinate in [0,1]
			IntegrationPoint ip;
			ip.x = xi;
			T->SetIntPoint(&ip);

			// Get physical location
			T->Transform(ip, x_phys);

			// Evaluate u at this point
			u.GetVectorValue(e, ip, u_val);

			out << x_phys(0) << " " << u_val(0) << "\n";
		}
	}
	out.close();
	if (plt)
	{
		system("gnuplot -persist -e \"plot 'temp/data.txt' with lines title 'u(x)'\"");
	}

	return 0;
}

int main(int argc, char *argv[])
{
	// 1. Parse command-line options

	int order = 3;
	Mesh mesh = Mesh::MakeCartesian1D(10, 1.0);
	Mesh mesh = Mesh::MakeCartesian1D(10, 1.0);
	mesh.SetCurvature(order, false); // Enable high-order geometry

	// Vector min, max;
	// mesh.GetBoundingBox(min, max);
	// cout
	// 	<< "x_min = " << min(0) << endl;
	// cout << "x_max = " << max(0) << endl;
	std::cout << "Press Enter to continue...";
	std::cin.get();

	{ // 3. Print general information
		cout << "Mesh Information:\n";
		cout << "------------------\n";
		// cout << "Mesh file: " << mesh_file << endl;
		cout << "Mesh dimension: mesh.Dimension()" << mesh.Dimension() << endl;
		std::cin.get();
		cout << "Embedding dimension (space dim): mesh.SpaceDimension() " << mesh.SpaceDimension() << endl;
		std::cin.get();
		cout << "Number of elements: mesh.GetNE()" << mesh.GetNE() << endl;
		std::cin.get();
		cout << "Number of boundary elements: mesh.GetNBE()" << mesh.GetNBE() << endl;
		std::cin.get();
		cout << "Number of vertices: mesh.GetNV() " << mesh.GetNV() << endl;
		std::cin.get();
		cout << "Number of edges: mesh.GetNEdges()" << mesh.GetNEdges() << endl;
		std::cin.get();
		cout << "Number of faces: mesh.GetNFaces()" << mesh.GetNFaces() << endl;

		// cout << "Number of boundary attributes: ";
		//  << mesh.GetNumBdrGeometries(mesh.Dimension()-1) << endl;
	}

	H1_FECollection fec(order, mesh.Dimension());
	FiniteElementSpace fespace(&mesh, &fec);
	FiniteElementSpace vfespace(&mesh, &fec, 2);
	GridFunction u_scalar(&fespace);
	GridFunction u_vector(&vfespace);
	const GridFunction *nodes = mesh.GetNodes();

	cout << "fespace.GetNDofs()\t" << fespace.GetNDofs() << endl;
	cout << "vfespace.GetNDofs()\t" << vfespace.GetNDofs() << endl;
	cout << "fespace.GetVDim()\t" << fespace.GetVDim() << endl;
	cout << "fespace.GetTrueVSize()\t" << fespace.GetTrueVSize() << endl;
	cout << "vfespace.GetVDim()\t" << vfespace.GetVDim() << endl;
	cout << "vfespace.GetTrueVSize()\t" << vfespace.GetTrueVSize() << endl;
	cout << "nodes->GetNDofs()\t" << nodes->Size() << endl;
	cout << "nodes->VectorDim()\t" << nodes->VectorDim() << endl;
	cout << "u_scalar.Size()\t" << u_scalar.Size() << endl;
	cout << "u_vector.Size()\t" << u_vector.Size() << endl;

	std::cin.get();

	GridFunction u(&fespace);
	FunctionCoefficient coeff([](const Vector &x)
							  { return sin(2 * M_PI * x(0)); });
	u.ProjectCoefficient(coeff);

	for (int n = 0; n < fespace.GetNDofs(); n++)
	{
		GridFunction phi_n(&fespace);
		phi_n = 0.0;
		phi_n(n) = 0.1;
		plot_fx_data(mesh, phi_n, "temp/phi" + std::to_string(n) + ".dat");
	}

	// 4. Print element types
	cout << "\nElement Types:\n";
	for (int i = 0; i < mesh.GetNE(); i++)
	{
		Geometry::Type geom = mesh.GetElementBaseGeometry(i);
		cout << "  Element " << i << ": " << Geometry::Name[geom] << endl;
	}

	if (nodes && nodes->Size() > 0)
	{
		cout << "nodes->Size()\t" << nodes->Size() << endl;
		int vdim = nodes->VectorDim(); // = space dimension
		cout << "nodes->VectorDim()\t" << vdim << endl;
		cout << "fespace.GetNDofs()\t" << fespace.GetNDofs() << endl;
		for (int i = 0; i < fespace.GetNDofs(); i++)
		{
			std::cout << "DOF " << i << ": (";
			for (int d = 0; d < vdim; d++)
			{
				std::cout << (*nodes)(i + d * fespace.GetNDofs());
				if (d < vdim - 1)
					std::cout << ", ";
			}
			std::cout << ")" << std::endl;
		}
	}
	else
	{
		std::cout << "Mesh is not curved. No node coordinates are stored." << std::endl;
	}

	WriteNodeCoordinates(mesh, "temp/nodes.dat");

	return 0;
}

// // 5. Print boundary element types
// 	cout << "\nBoundary Element Types:\n";
// 	for (int i = 0; i < mesh.GetNBE(); i++)
// 	{
// 		Geometry::Type geom = mesh.GetBdrElementBaseGeometry(i);
// 		cout << "  Bdr Element " << i << ": " << Geometry::Name[geom] << endl;
// 	}