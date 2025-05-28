
#include "save.hpp" // Always include your own header first

#include <iostream> // For std::cout, std::cerr, etc.
#include <fstream>	// For std::ofstream
#include <iomanip>	// For std::setw, std::setprecision
#include <sstream>	// For std::ostringstream
#include <cstdlib>	// For system()
#include <cmath>	// If you use sin, cos, etc.

using namespace mfem;
using namespace std;

void saveData(std::string folder, const int cycle, const Vector &u, const GridFunction &nodes2)
{

	std::ostringstream oss;
	oss << std::setw(5) << std::setfill('0') << cycle;
	std::string padded = oss.str();
	std::ofstream u_out(folder + "/u_data_t" + padded + ".dat");

	for (int i = 0; i < u.Size(); i++)
	{
		double x1 = nodes2(i);
		double val = u(i);
		u_out << x1 << " " << val << "\n";
	}
	u_out.close();

	return;
}

void PrintVectorFormatted(const mfem::Vector &v,
						  int entries_per_line = 3,
						  int precision = 8,
						  int starting_line_number = 0,
						  int line_increment = 1)
{
	std::ios old_state(nullptr);
	old_state.copyfmt(std::cout); // Save current formatting

	std::cout << std::fixed << std::setprecision(precision);

	int line_num = starting_line_number;
	int width = precision + 6;

	for (int i = 0; i < v.Size(); ++i)
	{
		if (i % entries_per_line == 0)
		{
			if (i > 0)
				std::cout << "\n";
			std::cout << std::setw(6) << line_num << ": ";
			line_num += line_increment;
		}

		std::cout << std::setw(width) << v[i];
	}
	std::cout << "\n";

	std::cout.copyfmt(old_state); // Restore formatting
}

void WriteDOFCoordinates(const Mesh &mesh, const std::string &filename)
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

void mesh_fespaceinfo(const Mesh &mesh, const FiniteElementSpace &fespace, const FiniteElementSpace &vfespace)
{
	const GridFunction *nodes = mesh.GetNodes();
	// GridFunction u_scalar(&fespace);
	// GridFunction u_vector(&vfespace);
	std::cout << "Press Enter to continue...";
	std::cin.get();

	cout << "Mesh Information:\n";
	cout << "------------------\n";
	// cout << "Mesh file: " << mesh_file << endl;
	cout << "Mesh dimension: mesh.Dimension() :" << mesh.Dimension() << endl;
	std::cin.get();
	cout << "Embedding dimension (space dim): mesh.SpaceDimension() :" << mesh.SpaceDimension() << endl;
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

	cout << "fespace.GetNDofs()\t" << fespace.GetNDofs() << endl;
	cout << "fespace.GetVDim()\t" << fespace.GetVDim() << endl;
	cout << "fespace.GetTrueVSize()\t" << fespace.GetTrueVSize() << endl;
	cout << endl;
	cout << "vfespace.GetNDofs()\t" << vfespace.GetNDofs() << endl;
	cout << "vfespace.GetVDim()\t" << vfespace.GetVDim() << endl;
	cout << "vfespace.GetTrueVSize()\t" << vfespace.GetTrueVSize() << endl;
	cout << endl;
	cout << endl;
	cout << "nodes->Size()\t" << nodes->Size() << endl;
	cout << "nodes->VectorDim()\t" << nodes->VectorDim() << endl;
	// cout << "u_scalar.Size()\t" << u_scalar.Size() << endl;
	// cout << "u_vector.Size()\t" << u_vector.Size() << endl;

	std::cout << "Press Enter to continue...";
	std::cin.get();

	return;
}

Vector GridFunctionValAtDOF(const GridFunction &vfield, int i, std::string PreName = "X", std::string PostName = "Vec")
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
