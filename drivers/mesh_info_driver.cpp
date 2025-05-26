#include "mfem.hpp"
#include <iostream>

using namespace mfem;
using namespace std;

int main(int argc, char *argv[])
{
	// 1. Parse command-line options
	const char *mesh_file = "../mfem/data/star.mesh";
	OptionsParser args(argc, argv);
	args.AddOption(&mesh_file, "-m", "--mesh",
				   "Mesh file to inspect.");
	args.Parse();
	if (!args.Good())
	{
		args.PrintUsage(cout);
		return 1;
	}
	args.PrintOptions(cout);

	// 2. Read the mesh
	Mesh mesh(mesh_file, 1, 1, true);

	// 3. Print general information
	cout << "Mesh Information:\n";
	cout << "------------------\n";
	cout << "Mesh file: " << mesh_file << endl;
	cout << "Mesh dimension: " << mesh.Dimension() << endl;
	cout << "Embedding dimension (space dim): " << mesh.SpaceDimension() << endl;
	cout << "Number of elements: " << mesh.GetNE() << endl;
	cout << "Number of boundary elements: " << mesh.GetNBE() << endl;
	cout << "Number of vertices: " << mesh.GetNV() << endl;
	cout << "Number of edges: " << mesh.GetNEdges() << endl;
	cout << "Number of faces: " << mesh.GetNFaces() << endl;
	cout << "Number of attributes (volume): " << mesh.GetNumGeometries(mesh.Dimension()) << endl;
	cout << "Number of boundary attributes: ";
	//  << mesh.GetNumBdrGeometries(mesh.Dimension()-1) << endl;

	// 4. Print element types
	cout << "\nElement Types:\n";
	for (int i = 0; i < mesh.GetNE(); i++)
	{
		Geometry::Type geom = mesh.GetElementBaseGeometry(i);
		cout << "  Element " << i << ": " << Geometry::Name[geom] << endl;
	}

	// 5. Print boundary element types
	cout << "\nBoundary Element Types:\n";
	for (int i = 0; i < mesh.GetNBE(); i++)
	{
		Geometry::Type geom = mesh.GetBdrElementBaseGeometry(i);
		cout << "  Bdr Element " << i << ": " << Geometry::Name[geom] << endl;
	}

	socketstream sol_sock("localhost", 19916, false);
	sol_sock.precision(8);
	sol_sock << "mesh\n"
			 << mesh << flush;

	ParaViewDataCollection pvdc("mesh_output", &mesh);
	pvdc.SetPrefixPath("temp"); // Output directory
	pvdc.SetLevelsOfDetail(1);	// For high-order viz
	pvdc.SetCycle(0);
	pvdc.SetTime(0.0);
	pvdc.Save(); // Will write mesh_output.pvd and .vtu/.vtm files

	return 0;
}