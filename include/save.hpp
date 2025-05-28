#pragma once
#include <string>
#include "mfem.hpp"

void saveData(std::string folder, const int cycle,
			  const mfem::Vector &u,
			  const mfem::GridFunction &nodes2);

// Nicely print an mfem::Vector with line numbers and formatting
void PrintVectorFormatted(const mfem::Vector &v,
						  int,
						  int,
						  int,
						  int);

// Write coordinates of all DOFs (nodes) to a text file
void WriteDOFCoordinates(const mfem::Mesh &mesh,
						 const std::string &filename);

// Sample a GridFunction on each element and save to file (optionally plot)
bool plot_fx_data(mfem::Mesh &mesh,
				  const mfem::GridFunction &u,
				  std::string datafile,
				  bool plt,
				  double a,
				  double b);

// Print detailed info about the mesh and FESpaces
void mesh_fespaceinfo(const mfem::Mesh &mesh,
					  const mfem::FiniteElementSpace &fespace,
					  const mfem::FiniteElementSpace &vfespace);

// Extract mfem::vector value of a GridFunction at a specific DOF index
mfem::Vector GridFunctionValAtDOF(const mfem::GridFunction &vfield,
								  int i,
								  std::string PreName,
								  std::string PostName);

// Get physical coordinates of the i-th DOF in the mesh
mfem::Vector CoordinatesOfDoF(mfem::Mesh &mesh,
							  const int i,
							  std::string PreName,
							  std::string PostName);