/*
* This file is part of the TQMesh library.  
* This code was written by Florian Setzwein in 2022, 
* and is covered under the MIT License
* Refer to the accompanying documentation for details
* on usage and license.
*/
#pragma once

#include <vector>
#include <string>

/*********************************************************************
* The main function to run the examples
*********************************************************************/
int run_examples(const std::string& library);

/*********************************************************************
* Test functions
*********************************************************************/
bool simple_triangular_mesh();
bool square_in_channel();
bool boundary_shapes();
bool fixed_vertices();
bool merge_meshes();
bool airfoil_from_csv();
bool multiple_meshes();
bool thin_fracture();
bool fixed_edges();
bool tqmesh_banner();
bool run_example_7();
bool run_example_8();
bool run_example_9();
