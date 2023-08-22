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
void simple_triangular_mesh();
void square_in_channel();
void boundary_shapes();
void fixed_vertices();
void merge_meshes();
void tqmesh_banner();
void run_example_7();
void run_example_8();
void run_example_9();
