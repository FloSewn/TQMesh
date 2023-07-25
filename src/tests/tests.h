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
* Adjust logging output stream
* -> Log debug messages to specified output-file
*********************************************************************/
void adjust_logging_output_stream(const std::string& file);

/*********************************************************************
* The main test function
*********************************************************************/
int run_tests(const std::string& library);

/*********************************************************************
* Test functions
*********************************************************************/
void run_tests_Vertex();
void run_tests_Triangle();
void run_tests_Quad();
void run_tests_Front();
void run_tests_EdgeList();
void run_tests_Boundary();
void run_tests_SizeFunction();
void run_tests_Smoother();
void run_tests_Mesh();
void run_tests_MeshGenerator();
void run_tests_Cleanup();
