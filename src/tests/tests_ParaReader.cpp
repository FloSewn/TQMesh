/*
* This file is part of the TQMesh library.  
* This code was written by Florian Setzwein in 2022, 
* and is covered under the MIT License
* Refer to the accompanying documentation for details
* on usage and license.
*/
#include <TQMeshConfig.h>
#include "STLHeaders.h"
#include "CppUtils.h"
#include "tests.h"
 
namespace ParaReaderTests 
{
using namespace CppUtils;

/*
 * -> query(Block) seems to work, since correct start and end line indices are set
 * -> Thus, the problem must occur during the interior boundary query
 * -> It is probably related to the while() loop
 * -> query(Param) is missing an if clause to check if found starting-line-index
 *    is behind ending-line-index of its current block
 *    -> This should be handeled by get_query_data() !
 * -> maybe also add an assertion for its ending-line-index and the 
 *    block's starting-line-index, since this should never occur
*/

/*********************************************************************
*
*********************************************************************/
template <typename T> 
void print_parameter(ParaReader& reader, const std::string& name)
{
  auto para = reader.get_parameter<T>(name);

  if (  para.type() == ParaType::scalar 
     || para.type() == ParaType::vector )
  {
    auto key = para.start_key();

    size_t n_values = para.columns();

    std::stringstream ss;

    for ( size_t i = 0; i < n_values; ++i )
    {
      ss << para.get_value(i);

      if ( i < n_values - 1 )
        ss << ", ";
    }

    LOG(INFO) << key << " " << ss.str();
    LOG(INFO) << "";

  }
  else
  {
    auto start_key = para.start_key();
    auto end_key = para.end_key();

    size_t n_rows = para.rows();
    size_t n_cols = para.columns();

    LOG(INFO) << start_key;

    for ( size_t j = 0; j < n_rows; ++j )
    {
      std::stringstream ss;

      for ( size_t i = 0; i < n_cols; ++i )
      {
        ss << para.get_value(i,j);

        if ( i < n_cols - 1 )
          ss << ", ";
      }
      
      LOG(INFO) << ss.str();
    }

    LOG(INFO) << end_key;
    LOG(INFO) << "";
  }

} // print_parameter()



/*********************************************************************
*
*********************************************************************/
void read_parameters()
{
  std::string source_directory { TQMESH_SOURCE_DIR };
  std::string filepath { source_directory + "/auxiliary/test_data/ParaReaderTests.input.para" };

  ParaReader reader { filepath };

  // Init mesh blocks
  reader.new_block_parameter(
      "mesh_reader", "Define mesh:", "End mesh");

  // Create sub-reader for mesh properties
  ParaReader& mesh_reader = reader.get_block("mesh_reader");

  mesh_reader.new_matrix_parameter<double>(
      "extr_bdry_coords", "Define exterior boundary:", 
      "End exterior boundary", 5);

  mesh_reader.new_matrix_parameter<double>(
      "intr_bdry_coords", "Define interior boundary:", 
      "End interior boundary", 5);



  int mesh_id = 0; 

  while( reader.query( "mesh_reader" ) )
  {
    ParaReader& mesh_block = reader.get_block("mesh_reader");

    LOG(INFO) << "";
    LOG(INFO) << "============== " << "Create mesh " << mesh_id 
              << " ==============";


    // External boundary definition via direct edge coordinates 
    if ( mesh_block.query<double>("extr_bdry_coords", true, -1.0) )
    {
      auto para_extr_bdry 
        = mesh_block.get_parameter<double>("extr_bdry_coords");

      print_parameter<double>(mesh_block, "extr_bdry_coords");
    }


    // Query and initialize interior boundary coordinate definitions
    while( mesh_block.query<double>("intr_bdry_coords", true, -1.0) )
    {
      auto para_intr_bdry 
        = mesh_block.get_parameter<double>("intr_bdry_coords");

      print_parameter<double>(mesh_block, "intr_bdry_coords");
    }


    ++mesh_id;
  }


} /* read_parameters() */

} // namespace ParaReaderTests

/*********************************************************************
* Run tests for: ParaReader.h
*********************************************************************/
void run_tests_ParaReader()
{
   
  // Reset debug logging ostream
  adjust_logging_output_stream("COUT");

  ParaReaderTests::read_parameters();

} // run_tests_ParaReader()
