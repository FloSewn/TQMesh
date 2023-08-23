/*
* This source file is part of the tqmesh library.  
* This code was written by Florian Setzwein in 2022, 
* and is covered under the MIT License
* Refer to the accompanying documentation for details
* on usage and license.
*/
#pragma once

#include "VecND.h"
#include "VtkIO.h"

#include "Mesh.h"
#include "MeshCleanup.h"

namespace TQMesh {
namespace TQAlgorithm {

using namespace CppUtils;


/*********************************************************************
* 
*********************************************************************/
enum class MeshExportType { 
  COUT, 
  TXT, 
  VTU 
};


/*********************************************************************
* Class for the export of meshes
*********************************************************************/
class MeshWriter
{
public:
  /*------------------------------------------------------------------
  | Constructor
  ------------------------------------------------------------------*/
  MeshWriter(Mesh& mesh, const Domain& domain)
  : mesh_ { &mesh }
  , domain_ { &domain }
  {}

  ~MeshWriter() {}

  /*------------------------------------------------------------------
  | 
  ------------------------------------------------------------------*/
  bool write(const std::string& filename, MeshExportType export_type)
  {
    MeshCleanup::assign_size_function_to_vertices(*mesh_, *domain_);
    MeshCleanup::assign_mesh_indices(*mesh_);
    MeshCleanup::setup_facet_connectivity(*mesh_);

    switch (export_type)
    {
      case MeshExportType::COUT:
        std::cout << (*mesh_);
        return true;

      case MeshExportType::TXT:
        return write_to_txt( filename );

      case MeshExportType::VTU:
        return write_to_vtu( filename );
    }

    return false;
  }

private:

  /*------------------------------------------------------------------
  | Export the mesh to a text file
  ------------------------------------------------------------------*/
  bool write_to_txt(const std::string& filepath)
  {
    std::string fullpath = filepath;

    if(fullpath.substr(fullpath.find_last_of(".") + 1) != "txt") 
      fullpath += ".txt";

    std::ofstream outfile;

    outfile.open( fullpath );

    outfile << (*mesh_);

    outfile.close();

    return true;

  } // MeshWriter::write_to_txt()

  /*------------------------------------------------------------------
  | Export the mesh to a vtu file
  ------------------------------------------------------------------*/
  bool write_to_vtu(const std::string& filepath)
  {
    std::string fullpath = filepath;

    if(fullpath.substr(fullpath.find_last_of(".") + 1) != "vtu") 
      fullpath += ".vtu";

    // Create data structure for VTU format
    std::vector<double> points {};
    std::vector<size_t> connectivity {};
    std::vector<size_t> offsets {};
    std::vector<size_t> types {};
    std::vector<double> size_function {};
    std::vector<int>    in_quad_layer {};
    std::vector<int>    is_fixed {};
    std::vector<int>    element_color {};
    std::vector<double> edge_length {};
    std::vector<double> max_angle {};
    std::vector<double> cell_quality {};

    size_t i_offset = 0;

    for ( const auto& v_ptr : mesh_->vertices() )
    {
      points.push_back( v_ptr->xy().x );
      points.push_back( v_ptr->xy().y );
      points.push_back( 0.0 );

      size_function.push_back( domain_->size_function(v_ptr->xy()) ); 
      in_quad_layer.push_back( static_cast<int>( v_ptr->in_quad_layer() ) );
      is_fixed.push_back( static_cast<int>( v_ptr->is_fixed() ) );
    }

    for ( const auto& q_ptr : mesh_->quads() )
    {
      connectivity.push_back( q_ptr->v1().index() );
      connectivity.push_back( q_ptr->v2().index() );
      connectivity.push_back( q_ptr->v3().index() );
      connectivity.push_back( q_ptr->v4().index() );

      i_offset += 4;
      offsets.push_back( i_offset );

      /// Type == 9 -> VTK_QUAD
      types.push_back( 9 );

      element_color.push_back( q_ptr->color() );

      edge_length.push_back( q_ptr->max_edge_length() );
      max_angle.push_back( q_ptr->max_angle() * 180. / M_PI );

      const double s = domain_->size_function(q_ptr->xy());
      cell_quality.push_back( q_ptr->quality(s) );
    }

    for ( const auto& t_ptr : mesh_->triangles() )
    {
      connectivity.push_back( t_ptr->v1().index() );
      connectivity.push_back( t_ptr->v2().index() );
      connectivity.push_back( t_ptr->v3().index() );

      i_offset += 3;
      offsets.push_back( i_offset );

      /// Type == 5 -> VTK_TRIANGLE
      types.push_back( 5 );

      element_color.push_back( t_ptr->color() );

      edge_length.push_back( t_ptr->max_edge_length() );
      max_angle.push_back( t_ptr->max_angle() * 180. / M_PI );

      const double s = domain_->size_function(t_ptr->xy());
      cell_quality.push_back( t_ptr->quality(s) );
    }


    VtuWriter writer { points, connectivity, offsets, types };

    writer.add_point_data( size_function, "size_function", 1 );
    writer.add_point_data( in_quad_layer, "in_quad_layer", 1 );
    writer.add_point_data( is_fixed, "fixed_vertices", 1 );
    writer.add_cell_data( element_color, "element_color", 1 );
    writer.add_cell_data( edge_length, "edge_length", 1 );
    writer.add_cell_data( max_angle, "max_angle", 1 );
    writer.add_cell_data( cell_quality, "cell_quality", 1 );

    writer.write( fullpath );

    return true;

  } // MeshWriter::write_to_vtu()

  /*------------------------------------------------------------------
  | Attributes
  ------------------------------------------------------------------*/
  Mesh*         mesh_;
  const Domain* domain_;

}; // MeshWriter


} // namespace TQAlgorithm
} // namespace TQMesh
