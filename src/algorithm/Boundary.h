/*
* This source file is part of the tqmesh library.  
* This code was written by Florian Setzwein in 2022, 
* and is covered under the MIT License
* Refer to the accompanying documentation for details
* on usage and license.
*/
#pragma once

#include <TQMeshConfig.h>
#include "STLHeaders.h"
#include "CppUtils.h"

#include "Error.h"
#include "EdgeList.h"
#include "Edge.h"
#include "Vertex.h"

namespace TQMesh {

using namespace CppUtils;

/********************************************************************* 
* Boundary type orientation
*********************************************************************/
enum class BdryType
{
  EXTERIOR,
  INTERIOR 
}; 

class Domain;


/********************************************************************* 
* A simple class to read mesh boundary defining CSV files 
*********************************************************************/
class CSVBoundaryReader
{
public:
  /*------------------------------------------------------------------
  | Constructor
  ------------------------------------------------------------------*/
  CSVBoundaryReader(const std::string& filename) 
  { read_csv_file(filename); }

  /*------------------------------------------------------------------
  | Getter
  ------------------------------------------------------------------*/
  std::vector<Vec2d>& vertex_coords() { return vertex_coords_; }
  std::vector<Vec2d>& vertex_props() { return vertex_props_; }
  std::vector<int>& edge_colors() { return edge_colors_; }
  const std::vector<Vec2d>& vertex_coords() const { return vertex_coords_; }
  const std::vector<Vec2d>& vertex_props() const { return vertex_props_; }
  const std::vector<int>& edge_colors() const { return edge_colors_; }

private:
  /*------------------------------------------------------------------
  | Attributes
  ------------------------------------------------------------------*/
  void read_csv_file(const std::string& filename)
  {
    vertex_coords_.clear();
    vertex_props_.clear();
    edge_colors_.clear();

    std::ifstream file { filename };

    if ( !file.is_open() )
      throw_error("Failed to open CSV file " + filename);

    std::string line;

    int i_line = 1;

    while ( std::getline(file, line) )
    {
      Vec2d coords { 0.0, 0.0 };
      Vec2d props { 0.0, 0.0 };
      int color { -1 };

      std::istringstream line_stream( line );
      std::string cell;

      int i_cell = 0;

      while ( std::getline(line_stream, cell, ',') )
      {
        try 
        {
          if ( i_cell == 0 )
            coords.x = std::stod(cell);
          if ( i_cell == 1 )
            coords.y = std::stod(cell);
          if ( i_cell == 2 )
            color = std::stoi(cell);
          if ( i_cell == 3 )
            props.x = std::stod(cell);
          if ( i_cell == 4 )
            props.y = std::stod(cell);
        }
        catch (const std::invalid_argument& e)
        {
          throw_error("Failed to read CSV file " + filename 
              + ": Invalid value in line " + std::to_string(i_line));
        }
        catch (const std::out_of_range& e)
        {
          throw_error("Failed to read CSV file " + filename 
              + ": Invalid value in line " + std::to_string(i_line));
        }

        ++i_cell;
      }

      if ( i_cell < 2 )
        throw_error("Failed to read CSV file " + filename 
            + ": Missing input data in in line " + std::to_string(i_line));

      if ( color < 0 )
        throw_error("Failed to read CSV file " + filename 
            + ": Invalid boundary color definition in line " + std::to_string(i_line));

      vertex_coords_.push_back( coords );
      vertex_props_.push_back( props );
      edge_colors_.push_back( color );
      ++i_line;
    }

  } // CSVBoundaryReader::read_csv_file()

  /*------------------------------------------------------------------
  | Attributes
  ------------------------------------------------------------------*/
  std::vector<Vec2d> vertex_coords_;
  std::vector<Vec2d> vertex_props_;
  std::vector<int>   edge_colors_;

}; // CSVBoundaryReader



/*********************************************************************
* A mesh boundary - defined by a list of edges
* > Interior boundaries are defined clockwise (CW)
* > Exterior boundaries are defined counter-clockwise (CCW)
*********************************************************************/
class Boundary : public EdgeList
{
public:

  /*------------------------------------------------------------------
  | Constructor
  ------------------------------------------------------------------*/
  Boundary(Vertices& domain_vertices, BdryType btype)
  : EdgeList( (btype == BdryType::EXTERIOR 
              ? Orientation::CCW : Orientation::CW) ) 
  , domain_vertices_ { &domain_vertices }
  , btype_ { btype }
  { }

  /*------------------------------------------------------------------
  | Getters
  ------------------------------------------------------------------*/
  bool is_exterior() const 
  { return (btype_ == BdryType::EXTERIOR); }
  bool is_interior() const
  { return (btype_ == BdryType::INTERIOR); }

  /*------------------------------------------------------------------
  | Override insert_edge() method of parent EdgeList, since all 
  | boundary edges must be defined with an appropriate boundary 
  | color (integers > 0)
  ------------------------------------------------------------------*/
  Edge& insert_edge(const_iterator pos, Vertex& v1, Vertex& v2, 
                    int color)
  { 
    ASSERT( (color >= 0), 
            "Boundary colors must be greater than zero");

    Edge& new_edge = EdgeList::insert_edge(pos, v1, v2, color); 

    return new_edge;
  }

  /*------------------------------------------------------------------
  | Override add_edge() method of parent EdgeList, since all 
  | boundary edges must be defined with an appropriate boundary 
  | color (integers > 0)
  ------------------------------------------------------------------*/
  Edge& add_edge(Vertex& v1, Vertex& v2, int color)
  { return this->insert_edge( edges_.end(), v1, v2, color ); }

  /*------------------------------------------------------------------
  | Set the boundary to a shape that is defined from a list of 
  | connected vertices
  ------------------------------------------------------------------*/
  void set_shape_from_coordinates(const std::vector<Vec2d>& v_coords,
                                  const std::vector<int>& colors, 
                                  const std::vector<Vec2d>& v_props)
  { create_boundary_shape(v_coords, v_props, colors); }

  void set_shape_from_coordinates(const std::vector<Vec2d>& v_coords,
                                  const std::vector<int>& colors) 
  { 
    std::vector<Vec2d> v_props ( v_coords.size(), {-1.0, -1.0} );
    create_boundary_shape(v_coords, v_props, colors); 
  }

  /*------------------------------------------------------------------
  | Set the boundary to a shape that is defined from a list of 
  | connected vertices, defined in a csv file
  ------------------------------------------------------------------*/
  void set_shape_from_csv(const std::string& filename)
  {
    CSVBoundaryReader reader { filename };

    create_boundary_shape(reader.vertex_coords(), 
                          reader.vertex_props(), 
                          reader.edge_colors());

  } // Boundary::set_shape_from_csv()

  /*------------------------------------------------------------------
  | Set the boundary to square shape
  ------------------------------------------------------------------*/
  void set_shape_square(int color, 
                        const Vec2d& xy, double w,
                        double mesh_size=0.0, double mesh_range=0.0)
  {
    set_shape_rectangle(color, xy, w, w, mesh_size, mesh_range);

  } // Boundary::set_shape_square()

  /*------------------------------------------------------------------
  | Set the boundary to equilateral triangle
  ------------------------------------------------------------------*/
  void set_shape_triangle(int color,
                          const Vec2d& xy, double a,
                          double mesh_size=0.0, double mesh_range=0.0)
  {
    set_shape_circle(color, xy, a/sqrt(3), 3, mesh_size, mesh_range);

  } // Boundary::set_shape_square()

  /*------------------------------------------------------------------
  | Set the boundary to rectangular shape
  ------------------------------------------------------------------*/
  void set_shape_rectangle(int color,
                           const Vec2d& xy, double w, double h,
                           double mesh_size=0.0, double mesh_range=0.0)
  {
    const double half_w = 0.5 * w;
    const double half_h = 0.5 * h;

    std::vector<Vec2d> v_shape {
      {xy[0] - half_w, xy[1] - half_h},
      {xy[0] + half_w, xy[1] - half_h},
      {xy[0] + half_w, xy[1] + half_h},
      {xy[0] - half_w, xy[1] + half_h},
    };

    std::vector<Vec2d> v_properties { 
      {mesh_size, mesh_range},
      {mesh_size, mesh_range},
      {mesh_size, mesh_range},
      {mesh_size, mesh_range},
    };

    std::vector<int> e_colors(4, color);

    create_boundary_shape(v_shape, v_properties, e_colors);

  } // Boundary::set_shape_rectangle()

  /*------------------------------------------------------------------
  | Set the boundary to circular shape
  ------------------------------------------------------------------*/
  void set_shape_circle(int color,
                        const Vec2d& xy, double r, size_t n=30,
                        double mesh_size=0.0, double mesh_range=0.0)
  {
    if ( n < 3 ) 
      return;

    std::vector<Vec2d> v_shape;
    std::vector<Vec2d> v_properties;
    std::vector<int>   e_colors;

    double delta_a = 2.0 * M_PI / static_cast<double>(n);

    for ( size_t i = 0; i < n; ++i )
    {
      double ang_i = i * delta_a;
      Vec2d d_xy = { cos(ang_i), sin(ang_i) }; 

      v_shape.push_back( xy + r * d_xy );
      v_properties.push_back( {mesh_size, mesh_range} );
      e_colors.push_back( color );
    }

    create_boundary_shape(v_shape, v_properties, e_colors);

  } // Boundary::set_shape_Circle()


private:

  /*------------------------------------------------------------------
  | Create the boundary from a given shape
  ------------------------------------------------------------------*/
  void create_boundary_shape(const std::vector<Vec2d>& v_shape,
                             const std::vector<Vec2d>& v_properties,
                             const std::vector<int>& e_colors)
  {
    // Do nothinng if the boundary contains no edges
    if ( edges_.size() > 0 )
      return;

    double area = polygon_area(v_shape);
    bool is_ccw = area > 0.0;

    // Create new vertices
    std::vector<Vertex*> new_verts {};

    int N = static_cast<int>(v_shape.size());

    for (int i = 0; i < N; ++i)
    {
      const Vec2d& xy = v_shape[i];
      const Vec2d& props = v_properties[i];

      // Search for already exiting vertex at this location
      Vertex* nearest = domain_vertices_->get_nearest(xy);

      if ( nearest && EQ0( ( nearest->xy() - xy ).norm_sqr() ) ) 
      {
        new_verts.push_back( nearest );
        continue;
      }

      // No nearest vertex found -> Generate new one
      Vertex& v_new = domain_vertices_->push_back(xy, props.x, props.y);
      new_verts.push_back( &v_new );
    }

    // CCW orientation for exteriror boundary
    // CW orientation for interior boundary
    if ( (btype_ == BdryType::EXTERIOR &&  is_ccw) ||
         (btype_ == BdryType::INTERIOR && !is_ccw) )
    {
      for (int i = 0; i < N; ++i)
      {
        int j = MOD(i+1, N);
        add_edge( *new_verts[i], *new_verts[j], e_colors[i] );
      }
    }
    else 
    {
      for (int i = N-1; i >= 0; --i)
      {
        int j = MOD(i+1, N);
        add_edge( *new_verts[j], *new_verts[i], e_colors[j] );
      }
    }
  };

  /*------------------------------------------------------------------
  | Boundary attributes
  ------------------------------------------------------------------*/
  BdryType  btype_;
  Vertices* domain_vertices_;


}; // Boundary

} // namespace TQMesh
