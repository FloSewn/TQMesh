/*
* This source file is part of the tqmesh library.  
* This code was written by Florian Setzwein in 2022, 
* and is covered under the MIT License
* Refer to the accompanying documentation for details
* on usage and license.
*/
#pragma once

#include <functional>
#include <cassert>
#include <float.h>
#include <limits.h>

#ifdef WIN32
#include <corecrt_math_defines.h>
#endif

#include "Helpers.h"

namespace TQMesh {
namespace TQAlgorithm {

/*********************************************************************
* This class is a container for all constants that are used 
* in the code for the mesh generation 
*********************************************************************/
class MeshingConstants
{
public:
  /*------------------------------------------------------------------
  | Default constructor
  ------------------------------------------------------------------*/
  MeshingConstants() {}

  /*------------------------------------------------------------------
  | Setters 
  ------------------------------------------------------------------*/
  void minimum_element_size(double s)  
  { minimum_element_size_ = s; }

  void minimum_element_scaling(double s)
  { minimum_element_scaling_ = s; }

  void quad_layer_angle(double s) 
  { quad_layer_angle_ = s; }

  void quad_layer_factor(double s) 
  { quad_layer_factor_ = s; }

  void quad_layer_range(double s) 
  { quad_layer_range_ = s; }

  void quad_range_factor(double s) 
  { quad_range_factor_ = s; }

  void mesh_range_factor(double s) 
  { mesh_range_factor_ = s; }

  void base_height_factor(double s) 
  { base_height_factor_ = s; }

  void base_vertex_factor(double s) 
  { base_vertex_factor_ = s; }

  void min_cell_quality(double s) 
  { min_cell_quality_ = s; }

  void max_cell_angle(double s) 
  { max_cell_angle_ = s; }

  void wide_search_factor(double s) 
  { wide_search_factor_ = s; }

  void edge_search_factor(double s) 
  { edge_search_factor_ = s; }

  void sort_triangulation_front(int s) 
  { sort_triangulation_front_ = s; }

  /*------------------------------------------------------------------
  | Getters 
  ------------------------------------------------------------------*/
  double minimum_element_size() const 
  { return minimum_element_size_; }

  double minimum_element_scaling() const
  { return minimum_element_scaling_; }

  double quad_layer_angle() const
  { return quad_layer_angle_; }

  double quad_layer_factor() const
  { return quad_layer_factor_; }

  double quad_layer_range() const
  { return quad_layer_range_; }

  double quad_range_factor() const
  { return quad_range_factor_; }

  double mesh_range_factor() const
  { return mesh_range_factor_; }

  double base_height_factor() const
  { return base_height_factor_; }

  double base_vertex_factor() const
  { return base_vertex_factor_; }

  double min_cell_quality() const
  { return min_cell_quality_; }
  
  double max_cell_angle() const
  { return max_cell_angle_; }

  double wide_search_factor() const
  { return wide_search_factor_; }

  double edge_search_factor() const
  { return edge_search_factor_; }

  int sort_triangulation_front() const
  { return sort_triangulation_front_; }


  int interior_edge_marker() const
  { return interior_edge_marker_; }

  int default_element_color() const 
  { return default_element_color_; }

  int default_mesh_id() const 
  { return default_mesh_id_; }


private:
  /*------------------------------------------------------------------
  | Adjustable attributes 
  ------------------------------------------------------------------*/
  double minimum_element_size_    = 0.001;
  double minimum_element_scaling_ = 0.001;

  double quad_layer_angle_        = 1.57079633; // = 1/2 pi
  double quad_layer_factor_       = 4.0;
  double quad_layer_range_        = 0.75;

  double quad_range_factor_       = 0.50;
  double mesh_range_factor_       = 1.0;
  double base_height_factor_      = 0.43; // ~ sqrt(3) / 4
  double base_vertex_factor_      = 1.00;

  double min_cell_quality_        = 0.0;
  double max_cell_angle_          = M_PI;

  double wide_search_factor_      = 10.0;

  // This value is used to enlarge the search radius for edges
  double edge_search_factor_      = 1.5;

  // This parameter controls, when the advancing front should be 
  // sorted in every iteration of the triangulation
  int sort_triangulation_front_   = -1;

  /*------------------------------------------------------------------
  | Fixed attributes 
  ------------------------------------------------------------------*/
  int interior_edge_marker_ = -1;
  int default_element_color_ = 0;
  int default_mesh_id_ = 0;

}; // MeshingConstants

inline MeshingConstants CONSTANTS;

/*********************************************************************
* CONSTANTS
*********************************************************************/
constexpr double TQ_SMALL  = 1.0E-13; //DBL_EPSILON;
constexpr double TQ_MAX    = DBL_MAX;
constexpr double TQ_MIN    = DBL_MIN;

/*********************************************************************
* 
*********************************************************************/
enum class ExportType { cout, txt, vtu };

} // namespace TQAlgorithm
} // namespace TQMesh 
