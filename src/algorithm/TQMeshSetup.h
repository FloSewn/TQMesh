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

namespace TQMesh {

/*********************************************************************
* 
*********************************************************************/
class TQMeshSetup {
  public:

  /*------------------------------------------------------------------
  | Constexprs 
  ------------------------------------------------------------------*/
  static constexpr int    interior_edge_color   { -1 };
  static constexpr int    default_element_color {  0 };
  static constexpr int    default_mesh_id       {  0 };
  static constexpr double dbl_small           { 1.0E-13 };
  static constexpr double dbl_min             { DBL_MIN };
  static constexpr double dbl_max             { DBL_MAX };

  /*------------------------------------------------------------------
  | Enforce singleton
  ------------------------------------------------------------------*/
  TQMeshSetup(const TQMeshSetup&) = delete;
  TQMeshSetup& operator=(const TQMeshSetup&) = delete;

  /*------------------------------------------------------------------
  | Get singleton instance 
  ------------------------------------------------------------------*/
  static TQMeshSetup& get_instance(void) 
  {
    static TQMeshSetup instance;
    return instance;
  }

  /*------------------------------------------------------------------
  | Setters
  ------------------------------------------------------------------*/
  TQMeshSetup& set_quadtree_scale(double s) 
  { quadtree_scale_ = s; return *this; }
  TQMeshSetup& set_quadtree_max_items(size_t n) 
  { quadtree_max_items_ = n; return *this; }
  TQMeshSetup& set_quadtree_max_depth(size_t n)
  { quadtree_max_depth_ = n; return *this; }

  /*------------------------------------------------------------------
  | Getters
  ------------------------------------------------------------------*/
  double get_quadtree_scale(void) const { return quadtree_scale_; }
  size_t get_quadtree_max_items(void) const { return quadtree_max_items_; }
  size_t get_quadtree_max_depth(void) const { return quadtree_max_depth_; }

  private:

  /*------------------------------------------------------------------
  | Enforce singleton
  ------------------------------------------------------------------*/
  TQMeshSetup() {}

  /*------------------------------------------------------------------
  | 
  ------------------------------------------------------------------*/
  double quadtree_scale_ { 10000.0 };
  size_t quadtree_max_items_ { 100 };
  size_t quadtree_max_depth_ { 25 };

}; // TQMeshSetup


/*********************************************************************
* 
*********************************************************************/
template <typename T>
class ContainerFactory {
public:

  /*------------------------------------------------------------------
  | Enforce singleton
  ------------------------------------------------------------------*/
  ContainerFactory(const ContainerFactory&) = delete;
  ContainerFactory& operator=(const ContainerFactory&) = delete;

  /*------------------------------------------------------------------
  | Build a new container, after initializing QuadTree properties
  ------------------------------------------------------------------*/
  static Container<T> build_container() 
  {
    TQMeshSetup& setup = TQMeshSetup::get_instance();

    QuadTreeBuilder<T, double>::get_instance()
      .set_scale( setup.get_quadtree_scale() )
      .set_max_items( setup.get_quadtree_max_items() )
      .set_max_depth( setup.get_quadtree_max_depth() );

    Container<T> container {};

    return container;
  }

private:

  /*------------------------------------------------------------------
  | Enforce singleton
  ------------------------------------------------------------------*/
  ContainerFactory() {}

}; // ContainerFactory 


} // namespace TQMesh 
