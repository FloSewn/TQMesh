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

#include "Vertex.h"
#include "Triangle.h"
#include "Quad.h"
#include "Facet.h"

namespace TQMesh {

using namespace CppUtils;

/*********************************************************************
* This class contains all functions that are required to verify 
* different entitis for the mesh generation process
*********************************************************************/
class EntityChecks
{
public:

  /*------------------------------------------------------------------
  | Check if the provided domain is valid for the mesh generation
  ------------------------------------------------------------------*/
  template<typename Domain>
  static inline bool check_domain_validity(const Domain& domain) 
  {
    // Check if any boundaries are defined 
    if ( domain.size() < 1 )
    {
      LOG(ERROR) << "Invalid domain: No boundaries are defined.";
      return false;
    }

    // Check if at least one exterior boundary exists
    bool extr_bdry_found = false;

    for ( const auto& boundary : domain )
    {
      if ( boundary->is_exterior() )
      {
        extr_bdry_found = true;
        break;
      }
    }

    if ( !extr_bdry_found )
    {
      LOG(ERROR) << "Invalid domain: No exterior boundary defined.";
      return false;
    }

    // Check if domain boundaries are traversable
    for ( const auto& boundary : domain )
    {
      Edge& e_start = boundary->edges()[0];

      if ( !boundary->is_traversable(e_start, e_start) )
      {
        LOG(ERROR) << "Invalid domain: Boundaries not traversable.";
        return false;
      }
    }

    // Check if domain boundary edges are intersecting
    for ( const auto& boundary : domain )
    {
      if ( boundary->intersects_self() )
      {
        LOG(ERROR) << "Invalid domain: Self-intersecting boundary.";
        return false;
      }

      for ( const auto& nbr_boundary : domain )
      {
        if (nbr_boundary == boundary)
          continue;

        if ( boundary->intersects_edgelist( *nbr_boundary ) )
        {
          LOG(ERROR) << "Invalid domain: Intersection between two boundaries.";
          return false;
        }
      }
    }

    return true;

  } // check_domain_validity()


  /*------------------------------------------------------------------
  | Check the facet-vertex-edge connectivtiy of a given mesh
  ------------------------------------------------------------------*/
  template <typename Mesh>
  static inline bool check_mesh_validity(Mesh& mesh)
  { 
    // Check connectivity for interior edges
    for ( const auto& e_ptr : mesh.interior_edges() )
    {
      const Vertex& v1 = e_ptr->v1();
      const Vertex& v2 = e_ptr->v2();

      bool check_1 = false; bool check_2 = false;

      // Traverse adjacents facets of both vertices 
      // and search for the current interior edge
      for ( auto f : v1.facets() )
      {
        check_1 = ( f->get_edge_index(v1, v2) >= 0 );
        if (check_1) break;
      }

      for ( auto f : v2.facets() )
      {
        check_2 = ( f->get_edge_index(v1, v2) >= 0 );
        if (check_2) break;
      }

      // Edge was not found 
      if (!check_1 || !check_2) 
        return false;
    }

    // Check connectivity for boundary edges
    for ( const auto& e_ptr : mesh.boundary_edges() )
    {
      const Vertex& v1 = e_ptr->v1();
      const Vertex& v2 = e_ptr->v2();

      ASSERT( v1.has_property( VertexProperty::on_boundary ),
        "EntityChecks::check_mesh_validity(): Missing "
        "vertex property \"on_boundary\".");
      ASSERT( v2.has_property( VertexProperty::on_boundary ),
        "EntityChecks::check_mesh_validity(): Missing "
        "vertex property \"on_boundary\".");

      bool check_1 = false; bool check_2 = false;

      // Traverse adjacents facets of both vertices 
      // and search for the current interior edge
      for ( auto f : v1.facets() )
      {
        check_1 = ( f->get_edge_index(v1, v2) >= 0 );
        if (check_1) break;
      }

      for ( auto f : v2.facets() )
      {
        check_2 = ( f->get_edge_index(v1, v2) >= 0 );
        if (check_2) break;
      }

      // Edge was not found 
      if (!check_1 || !check_2) 
        return false;
    }

    return true;

  } // EntityChecks::check_mesh_validity()



  /*------------------------------------------------------------------
  | Check a given advancing front structure
  ------------------------------------------------------------------*/
  template <typename Front>
  static inline bool check_front_validity(Front& front)
  {


    return true;

  } // EntityChecks::check_front_validity()


private:
  /*------------------------------------------------------------------
  | We hide the constructor, since this class acts only as container
  | for static inline functions
  ------------------------------------------------------------------*/
  EntityChecks() = default;
  ~EntityChecks() {};

}; // EntityChecks

} // namespace TQMesh
