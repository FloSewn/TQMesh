/*
* This source file is part of the tqmesh library.  
* This code was written by Florian Setzwein in 2022, 
* and is covered under the MIT License
* Refer to the accompanying documentation for details
* on usage and license.
*/
#pragma once

#include <algorithm>
#include <limits.h>

#include "VecND.h"
#include "utils.h"
#include "VtkIO.h"
#include "ProgressBar.h"

#include "Vertex.h"
#include "Edge.h"
#include "Facet.h"
#include "Boundary.h"
#include "Front.h"
#include "Domain.h"
#include "Mesh.h"
#include "MeshValidator.h"

namespace TQMesh {
namespace TQAlgorithm {

using namespace CppUtils;

/*********************************************************************
* 
*********************************************************************/
class MeshingFunctions
{
public:
  using VertexVector   = std::vector<Vertex*>;
  using TriVector      = std::vector<Triangle*>;

  /*------------------------------------------------------------------
  | 
  ------------------------------------------------------------------*/
  static inline 
  Vertex& create_base_vertex(const Edge& base_edge, 
                             const Domain& domain,
                             Mesh& mesh,
                             double base_vertex_factor)
  {
    // Half of the factor h for height of equlateral triangle
    // h := sqrt(3) / 2  -   h_fac := h / 2
    constexpr double h_fac = 0.4330127019; 
    const double v_fac = h_fac * base_vertex_factor;

    // Obtain size function value at the centroid of an equlateral
    // triangle, created from the current base edge
    Vec2d c = base_edge.xy() + base_edge.normal() * base_edge.length() * v_fac;
    const double rho = domain.size_function(c);

    // Coordinate of new vertex 
    Vec2d xy = base_edge.xy() + base_edge.normal() * rho;

    return mesh.add_vertex( xy );

  } // create_base_vertex() 

  /*------------------------------------------------------------------
  | This function is part of the advancing front loop.
  | We loop over a given set of <vertex_candidates> and check
  | if we can create a possible triangle with the current base edge
  | (<b1>,<b2>) and the given vertices.
  | If it is possible, new triangles are created and pushed back to 
  | the vector <new_triangles>.
  ------------------------------------------------------------------*/
  static inline 
  void check_vertex_candidates(const VertexVector& vertex_candidates,
                               Edge& base_edge, TriVector& new_triangles,
                               Mesh& mesh, MeshValidator& validator)
  {
    for ( Vertex* v : vertex_candidates )
    {
      // Skip vertices that are not located on the advancing front
      if ( !v->on_front() )
        continue;

      // Skip vertices that are colinear to the current base edge
      if ( orientation( base_edge.v1().xy(), base_edge.v2().xy(), v->xy() )
          == Orientation::CL )
        continue;

      // Create new potential triangle 
      Triangle& t_new = mesh.add_triangle( base_edge.v1(), base_edge.v2(), *v );

      // Check if new potential triangle is valid
      if ( !validator.remove_from_mesh_if_invalid(t_new) )
        new_triangles.push_back( &t_new );
    }

  } // check_vertex_candidates()


}; // MeshingFunctions

} // namespace TQAlgorithm
} // namespace TQMesh
