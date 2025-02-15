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

#include "utils.h"
#include "Vertex.h"
#include "Edge.h"
#include "EdgeList.h"

namespace EdgeListTests 
{
using namespace CppUtils;
using namespace TQMesh;

/*********************************************************************
* Test EdgeList insertion and removal
*
*   x---x---x
*   |       |
*   x---x---x
*
*********************************************************************/
void add_remove()
{
  Container<Vertex> vertices { ContainerFactory<Vertex>::build_container() };

  Vertex& v1 = vertices.push_back( 1.0, 1.0 );
  Vertex& v2 = vertices.push_back( 2.0, 1.0 );
  Vertex& v3 = vertices.push_back( 3.0, 1.0 );
  Vertex& v4 = vertices.push_back( 3.0, 2.0 );
  Vertex& v5 = vertices.push_back( 2.0, 2.0 );
  Vertex& v6 = vertices.push_back( 1.0, 2.0 );

  EdgeList edges{ Orientation::CCW };

  Edge& e1 = edges.add_edge(v1,v2,1);
  Edge& e2 = edges.add_edge(v2,v3,1);
  Edge& e3 = edges.add_edge(v3,v4,2);
  Edge& e4 = edges.add_edge(v4,v5,3);
  Edge& e5 = edges.add_edge(v5,v6,3);
  Edge& e6 = edges.add_edge(v6,v1,4);

  // Check edge adjacency
  int nv = static_cast<int>(vertices.size());

  for (int i = 0; i < nv; i++)
  {
    CHECK( vertices[i].is_adjacent( edges[ MOD(i-1, nv) ] ) );

    CHECK( vertices[i].is_adjacent( edges[ i ] ) );

    CHECK( (edges[i].v1() == vertices[i]) );

    CHECK( (edges[i].v2() == vertices[ MOD(i+1, nv) ]) );
  }

  // Check area
  CHECK( EQ(edges.area(), 2.0) );

  // Remove edges
  edges.remove( e3 );

  CHECK( v2.edges().size() == 2 );
  CHECK( v3.edges().size() == 1 );
  CHECK( v4.edges().size() == 1 );
  CHECK( edges.size() == 5 );

  (void) e1,e2,e3,e4,e5,e6;

} // add_remove()

/*********************************************************************
* Test EdgeList::is_inside() function
*
*   x---x---x
*   |       |
*   x---x---x
*
*********************************************************************/
void is_inside()
{
  Container<Vertex> vertices { ContainerFactory<Vertex>::build_container() };

  Vertex& v1 = vertices.push_back( 1.0, 1.0 );
  Vertex& v2 = vertices.push_back( 2.0, 1.0 );
  Vertex& v3 = vertices.push_back( 3.0, 1.0 );
  Vertex& v4 = vertices.push_back( 3.0, 2.0 );
  Vertex& v5 = vertices.push_back( 2.0, 2.0 );
  Vertex& v6 = vertices.push_back( 1.0, 2.0 );

  EdgeList edges{ Orientation::CCW };

  Edge& e1 = edges.add_edge(v1,v2,1);
  Edge& e2 = edges.add_edge(v2,v3,1);
  Edge& e3 = edges.add_edge(v3,v4,2);
  Edge& e4 = edges.add_edge(v4,v5,3);
  Edge& e5 = edges.add_edge(v5,v6,3);
  Edge& e6 = edges.add_edge(v6,v1,4);


  Vertex& v7 = vertices.push_back( 2.0, 2.0 );
  Vertex& v8 = vertices.push_back( 5.0, 5.0 );

  CHECK( edges.is_inside(v7) );
  
  CHECK( !(edges.is_inside(v8)) );

  CHECK( edges.is_inside(v2) );

  (void) e1,e2,e3,e4,e5,e6;
  (void) v7, v8;

} // is_inside()

/*********************************************************************
* Test EdgeList::split_edge() function
*
*   x---x---x
*   |       |
*   x---x---x
*
*********************************************************************/
void split_edge()
{
  Container<Vertex> vertices { ContainerFactory<Vertex>::build_container() };

  Vertex& v1 = vertices.push_back( 1.0, 1.0 );
  Vertex& v2 = vertices.push_back( 2.0, 1.0 );
  Vertex& v3 = vertices.push_back( 3.0, 1.0 );
  Vertex& v4 = vertices.push_back( 3.0, 2.0 );
  Vertex& v5 = vertices.push_back( 2.0, 2.0 );
  Vertex& v6 = vertices.push_back( 1.0, 2.0 );

  EdgeList edges{ Orientation::CCW };

  Edge& e1 = edges.add_edge(v1,v2,1);
  Edge& e2 = edges.add_edge(v2,v3,1);
  Edge& e3 = edges.add_edge(v3,v4,2);
  Edge& e4 = edges.add_edge(v4,v5,3);
  Edge& e5 = edges.add_edge(v5,v6,3);
  Edge& e6 = edges.add_edge(v6,v1,4);

  auto new_edges = edges.split_edge( e1, vertices, 0.5 );
  Edge* e1_new = new_edges.first;
  Edge* e2_new = new_edges.second;

  CHECK( (edges.size() == 7) );

  CHECK( (e1_new->v1() == v1) );
  CHECK( (e2_new->v2() == v2) );

  const Vec2d xy_new = 0.5 * (v1.xy() + v2.xy());

  CHECK( ((e1_new->v2().xy() - xy_new).norm() < TQ_SMALL) );

  CHECK( ((e2_new->v1().xy() - xy_new).norm() < TQ_SMALL) );


  // Check for correct order of new edges
  auto pos_e2 = e2.pos();
  auto pos_e5 = e5.pos();
  auto pos_e6 = e6.pos();

  auto pos_e1_new = e1_new->pos();
  auto pos_e2_new = e2_new->pos();

  CHECK( ((++pos_e5) == pos_e6) );

  CHECK( ((++pos_e1_new) == pos_e2_new) );

  CHECK( ((++pos_e2_new) == pos_e2) );

  (void) e1,e2,e3,e4,e5,e6;
  (void) xy_new;
  (void) new_edges;
  (void) e1_new, e2_new;
  (void) pos_e2, pos_e5, pos_e6, pos_e1_new, pos_e2_new;

} // split_edge()


} // namespace EdgeListTests


/*********************************************************************
* Run tests for: EdgeList.h
*********************************************************************/
void run_tests_EdgeList()
{
  EdgeListTests::add_remove();
  EdgeListTests::is_inside();
  EdgeListTests::split_edge();

} // run_tests_EdgeList()
