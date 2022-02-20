#include <iostream>
#include <vector>
#include <list>
#include <cassert>
#include <iomanip>   
#include <cmath>
#include <cstdlib>
#include <memory>       

#include "run_tests.h"

#include "utils.h"
#include "Vec2.h"
#include "Timer.h"
#include "Container.h"

namespace ContainerTests
{

using namespace TQMesh::TQUtils;

class Edge;

/*********************************************************************
* A spiral function for the generation of vertices in the QTree
* benchmark
*********************************************************************/
class SpiralFunction
{
public:
  SpiralFunction(double a, double b, double c)
  : a_ { a }, b_ { b }, c_ { c } {}

  inline Vec2d eval( const double t ) const
  {
    double x = (a_ + b_ * t) * cos(t) + c_ * sin(40.*t);
    double y = (a_ + b_ * t) * sin(t) + c_ * cos(40.*t);
    return { x,y };
  }

  inline void eval( const double t, double& x, double& y)
  {
    x = (a_ + b_ * t) * cos(t) + c_ * sin(40.*t);
    y = (a_ + b_ * t) * sin(t) + c_ * cos(40.*t);
  }

private:
  double a_  { 0.0 };
  double b_  { 0.0 };
  double c_  { 0.0 };

}; 

/*********************************************************************
* A simple vertex class that is stored in the Container
*********************************************************************/
class Vertex
{
public:

  friend Container<Vertex>;
  using List = Container<Vertex>::List; 
  using Iterator = List::iterator;
  using EdgeList = std::list<Edge*>;

  /*------------------------------------------------------------------
  | Constructor 
  ------------------------------------------------------------------*/
  Vertex(double x, double y) : xy_ {x, y} {}

  /*------------------------------------------------------------------
  | Getters 
  ------------------------------------------------------------------*/
  const Vec2d xy() const { return xy_; }
  const Iterator& pos() const { return pos_; }
  const EdgeList& edges() const { return edges_;}
  const Edge& edges(size_t i) const 
  { 
    ASSERT(!( i < 0 || i >= edges_.size() ),
            "Invalid access to vertex edge list." );
    auto iter = edges_.begin();
    std::advance( iter, i );
    ASSERT( *iter, "Edge vertex is nullptr.");
    return *(*iter);
  }

  void add_edge(Edge* e) { edges_.push_back(e); }
  void remove_edge(Edge* e) { edges_.remove(e); }

  /*------------------------------------------------------------------
  | Container functions 
  ------------------------------------------------------------------*/
  bool in_container() const { return in_container_; }

private:
  /*------------------------------------------------------------------
  | Container functions 
  ------------------------------------------------------------------*/
  void container_destructor() {}

  /*------------------------------------------------------------------
  | Vertex attributes 
  ------------------------------------------------------------------*/
  Vec2d         xy_;
  EdgeList      edges_ {};

  Iterator      pos_ {nullptr};
  bool          in_container_;

}; 


/*********************************************************************
* A simple edge class that is stored in the Container
*********************************************************************/
class Edge
{
public:

  friend Container<Edge>;
  using List = Container<Edge>::List;
  using Iterator = List::iterator;

  /*------------------------------------------------------------------
  | Constructor 
  ------------------------------------------------------------------*/
  Edge(Vertex& v1, Vertex& v2) 
  : v1_ {&v1}
  , v2_ {&v2}
  {
    ASSERT((&v1_ && &v2_),
        "Failed to create edge structure due to given nullptr." );
      
    const Vec2d d_xy = v2_->xy() - v1_->xy();

    xy_     = 0.5 * ( v1_->xy() + v2_->xy() );
    length_ = d_xy.length();
    tang_   = d_xy / length_;

    norm_.x = -tang_.y;
    norm_.y =  tang_.x;

    v1_->add_edge( this );
    v2_->add_edge( this );
  }

  /*------------------------------------------------------------------
  | Destructor 
  ------------------------------------------------------------------*/
  ~Edge() {}

  /*------------------------------------------------------------------
  | Getters 
  ------------------------------------------------------------------*/
  const Vec2d xy() const { return xy_; }
  const Iterator& pos() const { return pos_; }

  const Vertex& v1() const { return *v1_; };
  const Vertex& v2() const { return *v2_; };
  Vertex& v1() { return *v1_; };
  Vertex& v2() { return *v2_; };

  double length() const { return length_; }
  const Vec2d& normal() const { return norm_;}
  const Vec2d& tangent() const { return tang_;}

  /*------------------------------------------------------------------
  | Container functions 
  ------------------------------------------------------------------*/
  bool in_container() const { return in_container_; }

private:
  /*------------------------------------------------------------------
  | Container functions 
  ------------------------------------------------------------------*/
  void container_destructor()
  { 
    if (v1_) v1_->remove_edge( this );
    if (v2_) v2_->remove_edge( this );

    v1_ = nullptr;
    v2_ = nullptr;
  }

  /*------------------------------------------------------------------
  | Edge attributes 
  ------------------------------------------------------------------*/
  Vertex*     v1_;
  Vertex*     v2_;

  Vec2d       xy_;
  Vec2d       tang_;
  Vec2d       norm_;
  double      length_;

  Iterator    pos_ {nullptr};
  bool        in_container_; 

}; 

/*********************************************************************
* Vertex ostream operator overload
*********************************************************************/
static std::ostream& operator<<(std::ostream& os, 
                                const Vertex& v)
{ return os << v.xy(); }

/*********************************************************************
* Edge ostream operator overload
*********************************************************************/
static std::ostream& operator<<(std::ostream& os, 
                                const Edge& e)
{ return os << e.v1() << " --> "  << e.v2(); }

/***********************************************************
* Vertex equality operator 
***********************************************************/
static bool operator==(const Vertex& v1, const Vertex& v2)
{ return &v1 == &v2; }
static bool operator!=(const Vertex& v1, const Vertex& v2)
{ return !(v1 == v2); }

/***********************************************************
* Edge equality operator 
***********************************************************/
static bool operator==(const Edge& e1, const Edge& e2)
{ return &e1 == &e2; }
static bool operator!=(const Edge& e1, const Edge& e2)
{ return !(e1 == e2); }


/*********************************************************************
* Test Container::push_back()
*********************************************************************/
void Test_Container_push_back()
{
  /*------------------------------------------------------------------
  | Initialize structure, vertex container and add vertices
  ------------------------------------------------------------------*/
  Container<Vertex> vertices {};

  vertices.push_back( 0.0, 0.0 );
  Vertex& v1 = vertices.push_back( 1.0, 1.0 );
  vertices.push_back( 2.0, 1.0 );
  vertices.push_back( 3.0, 0.0 );
  vertices.push_back( 2.0, 0.0 );

  ASSERT( (v1.in_container()),
          "Container::push_back() failed." );

  ASSERT( (vertices.size() == 5),
          "Container::push_back() failed." );

  ASSERT( (vertices[1] == v1),
          "Container::push_back() failed." );

  // Check Container::back() function
  ASSERT( (vertices.back() == vertices[4]),
          "Container::push_back() failed." );

  // Check Container::front() function
  ASSERT( (vertices.front() == vertices[0]),
          "Container::push_back() failed." );

  (void) v1;

  DBG_MSG("Tests for Container::push_back() succeeded");

} // Test_Container_push_back()

/*********************************************************************
* Test Container::()
*********************************************************************/
void Test_Container_insert()
{
  /*------------------------------------------------------------------
  | Initialize structure, vertex container and add vertices
  ------------------------------------------------------------------*/
  Container<Vertex> vertices {};

  vertices.push_back( 0.0, 0.0 );

  Vertex& v3 = vertices.push_back( 3.0, 1.0 );

  // Insert before vertices[1]
  Vertex& v1 = vertices.insert( vertices[1].pos(), 1.0, 1.0 );
  Vertex& v2 = vertices.insert( v3.pos(), 2.0, 1.0 );

  vertices.push_back( 2.0, 0.0 );

  ASSERT( (vertices.size() == 5),
          "Container::push_back() failed." );

  ASSERT( (vertices[1] == v1),
          "Container::push_back() failed." );
  ASSERT( (vertices[1].xy().x == 1.0 && vertices[1].xy().y == 1.0 ),
          "Container::push_back() failed." );

  ASSERT( (vertices[2] == v2),
          "Container::push_back() failed." );
  ASSERT( (vertices[2].xy().x == 2.0 && vertices[2].xy().y == 1.0 ),
          "Container::push_back() failed." );

  (void) v1;
  (void) v2;
  (void) v3;

  DBG_MSG("Tests for Container::insert() succeeded");

} // Test_Container_insert()

/*********************************************************************
* Test Container::remove()
*********************************************************************/
void Test_Container_remove()
{
  /*------------------------------------------------------------------
  | Initialize structure, vertex container and add vertices
  ------------------------------------------------------------------*/
  Container<Vertex> vertices {};

  vertices.push_back( 0.0, 0.0 );
  Vertex& v1 = vertices.push_back( 1.0, 1.0 );
  vertices.push_back( 2.0, 1.0 );
  vertices.push_back( 3.0, 0.0 );
  Vertex& v4 = vertices.push_back( 2.0, 0.0 );


  vertices.remove( vertices[3] );

  ASSERT( (vertices.size() == 4),
          "Container::push_back() failed." );

  ASSERT( (vertices.waste().size() == 1),
          "Container::push_back() failed." );

  vertices.remove( v1 );

  ASSERT( (vertices.size() == 3),
          "Container::push_back() failed." );

  ASSERT( (vertices.waste().size() == 2),
          "Container::push_back() failed." );

  vertices.remove( v4 );

  ASSERT( (vertices.size() == 2),
          "Container::push_back() failed." );

  vertices.remove( vertices[0] );
  vertices.remove( vertices[0] );

  ASSERT( (vertices.size() == 0),
          "Container::push_back() failed." );

  // Free everything in the garbage collector
  vertices.clear_waste();
  ASSERT( (vertices.waste().size() == 0),
          "Container::push_back() failed." );

  (void) v1;
  (void) v4;

  DBG_MSG("Tests for Container::remove() succeeded");

} // Test_Container_remove()

/*********************************************************************
* Test Container of Vertices and Edges
*********************************************************************/
void Test_Container_Vertices_Edges()
{
  /*------------------------------------------------------------------
  | Initialize structure, vertex container and add vertices
  ------------------------------------------------------------------*/
  Container<Vertex> vertices {};
  Container<Edge>   edges {};

  vertices.push_back( 0.0, 0.0 );
  vertices.push_back( 1.0, 1.0 );
  vertices.push_back( 2.0, 1.0 );
  vertices.push_back( 3.0, 1.0 );
  vertices.push_back( 2.0, 0.0 );

  for (size_t i = 0; i< vertices.size(); i++)
    edges.push_back( vertices[i], vertices[(i+1)%vertices.size()] );

  ASSERT( (edges.size() == 5),
          "Test_Container_Vertices_Edges() failed." );
  ASSERT( (vertices[0].edges().size() == 2),
          "Test_Container_Vertices_Edges() failed." );
  ASSERT( (vertices[0].edges(0) == edges[0] ),
          "Test_Container_Vertices_Edges() failed." );
  ASSERT( (vertices[1].edges(1) == edges[1] ),
          "Test_Container_Vertices_Edges() failed." );
  ASSERT( (vertices[3].edges(0) == edges[2] ),
          "Test_Container_Vertices_Edges() failed." );

  edges.remove( edges[3] );

  ASSERT( (edges.size() == 4),
          "Test_Container_Vertices_Edges() failed." );
  ASSERT( (vertices[0].edges().size() == 2),
          "Test_Container_Vertices_Edges() failed." );
  ASSERT( (vertices[3].edges().size() == 1),
          "Test_Container_Vertices_Edges() failed." );
  ASSERT( (vertices[4].edges().size() == 1),
          "Test_Container_Vertices_Edges() failed." );

  ASSERT( (vertices[3].edges(0) == edges[2] ),
          "Test_Container_Vertices_Edges() failed." );
  ASSERT( (vertices[4].edges(0) == edges[3] ),
          "Test_Container_Vertices_Edges() failed." );

  // Check that pointer equality works
  Edge& e = edges.push_back( vertices[1], vertices[3] );
  ASSERT( (e == edges[4]),
          "Test_Container_Vertices_Edges() failed." );
  ASSERT( (&e == &edges[4]),
          "Test_Container_Vertices_Edges() failed." );

  (void) e;


  DBG_MSG("Tests for Container of Vertices and Edges succeeded");

} // Test_Container_Vertices_Edges()

/*********************************************************************
* Test Container::sort()
*********************************************************************/
void Test_Container_sort()
{
  Container<Vertex> vertices {};

  vertices.push_back( 0.0, 0.0 );
  vertices.push_back( 1.0, 1.0 );
  vertices.push_back( 2.0, 1.0 );
  vertices.push_back( 3.0, 1.0 );
  vertices.push_back( 2.0, 0.0 );

  // Sort vertices by ascending distance to point "x"
  const Vec2d x {1.0,1.0};

  vertices.sort(
  [x] ( Container<Vertex>::value_type& a, 
        Container<Vertex>::value_type& b )
  {
    const double l1 = ((*a).xy()-x).length_squared();
    const double l2 = ((*b).xy()-x).length_squared();
    return ( l1 < l2 );
  });

  double l_old = 0.0;
  for ( const auto& v_ptr : vertices )
  {
    const double l = ((*v_ptr).xy()-x).length();
    ASSERT( (l_old <= l), 
          "Container::sort() failed." );
    l_old = l;
  }

  (void) l_old;



  DBG_MSG("Tests for Container::sort() succeeded");

} // Test_Container_sort()


/*********************************************************************
* Benchmark Container functions
*********************************************************************/
void Container_benchmark(size_t n)
{
  /*------------------------------------------------------------------
  | Initialize structure, vertex container and add vertices
  ------------------------------------------------------------------*/
  Container<Vertex> vertices { 100.0 };
  Timer timer {};
  SpiralFunction fun( 2.5, -3.4, 6.5);

  constexpr double t0 =  0.0;
  constexpr double t1 =  5.0 * M_PI;
  const double dt = (t1-t0) / static_cast<double>(n);

  double x, y;


  /*------------------------------------------------------------------
  | Push back new vertices to the container
  ------------------------------------------------------------------*/
  timer.count();

  for (size_t i = 0; i < n; ++i)
  {
    const double t = t0 + i*dt;
    fun.eval(t, x, y);
    vertices.push_back( x, y );
  }

  /*------------------------------------------------------------------
  | Remove vertices from the container
  ------------------------------------------------------------------*/
  timer.count();

  for (size_t i = 0; i < n; ++i)
  {
    Vertex& v_last = vertices.back();
    vertices.remove( v_last );
  }

  /*------------------------------------------------------------------
  | Finalize benchmark - output times to user
  ------------------------------------------------------------------*/
  timer.count();

  // User output
  MSG( "\n----------------------------------" );
  MSG( "n=" << n );
  MSG( "Container vertex insertion          : " 
            << timer.delta(0) << "s" );
  MSG( "Container vertex removal            : " 
            << timer.delta(1) << "s" );


} // Container_benchmark()

} // Namespace ContainerTests

/*********************************************************************
* run container tests
*********************************************************************/
void run_container_tests(bool benchmark)
{
  MSG("\n#===== Container tests =====");

  ContainerTests::Test_Container_push_back();
  ContainerTests::Test_Container_insert();
  ContainerTests::Test_Container_remove();
  ContainerTests::Test_Container_Vertices_Edges();
  ContainerTests::Test_Container_sort();

  if ( benchmark )
  {
    MSG("\n#===== Container benchmarks =====");

    ContainerTests::Container_benchmark(1000);
    ContainerTests::Container_benchmark(10000);
    ContainerTests::Container_benchmark(100000);
    ContainerTests::Container_benchmark(1000000);
    ContainerTests::Container_benchmark(10000000);
  }


}
