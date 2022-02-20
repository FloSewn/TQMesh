#include <iostream>
#include <vector>
#include <cassert>
#include <iomanip>   
#include <cmath>
#include <cstdlib>

#include "run_tests.h"

#include "utils.h"
#include "Vec2.h"
#include "QTree.h"
#include "Timer.h"
#include "geometry.h"

namespace QTreeTests
{

using namespace TQMesh::TQUtils;

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
* A simple vertex class that is stored in the QuadTree structure
*********************************************************************/
template<typename T>
class VertexType
{
public:
  VertexType(T x, T y) : xy_ {x, y} {}

  const Vec2<T> xy() const { return xy_; }


private:
  Vec2<T>                   xy_;
}; 

using Vertex = VertexType<double>;
using std::unique_ptr;
using std::vector;
using std::list;
using std::make_unique;


/*********************************************************************
* A simple benchmark for the qtree structure
*********************************************************************/
void QTree_benchmark(int n, int r, size_t imax, size_t dmax, 
                     bool brute_force, bool export_qtree=false)
{
  /*------------------------------------------------------------------
  | Initialize structure, vertex container and add vertices
  ------------------------------------------------------------------*/
  list<unique_ptr<Vertex>> vertices;
  Timer timer {};
  SpiralFunction fun( 2.5, -3.4, 6.5);

  const Vec2d center { 0.0, 0.0 };
  double scale       { 200.0 };
  size_t max_item    { imax };
  size_t max_depth   { dmax };
  QTree<Vertex,double> 
    qtree { scale, max_item, max_depth, center };

  constexpr double t0 =  0.0;
  constexpr double t1 =  5.0 * M_PI;
  const double dt = (t1-t0) / static_cast<double>(n);

  double x, y;
  int item_count_qt, item_count_bf;

  // Rectangle size for vertex search (scales with number of vertices)
  double dx = scale / static_cast<double>(n) ;
  double dy = scale / static_cast<double>(n) ;

  
  /*------------------------------------------------------------------
  | Create vertices
  ------------------------------------------------------------------*/
  timer.count();

  for (int i = 0; i < n; ++i)
  {
    const double t = t0 + i*dt;
    fun.eval(t, x, y);
    vertices.push_back( make_unique<Vertex>( x, y ) );
  }


  /*------------------------------------------------------------------
  | Add vertices to qtree
  ------------------------------------------------------------------*/
  timer.count();

  for ( const auto& v_ptr : vertices )
    qtree.add( v_ptr.get() );

  /*------------------------------------------------------------------
  | Find objects within rectangle using qtree
  ------------------------------------------------------------------*/
  timer.count();

  item_count_qt = 0;

  for (int i = 0; i < n; ++i)
  {
    const double t = t0 + i*dt;
    fun.eval( t, x, y );

    vector<const Vertex*> items {};
    qtree.get_items( {x-dx,y-dy}, {x+dx,y+dy}, items );

    item_count_qt += items.size();
  }

  /*------------------------------------------------------------------
  | Find objects within rectangle using brute force
  ------------------------------------------------------------------*/
  timer.count();

  if (brute_force)
  {
    item_count_bf = 0;

    for (int i = 0; i < n; ++i)
    {
      const double t = t0 + i*dt;
      fun.eval( t, x, y );

      for ( const auto& v_ptr : vertices )
      {
        bool in_rect = TQGeom::in_on_rect( 
            (*v_ptr).xy(), {x-dx,y-dy}, {x+dx,y+dy} );

        if ( in_rect )
          ++item_count_bf;
      }
    }

    // Assert that qtree found all items
    ASSERT( (item_count_qt == item_count_bf),
            "QTree::get_items() failed.");
  }

  /*------------------------------------------------------------------
  | Randomly remove vertices 
  ------------------------------------------------------------------*/
  timer.count();

  if ( r > 0)
  {
    std::srand(123);

    for (int i = 0; i < r; ++i)
    {
      auto it = vertices.begin();
      int pick = (std::rand() % (n-i));
      std::advance(it, pick);

      qtree.remove( it->get() );
      vertices.erase( it );
    }

    ASSERT( (qtree.size() == n-r), 
            "QTree::remove() failed.");
  }


  /*------------------------------------------------------------------
  | Repeat serach benchmark, to verify that removal 
  | of objects succeeded
  | AGAIN: Find objects within rectangle using qtree 
  ------------------------------------------------------------------*/
  timer.count();
  
  if ( r > 0)
  {
    item_count_qt = 0;

    for (int i = 0; i < n; ++i)
    {
      const double t = t0 + i*dt;
      fun.eval( t, x, y );

      std::vector<const Vertex*> items {};
      qtree.get_items( {x-dx,y-dy}, {x+dx,y+dy}, items );

      item_count_qt += items.size();
    }
  }


  /*------------------------------------------------------------------
  | AGAIN: Find objects within rectangle using brute force
  ------------------------------------------------------------------*/
  timer.count();

  if (brute_force && r > 0)
  {
    item_count_bf = 0;

    for (int i = 0; i < n; ++i)
    {
      const double t = t0 + i*dt;
      fun.eval( t, x, y );

      for ( const auto& v_ptr : vertices )
      {
        bool in_rect = TQGeom::in_on_rect( 
            (*v_ptr).xy(), {x-dx,y-dy}, {x+dx,y+dy} );

        if ( in_rect )
          ++item_count_bf;
      }
    }
    // Assert that qtree found all items
    ASSERT( (item_count_qt == item_count_bf),
            "QTree::get_items() failed.");
  }


  /*------------------------------------------------------------------
  | Finalize benchmark - output times to user
  ------------------------------------------------------------------*/
  timer.count();

  // User output
  MSG( "\n----------------------------------" );
  MSG( "n=" << n << ", r=" << r 
      << ", max. items= " << imax << ", max. depth=" << dmax );
  MSG( std::setprecision(7) << std::fixed 
      << "dx=" << dx << ", dy=" << dy);
  MSG( "Vertex initialization            : " 
            << timer.delta(0) << "s" );
  MSG( "QTree initialization             : " 
            << timer.delta(1) << "s" );
  MSG( "QTree search                     : " 
            << timer.delta(2) << "s" );
  if (brute_force)
  {
    MSG( "BruteForce search                : " 
              << timer.delta(3) << "s" );
    MSG( "QTree speedup                    : "
              << timer.delta(3) / timer.delta(2) );
  }
  if ( r > 0 )
  {
    MSG( "Data removal                     : "
              << timer.delta(4) );
    MSG( "QTree search after removal       : " 
              << timer.delta(5) << "s" );
  }
  if ( brute_force && r > 0 )
  {
    MSG( "BruteForce search after removal  : " 
              << timer.delta(6) << "s" );
    MSG( "QTree speedup after removal      : "
              << timer.delta(6) / timer.delta(5) );
  }


  /*------------------------------------------------------------------
  | 
  ------------------------------------------------------------------*/
  if (export_qtree)
  {
    MSG("EXPORT QTREE STRUCTURE");
    std::cout << "VERTICES " << vertices.size() << std::endl;
    for ( const auto& v_ptr : vertices )
      std::cout << std::setprecision(5) << std::fixed 
                << (*v_ptr).xy().x << "," 
                << (*v_ptr).xy().y << std::endl;

    std::cout << "QTREE-LEAFS " << qtree.n_leafs() << std::endl;
    std::cout << qtree;
  }

  
} // QTree_benchmark()



/*********************************************************************
* Test QTree Constructor
*********************************************************************/
void Test_QTree_constructor()
{
  const Vec2d center { 0.0, 0.0 };
  double scale       { 1.0 };
  size_t max_item    { 10 };
  size_t max_depth   { 3 };

  QTree<Vertex,double> 
    qtree { scale, max_item, max_depth, center };

  ASSERT( ( EQ(qtree.lowleft().x, -0.5) ),
          "QTree constructor failed." );
  ASSERT( ( EQ(qtree.lowleft().y, -0.5) ),
          "QTree constructor failed." );

  ASSERT( ( EQ(qtree.upright().x,  0.5) ),
          "QTree constructor failed." );
  ASSERT( ( EQ(qtree.upright().y,  0.5) ),
          "QTree constructor failed." );

  ASSERT( ( qtree.size() == 0 ),
          "QTree constructor failed." );

  ASSERT(!( qtree.children()[0] ),
          "QTree constructor failed." );
  ASSERT(!( qtree.children()[1] ),
          "QTree constructor failed." );
  ASSERT(!( qtree.children()[2] ),
          "QTree constructor failed." );
  ASSERT(!( qtree.children()[3] ),
          "QTree constructor failed." );

  DBG_MSG("Tests for QTree constructor succeeded");

} // Test_QTree_constructor()

/*********************************************************************
* Test QTree::add()
*********************************************************************/
void Test_QTree_add(bool export_qtree=false)
{
  /*------------------------------------------------------------------
  | Initialize structure and vertex container
  ------------------------------------------------------------------*/
  const Vec2d center { 0.0, 0.0 };
  double scale       { 1.0 };
  size_t max_item    { 3 };
  size_t max_depth   { 5 };

  QTree<Vertex,double> 
    qtree { scale, max_item, max_depth, center };

  // Vertex-container
  vector<unique_ptr<Vertex>> vertices;

  vertices.push_back( make_unique<Vertex>( 0.1, 0.1 ) );
  vertices.push_back( make_unique<Vertex>( 0.2, 0.2 ) );
  vertices.push_back( make_unique<Vertex>( 0.4, 0.2 ) );
  vertices.push_back( make_unique<Vertex>( 0.4, 0.3 ) );

  /*------------------------------------------------------------------
  | Add vertices
  ------------------------------------------------------------------*/
  qtree.add( vertices[0].get() );
  qtree.add( vertices[1].get() );
  qtree.add( vertices[2].get() );

  ASSERT(!( qtree.split() ),
          "QTree::add() failed." );

  /*------------------------------------------------------------------
  | QTree should split 
  ------------------------------------------------------------------*/
  qtree.add( vertices[3].get() );

  ASSERT( ( qtree.split() ),
          "QTree::add() failed." );
  ASSERT( ( qtree.size() == 4 ),
          "QTree::add() failed." );

  ASSERT( ( qtree.children()[0]->split() ),
          "QTree::add() failed." );

  ASSERT( ( qtree.children()[0]->children()[0]->size() == 1 ),
          "QTree::add() failed." );
  ASSERT( ( qtree.children()[0]->children()[1]->size() == 0 ),
          "QTree::add() failed." );
  ASSERT( ( qtree.children()[0]->children()[2]->size() == 2 ),
          "QTree::add() failed." );
  ASSERT( ( qtree.children()[0]->children()[3]->size() == 1 ),
          "QTree::add() failed." );

  ASSERT( ( qtree.n_leafs() == 7 ),
          "QTree::add() failed." );

  /*------------------------------------------------------------------
  | Check splitting behaviour for collision 
  ------------------------------------------------------------------*/
  vertices.push_back( make_unique<Vertex>( 0.5,-0.5 ) );
  vertices.push_back( make_unique<Vertex>( 0.5,-0.5 ) );
  vertices.push_back( make_unique<Vertex>( 0.5,-0.5 ) );
  vertices.push_back( make_unique<Vertex>( 0.5,-0.5 ) );
  vertices.push_back( make_unique<Vertex>( 0.5,-0.5 ) );

  qtree.add( vertices[4].get() );
  qtree.add( vertices[5].get() );
  qtree.add( vertices[6].get() );
  qtree.add( vertices[7].get() );
  qtree.add( vertices[8].get() );

  // Container in lower right corner should contain more elements
  // than given as maximum number
  ASSERT( ( qtree.children()[3]->children()[3]->children()[3]->
            children()[3]->children()[3]->size() == 5 ),
          "QTree::add() failed." );

  ASSERT( ( qtree.size() == 9 ),
          "QTree::add() failed." );

  if (export_qtree)
  {
    MSG("EXPORT QTREE STRUCTURE");
    std::cout << "VERTICES " << vertices.size() << std::endl;
    for (size_t i = 0; i < vertices.size(); ++i)
      std::cout << std::setprecision(5) << std::fixed 
                << vertices[i]->xy().x << "," 
                << vertices[i]->xy().y << std::endl;

    std::cout << "QTREE-LEAFS " << qtree.n_leafs() << std::endl;
    std::cout << qtree;
  }



  DBG_MSG("Tests for QTree::add() succeeded");

} // Test_QTree_add()

/*********************************************************************
* Test QTree::remove()
*********************************************************************/
void Test_QTree_remove(bool export_qtree=false)
{
  /*------------------------------------------------------------------
  | Initialize structure, vertex container and add vertices
  ------------------------------------------------------------------*/
  const Vec2d center { 0.0, 0.0 };
  double scale       { 1.0 };
  size_t max_item    { 3 };
  size_t max_depth   { 5 };

  QTree<Vertex,double> 
    qtree { scale, max_item, max_depth, center };

  // Vertex-container
  vector<unique_ptr<Vertex>> vertices;

  vertices.push_back( make_unique<Vertex>( 0.1, 0.1 ) );
  qtree.add( vertices[0].get() );

  vertices.push_back( make_unique<Vertex>( 0.2, 0.2 ) );
  qtree.add( vertices[1].get() );

  vertices.push_back( make_unique<Vertex>( 0.4, 0.2 ) );
  qtree.add( vertices[2].get() );

  vertices.push_back( make_unique<Vertex>( 0.4, 0.3 ) );
  qtree.add( vertices[3].get() );

  /*------------------------------------------------------------------
  | Remove / add / remove vertices
  ------------------------------------------------------------------*/
  qtree.remove( vertices[2].get() );
  qtree.add( vertices[2].get() );
  qtree.remove( vertices[2].get() );

  ASSERT( ( qtree.size() == 3 ),
          "QTree::add() failed." );
  ASSERT(!( qtree.split() ),
          "QTree::add() failed." );
  ASSERT( ( qtree.n_leafs() == 1 ),
          "QTree::add() failed." );

  qtree.remove( vertices[1].get() );
  qtree.add( vertices[1].get() );
  qtree.remove( vertices[1].get() );

  qtree.remove( vertices[3].get() );
  qtree.add( vertices[3].get() );
  qtree.remove( vertices[3].get() );

  qtree.remove( vertices[0].get() );
  qtree.add( vertices[0].get() );
  qtree.remove( vertices[0].get() );

  ASSERT( ( qtree.size() == 0 ),
          "QTree::add() failed." );


  if (export_qtree)
  {
    MSG("EXPORT QTREE STRUCTURE");
    std::cout << "VERTICES " << vertices.size() << std::endl;
    for ( const auto& v_ptr : vertices )
      std::cout << std::setprecision(5) << std::fixed 
                << (*v_ptr).xy().x << "," 
                << (*v_ptr).xy().y << std::endl;

    std::cout << "QTREE-LEAFS " << qtree.n_leafs() << std::endl;
    std::cout << qtree;
  }

  DBG_MSG("Tests for QTree::remove() succeeded");

} // Test_QTree_remove()

/*********************************************************************
* Test QTree::get_items()
*********************************************************************/
void Test_QTree_get_items(bool export_qtree=false)
{
  /*------------------------------------------------------------------
  | Initialize structure, vertex container and add vertices
  ------------------------------------------------------------------*/
  const Vec2d center { 0.0, 0.0 };
  double scale       { 1.0 };
  size_t max_item    { 3 };
  size_t max_depth   { 5 };

  QTree<Vertex,double> 
    qtree { scale, max_item, max_depth, center };

  // Vertex-container
  vector<unique_ptr<Vertex>> vertices;

  vertices.push_back( make_unique<Vertex>( 0.1, 0.1 ) );
  qtree.add( vertices[0].get() );

  vertices.push_back( make_unique<Vertex>( 0.2, 0.2 ) );
  qtree.add( vertices[1].get() );

  vertices.push_back( make_unique<Vertex>( 0.4, 0.2 ) );
  qtree.add( vertices[2].get() );

  vertices.push_back( make_unique<Vertex>( 0.4, 0.3 ) );
  qtree.add( vertices[3].get() );

  vertices.push_back( make_unique<Vertex>(-0.2,-0.2 ) );
  qtree.add( vertices[4].get() );

  vertices.push_back( make_unique<Vertex>(-0.4,-0.4 ) );
  qtree.add( vertices[5].get() );

  vertices.push_back( make_unique<Vertex>( 0.0,-0.3 ) );
  qtree.add( vertices[6].get() );

  vertices.push_back( make_unique<Vertex>(-0.1,-0.2 ) );
  qtree.add( vertices[7].get() );

  vertices.push_back( make_unique<Vertex>( 0.1,-0.1 ) );
  qtree.add( vertices[8].get() );


  /*------------------------------------------------------------------
  | Get items in rectangle
  ------------------------------------------------------------------*/
  vector<const Vertex*> r_found {};
  
  size_t n_r = qtree.get_items( {-0.25,-0.25}, {0.25,0.25}, r_found );

  ASSERT( ( n_r == 5 ),
          "QTree::get_items() failed." );
  ASSERT( ( n_r == r_found.size() ),
          "QTree::get_items() failed." );

  /*------------------------------------------------------------------
  | Get items in circle
  ------------------------------------------------------------------*/
  vector<const Vertex*> c_found {};
  
  size_t n_c = qtree.get_items( {0.25,0.25}, 0.25, c_found );

  ASSERT( ( n_c == 4 ),
          "QTree::get_items() failed." );
  ASSERT( ( n_c == c_found.size() ),
          "QTree::get_items() failed." );


  (void) n_r;
  (void) n_c;

  /*------------------------------------------------------------------
  | 
  ------------------------------------------------------------------*/
  if (export_qtree)
  {
    MSG("EXPORT QTREE STRUCTURE");
    std::cout << "VERTICES " << vertices.size() << std::endl;
    for (size_t i = 0; i < vertices.size(); ++i)
      std::cout << std::setprecision(5) << std::fixed 
                << vertices[i]->xy().x << "," 
                << vertices[i]->xy().y << std::endl;

    std::cout << "QTREE-LEAFS " << qtree.n_leafs() << std::endl;
    std::cout << qtree;
  }

  DBG_MSG("Tests for QTree::get_items() succeeded");

} // Test_QTree_get_items()

/*********************************************************************
* Test QTree::add_remove()
*********************************************************************/
void Test_QTree_add_remove(bool export_qtree=false)
{
  /*------------------------------------------------------------------
  | Initialize structure, vertex container and add vertices
  ------------------------------------------------------------------*/
  const Vec2d center { 0.0, 0.0 };
  double scale       { 2.0 };
  size_t max_item    { 3 };
  size_t max_depth   { 5 };

  QTree<Vertex,double> 
    qtree { scale, max_item, max_depth, center };

  // Vertex-container
  vector<unique_ptr<Vertex>> vertices;

  vertices.push_back( make_unique<Vertex>( 0.1, 0.1 ) );
  qtree.add( vertices[0].get() );

  vertices.push_back( make_unique<Vertex>( 0.2, 0.2 ) );
  qtree.add( vertices[1].get() );

  vertices.push_back( make_unique<Vertex>( 0.3, 0.3 ) );
  qtree.add( vertices[2].get() );

  vertices.push_back( make_unique<Vertex>( 0.3, 0.3 ) );
  qtree.add( vertices[3].get() );

  ASSERT( ( qtree.split() ),  "QTree::add() failed." );

  qtree.remove( vertices[3].get() );
  ASSERT( !( qtree.split() ),  "QTree::add() failed." );

  qtree.remove( vertices[2].get() );
  qtree.remove( vertices[1].get() );

  qtree.add( vertices[3].get() );
  qtree.add( vertices[2].get() );
  qtree.add( vertices[1].get() );
  ASSERT( ( qtree.split() ),  "QTree::add() failed." );

  /*------------------------------------------------------------------
  | 
  ------------------------------------------------------------------*/
  if (export_qtree)
  {
    MSG("EXPORT QTREE STRUCTURE");
    std::cout << "VERTICES " << vertices.size() << std::endl;
    for (size_t i = 0; i < vertices.size(); ++i)
      std::cout << std::setprecision(5) << std::fixed 
                << vertices[i]->xy().x << "," 
                << vertices[i]->xy().y << std::endl;

    std::cout << "QTREE-LEAFS " << qtree.n_leafs() << std::endl;
    std::cout << qtree;
  }

  DBG_MSG("Tests for QTree add /remove succeeded");

} // Test_QTree_add_remove()

} // Namespace QTreeTests

/*********************************************************************
* Run QTree tests
*********************************************************************/
void run_qtree_tests(bool benchmark)
{
  MSG("\n#===== QTree tests =====");

  QTreeTests::Test_QTree_constructor();
  QTreeTests::Test_QTree_add( false );
  QTreeTests::Test_QTree_remove( false );
  QTreeTests::Test_QTree_get_items( false );
  QTreeTests::Test_QTree_add_remove( false );

  if ( benchmark )
  {
    MSG("\n#===== QTree benchmarks =====");

    int n = 100;
    size_t e = 5;
    for (size_t i = 0; i < e; ++i)
    {
      QTreeTests::QTree_benchmark(n, n/2, 100, 25, true, (i==8));
      n *= 2;
    }
  }

  /*
  QTreeTests::QTree_benchmark(120000, 0, 5, 10, false, false);
  QTreeTests::QTree_benchmark(120000, 0, 20, 10, false, false);
  QTreeTests::QTree_benchmark(120000, 0, 80, 10, false, false);
  QTreeTests::QTree_benchmark(120000, 0, 240, 10, false, false);
  QTreeTests::QTree_benchmark(120000, 0, 960, 10, false, false);
  */
}
