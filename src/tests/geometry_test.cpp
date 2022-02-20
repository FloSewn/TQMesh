#include <iostream>
#include <vector>
#include <cassert>
#include <iomanip>   
#include <cmath>
#include <cstdlib>

#include "run_tests.h"

#include "utils.h"
#include "geometry.h"
#include "Vec2.h"

namespace GeometryTests
{

using namespace TQMesh::TQUtils;

/*********************************************************************
* Test TQGeom::orientation()
*********************************************************************/
void Test_TQGeom_orientation()
{
  /*------------------------------------------------------------------
  | Double
  ------------------------------------------------------------------*/
  const Vec2d p1 { -0.1, -0.1 };
  const Vec2d q1 { -0.2, -0.3 };
  const Vec2d r1 {  0.5,  0.2 };

  ASSERT( (TQGeom::orientation(p1,q1,r1) == TQGeom::Orientation::CCW),
          "TQGeom::orientation() failed." );

  ASSERT( (TQGeom::orientation(p1,r1,q1) == TQGeom::Orientation::CW),
          "TQGeom::orientation() failed." );

  ASSERT( (TQGeom::orientation(p1,r1,r1) == TQGeom::Orientation::CL),
          "TQGeom::orientation() failed." );


  /*------------------------------------------------------------------
  | Float
  ------------------------------------------------------------------*/
  const Vec2f p2 { -0.1f, -0.1f };
  const Vec2f q2 { -0.2f, -0.3f };
  const Vec2f r2 {  0.5f,  0.2f };

  ASSERT( (TQGeom::orientation(p2,q2,r2) == TQGeom::Orientation::CCW),
          "TQGeom::orientation() failed." );

  ASSERT( (TQGeom::orientation(p2,r2,q2) == TQGeom::Orientation::CW),
          "TQGeom::orientation() failed." );

  ASSERT( (TQGeom::orientation(p2,r2,r2) == TQGeom::Orientation::CL),
          "TQGeom::orientation() failed." );


  /*------------------------------------------------------------------
  | Int
  ------------------------------------------------------------------*/
  const Vec2i p3 { -1, -1 };
  const Vec2i q3 { -2, -3 };
  const Vec2i r3 {  5,  2 };

  ASSERT( (TQGeom::orientation(p3,q3,r3) == TQGeom::Orientation::CCW),
          "TQGeom::orientation() failed." );

  ASSERT( (TQGeom::orientation(p3,r3,q3) == TQGeom::Orientation::CW),
          "TQGeom::orientation() failed." );

  ASSERT( (TQGeom::orientation(p3,r3,r3) == TQGeom::Orientation::CL),
          "TQGeom::orientation() failed." );

  DBG_MSG("Tests for TQGeom::orientation() succeeded");

} // Test_TQGeom_orientation()

/*********************************************************************
* Test TQGeom::is_left()
*********************************************************************/
void Test_TQGeom_is_left()
{
  /*------------------------------------------------------------------
  | Double
  ------------------------------------------------------------------*/
  const Vec2d p1 { -0.1, -0.1 };
  const Vec2d q1 {  0.1,  0.1 };
  const Vec2d r1 { -0.1,  0.1 };
  const Vec2d s1 {  0.1, -0.1 };
  const Vec2d t1 {  0.0, -0.0 };

  ASSERT( (TQGeom::is_left(p1,q1,r1) ),
          "TQGeom::is_left() failed." );

  ASSERT( !(TQGeom::is_left(p1,q1,s1) ),
          "TQGeom::is_left() failed." );

  ASSERT( !(TQGeom::is_left(p1,q1,t1) ),
          "TQGeom::is_left() failed." );

  /*------------------------------------------------------------------
  | Float
  ------------------------------------------------------------------*/
  const Vec2f p2 { -0.1f, -0.1f };
  const Vec2f q2 {  0.1f,  0.1f };
  const Vec2f r2 { -0.1f,  0.1f };
  const Vec2f s2 {  0.1f, -0.1f };
  const Vec2f t2 {  0.0f, -0.0f };

  ASSERT( (TQGeom::is_left(p2,q2,r2) ),
          "TQGeom::is_left() failed." );

  ASSERT( !(TQGeom::is_left(p2,q2,s2) ),
          "TQGeom::is_left() failed." );

  ASSERT( !(TQGeom::is_left(p2,q2,t2) ),
          "TQGeom::is_left() failed." );

  /*------------------------------------------------------------------
  | Integer
  ------------------------------------------------------------------*/
  const Vec2i p3 { -1, -1 };
  const Vec2i q3 {  1,  1 };
  const Vec2i r3 { -1,  1 };
  const Vec2i s3 {  1, -1 };
  const Vec2i t3 {  0, -0 };

  ASSERT( (TQGeom::is_left(p3,q3,r3) ),
          "TQGeom::is_left() failed." );

  ASSERT( !(TQGeom::is_left(p3,q3,s3) ),
          "TQGeom::is_left() failed." );

  ASSERT( !(TQGeom::is_left(p3,q3,t3) ),
          "TQGeom::is_left() failed." );

  DBG_MSG("Tests for TQGeom::is_left() succeeded");

} // Test_TQGeom_is_left()

/*********************************************************************
* Test TQGeom::is_lefton()
*********************************************************************/
void Test_TQGeom_is_lefton()
{
  /*------------------------------------------------------------------
  | Double
  ------------------------------------------------------------------*/
  const Vec2d p1 { -0.1, -0.1 };
  const Vec2d q1 {  0.1,  0.1 };
  const Vec2d r1 { -0.1,  0.1 };
  const Vec2d s1 {  0.1, -0.1 };
  const Vec2d t1 {  0.0,  0.0 };

  ASSERT( (TQGeom::is_lefton(p1,q1,r1) ),
          "TQGeom::is_lefton() failed." );

  ASSERT( !(TQGeom::is_lefton(p1,q1,s1) ),
          "TQGeom::is_lefton() failed." );

  ASSERT( (TQGeom::is_lefton(p1,q1,t1) ),
          "TQGeom::is_lefton() failed." );

  /*------------------------------------------------------------------
  | Float
  ------------------------------------------------------------------*/
  const Vec2d p2 { -0.1f, -0.1f };
  const Vec2d q2 {  0.1f,  0.1f };
  const Vec2d r2 { -0.1f,  0.1f };
  const Vec2d s2 {  0.1f, -0.1f };
  const Vec2d t2 {  0.0f,  0.0f };

  ASSERT( (TQGeom::is_lefton(p2,q2,r2) ),
          "TQGeom::is_lefton() failed." );

  ASSERT( !(TQGeom::is_lefton(p2,q2,s2) ),
          "TQGeom::is_lefton() failed." );

  ASSERT( (TQGeom::is_lefton(p2,q2,t2) ),
          "TQGeom::is_lefton() failed." );

  /*------------------------------------------------------------------
  | Integer
  ------------------------------------------------------------------*/
  const Vec2d p3 { -1, -1 };
  const Vec2d q3 {  1,  1 };
  const Vec2d r3 { -1,  1 };
  const Vec2d s3 {  1, -1 };
  const Vec2d t3 {  0,  0 };

  ASSERT( (TQGeom::is_lefton(p3,q3,r3) ),
          "TQGeom::is_lefton() failed." );

  ASSERT( !(TQGeom::is_lefton(p3,q3,s3) ),
          "TQGeom::is_lefton() failed." );

  ASSERT( (TQGeom::is_lefton(p3,q3,t3) ),
          "TQGeom::is_lefton() failed." );

  DBG_MSG("Tests for TQGeom::is_lefton() succeeded");

} // Test_TQGeom_is_lefton()

/*********************************************************************
* Test TQGeom::in_segment()
*********************************************************************/
void Test_TQGeom_in_segment()
{
  /*------------------------------------------------------------------
  | Double
  ------------------------------------------------------------------*/
  const Vec2d p1 { -0.1, -0.1 };
  const Vec2d q1 {  0.1,  0.1 };
  const Vec2d r1 { -0.1,  0.1 };
  const Vec2d s1 {  0.1, -0.1 };
  const Vec2d t1 {  0.2,  0.2 };
  const Vec2d u1 { -0.2, -0.2 };
  const Vec2d v1 {  0.0,  0.0 };

  ASSERT( !(TQGeom::in_segment(p1,q1,r1) ),
          "TQGeom::in_segment() failed." );
  ASSERT( !(TQGeom::in_segment(p1,q1,s1) ),
          "TQGeom::in_segment() failed." );
  ASSERT( !(TQGeom::in_segment(p1,q1,t1) ),
          "TQGeom::in_segment() failed." );
  ASSERT( !(TQGeom::in_segment(p1,q1,u1) ),
          "TQGeom::in_segment() failed." );
  ASSERT( (TQGeom::in_segment(p1,q1,v1) ),
          "TQGeom::in_segment() failed." );
  ASSERT( !(TQGeom::in_segment(p1,q1,p1) ),
          "TQGeom::in_segment() failed." );
  ASSERT( !(TQGeom::in_segment(p1,q1,q1) ),
          "TQGeom::in_segment() failed." );

  /*------------------------------------------------------------------
  | Float
  ------------------------------------------------------------------*/
  const Vec2f p2 { -0.1f, -0.1f };
  const Vec2f q2 {  0.1f,  0.1f };
  const Vec2f r2 { -0.1f,  0.1f };
  const Vec2f s2 {  0.1f, -0.1f };
  const Vec2f t2 {  0.2f,  0.2f };
  const Vec2f u2 { -0.2f, -0.2f };
  const Vec2f v2 {  0.0f,  0.0f };

  ASSERT( !(TQGeom::in_segment(p2,q2,r2) ),
          "TQGeom::in_segment() failed." );
  ASSERT( !(TQGeom::in_segment(p2,q2,s2) ),
          "TQGeom::in_segment() failed." );
  ASSERT( !(TQGeom::in_segment(p2,q2,t2) ),
          "TQGeom::in_segment() failed." );
  ASSERT( !(TQGeom::in_segment(p2,q2,u2) ),
          "TQGeom::in_segment() failed." );
  ASSERT( (TQGeom::in_segment(p2,q2,v2) ),
          "TQGeom::in_segment() failed." );
  ASSERT( !(TQGeom::in_segment(p2,q2,p2) ),
          "TQGeom::in_segment() failed." );
  ASSERT( !(TQGeom::in_segment(p2,q2,q2) ),
          "TQGeom::in_segment() failed." );

  /*------------------------------------------------------------------
  | Integer
  ------------------------------------------------------------------*/
  const Vec2i p3 { -1, -1 };
  const Vec2i q3 {  1,  1 };
  const Vec2i r3 { -1,  1 };
  const Vec2i s3 {  1, -1 };
  const Vec2i t3 {  2,  2 };
  const Vec2i u3 { -2, -2 };
  const Vec2i v3 {  0,  0 };

  ASSERT( !(TQGeom::in_segment(p3,q3,r3) ),
          "TQGeom::in_segment() failed." );
  ASSERT( !(TQGeom::in_segment(p3,q3,s3) ),
          "TQGeom::in_segment() failed." );
  ASSERT( !(TQGeom::in_segment(p3,q3,t3) ),
          "TQGeom::in_segment() failed." );
  ASSERT( !(TQGeom::in_segment(p3,q3,u3) ),
          "TQGeom::in_segment() failed." );
  ASSERT( (TQGeom::in_segment(p3,q3,v3) ),
          "TQGeom::in_segment() failed." );
  ASSERT( !(TQGeom::in_segment(p3,q3,p3) ),
          "TQGeom::in_segment() failed." );
  ASSERT( !(TQGeom::in_segment(p3,q3,q3) ),
          "TQGeom::in_segment() failed." );

  DBG_MSG("Tests for TQGeom::in_segment() succeeded");

} // Test_TQGeom_in_segment()

/*********************************************************************
* Test TQGeom::in_on_segment()
*********************************************************************/
void Test_TQGeom_in_on_segment()
{
  /*------------------------------------------------------------------
  | Double
  ------------------------------------------------------------------*/
  const Vec2d p1 { -0.1, -0.1 };
  const Vec2d q1 {  0.1,  0.1 };
  const Vec2d r1 { -0.1,  0.1 };
  const Vec2d s1 {  0.1, -0.1 };
  const Vec2d t1 {  0.2,  0.2 };
  const Vec2d u1 { -0.2, -0.2 };
  const Vec2d v1 {  0.0,  0.0 };

  ASSERT( !(TQGeom::in_on_segment(p1,q1,r1) ),
          "TQGeom::in_on_segment() failed." );
  ASSERT( !(TQGeom::in_on_segment(p1,q1,s1) ),
          "TQGeom::in_on_segment() failed." );
  ASSERT( !(TQGeom::in_on_segment(p1,q1,t1) ),
          "TQGeom::in_on_segment() failed." );
  ASSERT( !(TQGeom::in_on_segment(p1,q1,u1) ),
          "TQGeom::in_on_segment() failed." );
  ASSERT( (TQGeom::in_on_segment(p1,q1,v1) ),
          "TQGeom::in_on_segment() failed." );
  ASSERT( (TQGeom::in_on_segment(p1,q1,p1) ),
          "TQGeom::in_on_segment() failed." );
  ASSERT( (TQGeom::in_on_segment(p1,q1,q1) ),
          "TQGeom::in_on_segment() failed." );

  /*------------------------------------------------------------------
  | Float
  ------------------------------------------------------------------*/
  const Vec2f p2 { -0.1f, -0.1f };
  const Vec2f q2 {  0.1f,  0.1f };
  const Vec2f r2 { -0.1f,  0.1f };
  const Vec2f s2 {  0.1f, -0.1f };
  const Vec2f t2 {  0.2f,  0.2f };
  const Vec2f u2 { -0.2f, -0.2f };
  const Vec2f v2 {  0.0f,  0.0f };

  ASSERT( !(TQGeom::in_on_segment(p2,q2,r2) ),
          "TQGeom::in_on_segment() failed." );
  ASSERT( !(TQGeom::in_on_segment(p2,q2,s2) ),
          "TQGeom::in_on_segment() failed." );
  ASSERT( !(TQGeom::in_on_segment(p2,q2,t2) ),
          "TQGeom::in_on_segment() failed." );
  ASSERT( !(TQGeom::in_on_segment(p2,q2,u2) ),
          "TQGeom::in_on_segment() failed." );
  ASSERT( (TQGeom::in_on_segment(p2,q2,v2) ),
          "TQGeom::in_on_segment() failed." );
  ASSERT( (TQGeom::in_on_segment(p2,q2,p2) ),
          "TQGeom::in_on_segment() failed." );
  ASSERT( (TQGeom::in_on_segment(p2,q2,q2) ),
          "TQGeom::in_on_segment() failed." );

  /*------------------------------------------------------------------
  | Integer
  ------------------------------------------------------------------*/
  const Vec2i p3 { -1, -1 };
  const Vec2i q3 {  1,  1 };
  const Vec2i r3 { -1,  1 };
  const Vec2i s3 {  1, -1 };
  const Vec2i t3 {  2,  2 };
  const Vec2i u3 { -2, -2 };
  const Vec2i v3 {  0,  0 };

  ASSERT( !(TQGeom::in_on_segment(p3,q3,r3) ),
          "TQGeom::in_on_segment() failed." );
  ASSERT( !(TQGeom::in_on_segment(p3,q3,s3) ),
          "TQGeom::in_on_segment() failed." );
  ASSERT( !(TQGeom::in_on_segment(p3,q3,t3) ),
          "TQGeom::in_on_segment() failed." );
  ASSERT( !(TQGeom::in_on_segment(p3,q3,u3) ),
          "TQGeom::in_on_segment() failed." );
  ASSERT( (TQGeom::in_on_segment(p3,q3,v3) ),
          "TQGeom::in_on_segment() failed." );
  ASSERT( (TQGeom::in_on_segment(p3,q3,p3) ),
          "TQGeom::in_on_segment() failed." );
  ASSERT( (TQGeom::in_on_segment(p3,q3,q3) ),
          "TQGeom::in_on_segment() failed." );


  DBG_MSG("Tests for TQGeom::in_on_segment() succeeded");

} // Test_TQGeom_in_on_segment()

/*********************************************************************
* Test TQGeom::line_line_intersection
*********************************************************************/
void Test_TQGeom_line_line_intersection()
{
  /*------------------------------------------------------------------
  | Double
  ------------------------------------------------------------------*/
  const Vec2d p1 {-0.1, -0.1};
  const Vec2d q1 { 0.1,  0.1};

  const Vec2d r1 {-0.1,  0.1};
  const Vec2d s1 { 0.1, -0.1};

  const Vec2d t1 { 0.0,  0.0};

  const Vec2d u1 {-0.1,  0.0};
  const Vec2d v1 { 0.1,  0.2};


  ASSERT( !( TQGeom::line_line_intersection(p1,q1,p1,q1) ),
          "TQGeom::line_line_intersection() failed." );
  ASSERT(  ( TQGeom::line_line_intersection(p1,q1,r1,s1) ),
          "TQGeom::line_line_intersection() failed." );
  ASSERT(  ( TQGeom::line_line_intersection(p1,q1,t1,q1) ),
          "TQGeom::line_line_intersection() failed." );
  ASSERT( !( TQGeom::line_line_intersection(p1,q1,u1,v1) ),
          "TQGeom::line_line_intersection() failed." );
  ASSERT( !( TQGeom::line_line_intersection(p1,q1,q1,v1) ),
          "TQGeom::line_line_intersection() failed." );


  /*------------------------------------------------------------------
  | Float
  ------------------------------------------------------------------*/
  const Vec2f p2 {-0.1f, -0.1f};
  const Vec2f q2 { 0.1f,  0.1f};

  const Vec2f r2 {-0.1f,  0.1f};
  const Vec2f s2 { 0.1f, -0.1f};

  const Vec2f t2 { 0.0f,  0.0f};

  const Vec2f u2 {-0.1f,  0.0f};
  const Vec2f v2 { 0.1f,  0.2f};


  ASSERT( !( TQGeom::line_line_intersection(p2,q2,p2,q2) ),
          "TQGeom::line_line_intersection() failed." );
  ASSERT(  ( TQGeom::line_line_intersection(p2,q2,r2,s2) ),
          "TQGeom::line_line_intersection() failed." );
  ASSERT(  ( TQGeom::line_line_intersection(p2,q2,t2,q2) ),
          "TQGeom::line_line_intersection() failed." );
  ASSERT( !( TQGeom::line_line_intersection(p2,q2,u2,v2) ),
          "TQGeom::line_line_intersection() failed." );
  ASSERT( !( TQGeom::line_line_intersection(p2,q2,q2,v2) ),
          "TQGeom::line_line_intersection() failed." );


  /*------------------------------------------------------------------
  | Integer
  ------------------------------------------------------------------*/
  const Vec2i p3 {-1, -1};
  const Vec2i q3 { 1,  1};

  const Vec2i r3 {-1,  1};
  const Vec2i s3 { 1, -1};

  const Vec2i t3 { 0,  0};

  const Vec2i u3 {-1,  0};
  const Vec2i v3 { 1,  2};


  ASSERT( !( TQGeom::line_line_intersection(p3,q3,p3,q3) ),
          "TQGeom::line_line_intersection() failed." );
  ASSERT(  ( TQGeom::line_line_intersection(p3,q3,r3,s3) ),
          "TQGeom::line_line_intersection() failed." );
  ASSERT(  ( TQGeom::line_line_intersection(p3,q3,t3,q3) ),
          "TQGeom::line_line_intersection() failed." );
  ASSERT( !( TQGeom::line_line_intersection(p3,q3,u3,v3) ),
          "TQGeom::line_line_intersection() failed." );
  ASSERT( !( TQGeom::line_line_intersection(p3,q3,q3,v3) ),
          "TQGeom::line_line_intersection() failed." );

  DBG_MSG("Tests for TQGeom::line_line_intersection() succeeded");

} // Test_TQGeom_line_line_intersection() */

/*********************************************************************
* Test TQGeom::line_tri_intersection
*********************************************************************/
void Test_TQGeom_line_tri_intersection()
{
  /*------------------------------------------------------------------
  | Double
  ------------------------------------------------------------------*/
  const Vec2d p1 {-0.1, -0.1};
  const Vec2d q1 { 0.1,  0.1};

  const Vec2d r1 {-0.1,  0.1};
  const Vec2d s1 { 0.1, -0.1};
  const Vec2d t1 { 0.2,  0.2};

  const Vec2d u1 { 0.2, -0.1};
  const Vec2d v1 { 0.3,  0.1};
  const Vec2d w1 { 0.0,  0.0};

  ASSERT(  ( TQGeom::line_tri_intersection(p1,q1, r1,s1,t1) ),
          "TQGeom::line_tri_intersection() failed." );
  ASSERT(  ( TQGeom::line_tri_intersection(p1,q1, r1,t1,s1) ),
          "TQGeom::line_tri_intersection() failed." );
  ASSERT( !( TQGeom::line_tri_intersection(u1,v1, r1,t1,s1) ),
          "TQGeom::line_tri_intersection() failed." );
  ASSERT( !( TQGeom::line_tri_intersection(r1,s1, r1,t1,s1) ),
          "TQGeom::line_tri_intersection() failed." );
  ASSERT( !( TQGeom::line_tri_intersection(t1,v1, r1,t1,s1) ),
          "TQGeom::line_tri_intersection() failed." );
  ASSERT(  ( TQGeom::line_tri_intersection(t1,p1, r1,t1,s1) ),
          "TQGeom::line_tri_intersection() failed." );
  ASSERT(  ( TQGeom::line_tri_intersection(r1,w1, r1,t1,s1) ),
          "TQGeom::line_tri_intersection() failed." );

  /*------------------------------------------------------------------
  | Float
  ------------------------------------------------------------------*/
  const Vec2f p2 {-0.1f, -0.1f};
  const Vec2f q2 { 0.1f,  0.1f};

  const Vec2f r2 {-0.1f,  0.1f};
  const Vec2f s2 { 0.1f, -0.1f};
  const Vec2f t2 { 0.2f,  0.2f};

  const Vec2f u2 { 0.2f, -0.1f};
  const Vec2f v2 { 0.3f,  0.1f};
  const Vec2f w2 { 0.0f,  0.0f};

  ASSERT(  ( TQGeom::line_tri_intersection(p2,q2, r2,s2,t2) ),
          "TQGeom::line_tri_intersection() failed." );
  ASSERT(  ( TQGeom::line_tri_intersection(p2,q2, r2,t2,s2) ),
          "TQGeom::line_tri_intersection() failed." );
  ASSERT( !( TQGeom::line_tri_intersection(u2,v2, r2,t2,s2) ),
          "TQGeom::line_tri_intersection() failed." );
  ASSERT( !( TQGeom::line_tri_intersection(r2,s2, r2,t2,s2) ),
          "TQGeom::line_tri_intersection() failed." );
  ASSERT( !( TQGeom::line_tri_intersection(t2,v2, r2,t2,s2) ),
          "TQGeom::line_tri_intersection() failed." );
  ASSERT(  ( TQGeom::line_tri_intersection(t2,p2, r2,t2,s2) ),
          "TQGeom::line_tri_intersection() failed." );
  ASSERT(  ( TQGeom::line_tri_intersection(r2,w2, r2,t2,s2) ),
          "TQGeom::line_tri_intersection() failed." );

  /*------------------------------------------------------------------
  | Integer
  ------------------------------------------------------------------*/
  const Vec2i p3 {-1, -1};
  const Vec2i q3 { 1,  1};

  const Vec2i r3 {-1,  1};
  const Vec2i s3 { 1, -1};
  const Vec2i t3 { 2,  2};

  const Vec2i u3 { 2, -1};
  const Vec2i v3 { 3,  1};
  const Vec2i w3 { 0,  0};

  ASSERT(  ( TQGeom::line_tri_intersection(p3,q3, r3,s3,t3) ),
          "TQGeom::line_tri_intersection() failed." );
  ASSERT(  ( TQGeom::line_tri_intersection(p3,q3, r3,t3,s3) ),
          "TQGeom::line_tri_intersection() failed." );
  ASSERT( !( TQGeom::line_tri_intersection(u3,v3, r3,t3,s3) ),
          "TQGeom::line_tri_intersection() failed." );
  ASSERT( !( TQGeom::line_tri_intersection(r3,s3, r3,t3,s3) ),
          "TQGeom::line_tri_intersection() failed." );
  ASSERT( !( TQGeom::line_tri_intersection(t3,v3, r3,t3,s3) ),
          "TQGeom::line_tri_intersection() failed." );
  ASSERT(  ( TQGeom::line_tri_intersection(t3,p3, r3,t3,s3) ),
          "TQGeom::line_tri_intersection() failed." );
  ASSERT(  ( TQGeom::line_tri_intersection(r3,w3, r3,t3,s3) ),
          "TQGeom::line_tri_intersection() failed." );

  DBG_MSG("Tests for TQGeom::line_tri_intersection() succeeded");

} // Test_TQGeom_line_tri_intersection()

/*********************************************************************
* Test TQGeom::line_quad_intersection
*********************************************************************/
void Test_TQGeom_line_quad_intersection()
{
  /*------------------------------------------------------------------
  | Double
  ------------------------------------------------------------------*/
  const Vec2d p1 {-0.1, -0.1};
  const Vec2d q1 { 0.1, -0.1};
  const Vec2d r1 { 0.1,  0.1};
  const Vec2d s1 {-0.1,  0.1};

  const Vec2d t1 { 0.0,  0.0};
  const Vec2d u1 { 0.2,  0.1};
  const Vec2d v1 { 0.2, -0.1};

  ASSERT(  ( TQGeom::line_quad_intersection(t1,u1, p1,q1,r1,s1) ),
          "TQGeom::line_tri_intersection() failed." );
  ASSERT(  ( TQGeom::line_quad_intersection(t1,u1, s1,r1,q1,p1) ),
          "TQGeom::line_tri_intersection() failed." );
  ASSERT( !( TQGeom::line_quad_intersection(u1,v1, p1,q1,r1,s1) ),
          "TQGeom::line_tri_intersection() failed." );
  ASSERT( !( TQGeom::line_quad_intersection(r1,u1, p1,q1,r1,s1) ),
          "TQGeom::line_tri_intersection() failed." );
  ASSERT( !( TQGeom::line_quad_intersection(r1,s1, p1,q1,r1,s1) ),
          "TQGeom::line_tri_intersection() failed." );

  /*------------------------------------------------------------------
  | Float
  ------------------------------------------------------------------*/
  const Vec2f p2 {-0.1f, -0.1f};
  const Vec2f q2 { 0.1f, -0.1f};
  const Vec2f r2 { 0.1f,  0.1f};
  const Vec2f s2 {-0.1f,  0.1f};

  const Vec2f t2 { 0.0f,  0.0f};
  const Vec2f u2 { 0.2f,  0.1f};
  const Vec2f v2 { 0.2f, -0.1f};

  ASSERT(  ( TQGeom::line_quad_intersection(t2,u2, p2,q2,r2,s2) ),
          "TQGeom::line_tri_intersection() failed." );
  ASSERT(  ( TQGeom::line_quad_intersection(t2,u2, s2,r2,q2,p2) ),
          "TQGeom::line_tri_intersection() failed." );
  ASSERT( !( TQGeom::line_quad_intersection(u2,v2, p2,q2,r2,s2) ),
          "TQGeom::line_tri_intersection() failed." );
  ASSERT( !( TQGeom::line_quad_intersection(r2,u2, p2,q2,r2,s2) ),
          "TQGeom::line_tri_intersection() failed." );
  ASSERT( !( TQGeom::line_quad_intersection(r2,s2, p2,q2,r2,s2) ),
          "TQGeom::line_tri_intersection() failed." );

  /*------------------------------------------------------------------
  | Integer
  ------------------------------------------------------------------*/
  const Vec2i p3 {-1, -1};
  const Vec2i q3 { 1, -1};
  const Vec2i r3 { 1,  1};
  const Vec2i s3 {-1,  1};

  const Vec2i t3 { 0,  0};
  const Vec2i u3 { 2,  1};
  const Vec2i v3 { 2, -1};

  ASSERT(  ( TQGeom::line_quad_intersection(t3,u3, p3,q3,r3,s3) ),
          "TQGeom::line_tri_intersection() failed." );
  ASSERT(  ( TQGeom::line_quad_intersection(t3,u3, s3,r3,q3,p3) ),
          "TQGeom::line_tri_intersection() failed." );
  ASSERT( !( TQGeom::line_quad_intersection(u3,v3, p3,q3,r3,s3) ),
          "TQGeom::line_tri_intersection() failed." );
  ASSERT( !( TQGeom::line_quad_intersection(r3,u3, p3,q3,r3,s3) ),
          "TQGeom::line_tri_intersection() failed." );
  ASSERT( !( TQGeom::line_quad_intersection(r3,s3, p3,q3,r3,s3) ),
          "TQGeom::line_tri_intersection() failed." );

  DBG_MSG("Tests for TQGeom::line_quad_intersection() succeeded");

} // Test_TQGeom_line_quad_intersection()

/*********************************************************************
* Test TQGeom::tri_tri_intersection
*********************************************************************/
void Test_TQGeom_tri_tri_intersection()
{
  /*------------------------------------------------------------------
  | Double
  ------------------------------------------------------------------*/
  const Vec2d p1 { 3.0, 2.0 };
  const Vec2d q1 { 9.0, 2.0 };
  const Vec2d r1 { 6.0, 6.0 };

  const Vec2d s1 { 3.0, 2.0 };
  const Vec2d t1 { 6.0, 2.0 };
  const Vec2d u1 { 5.0,-1.0 };

  const Vec2d v1 { 6.0, 4.0 };
  const Vec2d w1 {11.0, 4.0 };
  const Vec2d x1 {10.0, 8.0 };

  const Vec2d a1 { 7.0, 2.0 };
  const Vec2d b1 { 7.0,-1.0 };
  const Vec2d c1 { 9.0,-1.0 };

  const Vec2d d1 {-4.0, 2.0 };
  const Vec2d e1 {-2.0, 2.0 };
  const Vec2d f1 {-3.0, 4.0 };

  const Vec2d g1 { 3.0, 2.0 };
  const Vec2d h1 { 6.0, 6.0 };
  const Vec2d i1 { 2.0, 6.0 };

  ASSERT( (TQGeom::tri_tri_intersection( p1,q1,r1, s1,t1,u1 )),
          "TQGeom::tri_tri_intersection() failed." );
  ASSERT( (TQGeom::tri_tri_intersection( p1,q1,r1, v1,w1,x1 )),
          "TQGeom::tri_tri_intersection() failed." );
  ASSERT( (TQGeom::tri_tri_intersection( p1,q1,r1, a1,b1,c1 )),
          "TQGeom::tri_tri_intersection() failed." );
  ASSERT(!(TQGeom::tri_tri_intersection( p1,q1,r1, d1,e1,f1 )),
          "TQGeom::tri_tri_intersection() failed." );
  ASSERT(!(TQGeom::tri_tri_intersection( p1,q1,r1, g1,h1,i1 )),
          "TQGeom::tri_tri_intersection() failed." );
  ASSERT(!(TQGeom::tri_tri_intersection( p1,q1,r1, p1,q1,r1 )),
          "TQGeom::tri_tri_intersection() failed." );

  /*------------------------------------------------------------------
  | Float
  ------------------------------------------------------------------*/
  const Vec2f p2 { 3.0f, 2.0f };
  const Vec2f q2 { 9.0f, 2.0f };
  const Vec2f r2 { 6.0f, 6.0f };

  const Vec2f s2 { 3.0f, 2.0f };
  const Vec2f t2 { 6.0f, 2.0f };
  const Vec2f u2 { 5.0f,-1.0f };

  const Vec2f v2 { 6.0f, 4.0f };
  const Vec2f w2 {11.0f, 4.0f };
  const Vec2f x2 {10.0f, 8.0f };

  const Vec2f a2 { 7.0f, 2.0f };
  const Vec2f b2 { 7.0f,-1.0f };
  const Vec2f c2 { 9.0f,-1.0f };

  const Vec2f d2 {-4.0f, 2.0f };
  const Vec2f e2 {-2.0f, 2.0f };
  const Vec2f f2 {-3.0f, 4.0f };

  const Vec2f g2 { 3.0f, 2.0f };
  const Vec2f h2 { 6.0f, 6.0f };
  const Vec2f i2 { 2.0f, 6.0f };

  ASSERT( (TQGeom::tri_tri_intersection( p2,q2,r2, s2,t2,u2 )),
          "TQGeom::tri_tri_intersection() failed." );
  ASSERT( (TQGeom::tri_tri_intersection( p2,q2,r2, v2,w2,x2 )),
          "TQGeom::tri_tri_intersection() failed." );
  ASSERT( (TQGeom::tri_tri_intersection( p2,q2,r2, a2,b2,c2 )),
          "TQGeom::tri_tri_intersection() failed." );
  ASSERT(!(TQGeom::tri_tri_intersection( p2,q2,r2, d2,e2,f2 )),
          "TQGeom::tri_tri_intersection() failed." );
  ASSERT(!(TQGeom::tri_tri_intersection( p2,q2,r2, g2,h2,i2 )),
          "TQGeom::tri_tri_intersection() failed." );
  ASSERT(!(TQGeom::tri_tri_intersection( p2,q2,r2, p2,q2,r2 )),
          "TQGeom::tri_tri_intersection() failed." );

  /*------------------------------------------------------------------
  | Integer
  ------------------------------------------------------------------*/
  const Vec2i p3 { 3, 2 };
  const Vec2i q3 { 9, 2 };
  const Vec2i r3 { 6, 6 };

  const Vec2i s3 { 3, 2 };
  const Vec2i t3 { 6, 2 };
  const Vec2i u3 { 5,-1 };

  const Vec2i v3 { 6, 4 };
  const Vec2i w3 {11, 4 };
  const Vec2i x3 {10, 8 };

  const Vec2i a3 { 7, 2 };
  const Vec2i b3 { 7,-1 };
  const Vec2i c3 { 9,-1 };

  const Vec2i d3 {-4, 2 };
  const Vec2i e3 {-2, 2 };
  const Vec2i f3 {-3, 4 };

  const Vec2i g3 { 3, 2 };
  const Vec2i h3 { 6, 6 };
  const Vec2i i3 { 2, 6 };

  ASSERT( (TQGeom::tri_tri_intersection( p3,q3,r3, s3,t3,u3 )),
          "TQGeom::tri_tri_intersection() failed." );
  ASSERT( (TQGeom::tri_tri_intersection( p3,q3,r3, v3,w3,x3 )),
          "TQGeom::tri_tri_intersection() failed." );
  ASSERT( (TQGeom::tri_tri_intersection( p3,q3,r3, a3,b3,c3 )),
          "TQGeom::tri_tri_intersection() failed." );
  ASSERT(!(TQGeom::tri_tri_intersection( p3,q3,r3, d3,e3,f3 )),
          "TQGeom::tri_tri_intersection() failed." );
  ASSERT(!(TQGeom::tri_tri_intersection( p3,q3,r3, g3,h3,i3 )),
          "TQGeom::tri_tri_intersection() failed." );
  ASSERT(!(TQGeom::tri_tri_intersection( p3,q3,r3, p3,q3,r3 )),
          "TQGeom::tri_tri_intersection() failed." );


  DBG_MSG("Tests for TQGeom::tri_tri_intersection() succeeded");

} // Test_TQGeom_tri_tri_intersection()

/*********************************************************************
* Test TQGeom::quad_quad_intersection
*********************************************************************/
void Test_TQGeom_quad_quad_intersection()
{
  /*------------------------------------------------------------------
  | Double
  ------------------------------------------------------------------*/
  const Vec2d p1d {-0.5,-0.5 };
  const Vec2d q1d { 0.5,-0.5 };
  const Vec2d r1d { 0.5, 0.5 };
  const Vec2d s1d {-0.5, 0.5 };

  // Regular intersecion
  const Vec2d p2d { 0.0,-1.0 };
  const Vec2d q2d { 1.0,-1.0 };
  const Vec2d r2d { 1.0, 0.0 };
  const Vec2d s2d { 0.0, 0.0 };

  // No intersection
  const Vec2d p3d {-2.5, 0.5 };
  const Vec2d q3d {-2.0, 0.0 };
  const Vec2d r3d {-1.5, 0.0 };
  const Vec2d s3d {-1.5, 1.0 };

  // Adjacent edge --> no intersection
  const Vec2d p4d {-0.5, 0.5 };
  const Vec2d q4d { 0.5, 0.5 };
  const Vec2d r4d { 0.5, 1.5 };
  const Vec2d s4d {-0.5, 1.5 };

  // Adjacent vertex --> no intersection
  const Vec2d p5d { 0.5, 0.5 };
  const Vec2d q5d { 1.5, 0.5 };
  const Vec2d r5d { 1.5, 1.5 };
  const Vec2d s5d { 0.5, 1.5 };

  // Half-Edge intersection
  const Vec2d p6d { 0.5,-0.5 };
  const Vec2d q6d { 1.0,-0.5 };
  const Vec2d r6d { 1.0, 0.0 };
  const Vec2d s6d { 0.5, 0.0 };

  ASSERT( (TQGeom::quad_quad_intersection(p1d, q1d, r1d, s1d, p2d, q2d, r2d, s2d)),
          "TQGeom::quad_quad_intersection() failed." );
  ASSERT(!(TQGeom::quad_quad_intersection(p1d, q1d, r1d, s1d, p3d, q3d, r3d, s3d)),
          "TQGeom::quad_quad_intersection() failed." );
  ASSERT(!(TQGeom::quad_quad_intersection(p1d, q1d, r1d, s1d, p4d, q4d, r4d, s4d)),
          "TQGeom::quad_quad_intersection() failed." );
  ASSERT(!(TQGeom::quad_quad_intersection(p1d, q1d, r1d, s1d, p5d, q5d, r5d, s5d)),
          "TQGeom::quad_quad_intersection() failed." );
  ASSERT( (TQGeom::quad_quad_intersection(p1d, q1d, r1d, s1d, p6d, q6d, r6d, s6d)),
          "TQGeom::quad_quad_intersection() failed." );

  /*------------------------------------------------------------------
  | Float
  ------------------------------------------------------------------*/
  const Vec2f p1f {-0.5f,-0.5f };
  const Vec2f q1f { 0.5f,-0.5f };
  const Vec2f r1f { 0.5f, 0.5f };
  const Vec2f s1f {-0.5f, 0.5f };

  // Regular intersecion
  const Vec2f p2f { 0.0f,-1.0f };
  const Vec2f q2f { 1.0f,-1.0f };
  const Vec2f r2f { 1.0f, 0.0f };
  const Vec2f s2f { 0.0f, 0.0f };

  // No intersection
  const Vec2f p3f {-2.5f, 0.5f };
  const Vec2f q3f {-2.0f, 0.0f };
  const Vec2f r3f {-1.5f, 0.0f };
  const Vec2f s3f {-1.5f, 1.0f };

  // Adjacent edge --> no intersection
  const Vec2f p4f {-0.5f, 0.5f };
  const Vec2f q4f { 0.5f, 0.5f };
  const Vec2f r4f { 0.5f, 1.5f };
  const Vec2f s4f {-0.5f, 1.5f };

  // Adjacent vertex --> no intersection
  const Vec2f p5f { 0.5f, 0.5f };
  const Vec2f q5f { 1.5f, 0.5f };
  const Vec2f r5f { 1.5f, 1.5f };
  const Vec2f s5f { 0.5f, 1.5f };

  // Half-Edge intersection
  const Vec2f p6f { 0.5f,-0.5f };
  const Vec2f q6f { 1.0f,-0.5f };
  const Vec2f r6f { 1.0f, 0.0f };
  const Vec2f s6f { 0.5f, 0.0f };

  ASSERT( (TQGeom::quad_quad_intersection(p1f, q1f, r1f, s1f, p2f, q2f, r2f, s2f)),
          "TQGeom::quad_quad_intersection() failed." );
  ASSERT(!(TQGeom::quad_quad_intersection(p1f, q1f, r1f, s1f, p3f, q3f, r3f, s3f)),
          "TQGeom::quad_quad_intersection() failed." );
  ASSERT(!(TQGeom::quad_quad_intersection(p1f, q1f, r1f, s1f, p4f, q4f, r4f, s4f)),
          "TQGeom::quad_quad_intersection() failed." );
  ASSERT(!(TQGeom::quad_quad_intersection(p1f, q1f, r1f, s1f, p5f, q5f, r5f, s5f)),
          "TQGeom::quad_quad_intersection() failed." );
  ASSERT( (TQGeom::quad_quad_intersection(p1f, q1f, r1f, s1f, p6f, q6f, r6f, s6f)),
          "TQGeom::quad_quad_intersection() failed." );

  /*------------------------------------------------------------------
  | Integer
  ------------------------------------------------------------------*/
  const Vec2i p1i {-5,-5 };
  const Vec2i q1i { 5,-5 };
  const Vec2i r1i { 5, 5 };
  const Vec2i s1i {-5, 5 };

  // Regular intersecion
  const Vec2i p2i { 0,-10 };
  const Vec2i q2i { 10,-10 };
  const Vec2i r2i { 10, 0 };
  const Vec2i s2i { 0, 0 };

  // No intersection
  const Vec2i p3i {-25, 5 };
  const Vec2i q3i {-20, 0 };
  const Vec2i r3i {-15, 0 };
  const Vec2i s3i {-15, 10 };

  // Adjacent edge --> no intersection
  const Vec2i p4i {-5, 5 };
  const Vec2i q4i { 5, 5 };
  const Vec2i r4i { 5, 15 };
  const Vec2i s4i {-5, 15 };

  // Adjacent vertex --> no intersection
  const Vec2i p5i { 5, 5 };
  const Vec2i q5i { 15, 5 };
  const Vec2i r5i { 15, 15 };
  const Vec2i s5i { 5, 15 };

  // Half-Edge intersection
  const Vec2i p6i { 5,-5 };
  const Vec2i q6i { 10,-5 };
  const Vec2i r6i { 10, 0 };
  const Vec2i s6i { 5, 0 };

  ASSERT( (TQGeom::quad_quad_intersection(p1i, q1i, r1i, s1i, p2i, q2i, r2i, s2i)),
          "TQGeom::quad_quad_intersection() failed." );
  ASSERT(!(TQGeom::quad_quad_intersection(p1i, q1i, r1i, s1i, p3i, q3i, r3i, s3i)),
          "TQGeom::quad_quad_intersection() failed." );
  ASSERT(!(TQGeom::quad_quad_intersection(p1i, q1i, r1i, s1i, p4i, q4i, r4i, s4i)),
          "TQGeom::quad_quad_intersection() failed." );
  ASSERT(!(TQGeom::quad_quad_intersection(p1i, q1i, r1i, s1i, p5i, q5i, r5i, s5i)),
          "TQGeom::quad_quad_intersection() failed." );
  ASSERT( (TQGeom::quad_quad_intersection(p1i, q1i, r1i, s1i, p6i, q6i, r6i, s6i)),
          "TQGeom::quad_quad_intersection() failed." );

  DBG_MSG("Tests for TQGeom::quad_quad_intersection() succeeded");

} // Test_TQGeom_quad_quad_intersection()

/*********************************************************************
* Test TQGeom::tri_quad_intersection
*********************************************************************/
void Test_TQGeom_tri_quad_intersection()
{
  /*------------------------------------------------------------------
  | Double
  ------------------------------------------------------------------*/
  const Vec2d pd {-1.0,-0.5 };
  const Vec2d qd { 0.5,-0.5 };
  const Vec2d rd { 0.5, 1.5 };
  const Vec2d sd {-1.0, 1.5 };

  // Regular intersection
  const Vec2d m1d {-2.0, 0.5 };
  const Vec2d n1d {-1.5,-1.5 };
  const Vec2d o1d {-0.5, 0.0 };

  // No intersection
  const Vec2d m2d { 1.5,-1.0 };
  const Vec2d n2d { 2.5,-1.0 };
  const Vec2d o2d { 2.0,-0.5 };

  // Adjacent edges --> no intersection
  const Vec2d m3d { 0.5,-0.5 };
  const Vec2d n3d { 1.5, 0.5 };
  const Vec2d o3d { 0.5,-1.0 };

  // Adjacent vertex --> no intersection
  const Vec2d m4d { 0.5, 1.5 };
  const Vec2d n4d { 1.5, 1.0 };
  const Vec2d o4d { 1.5, 1.5 };

  // Half edge intersection
  const Vec2d m5d { 0.5,-0.5 };
  const Vec2d n5d { 1.5,-0.5 };
  const Vec2d o5d { 0.5, 0.5 };

  ASSERT( (TQGeom::tri_quad_intersection(m1d, n1d, o1d, pd, qd, rd, sd)),
          "TQGeom::tri_quad_intersection() failed." );
  ASSERT(!(TQGeom::tri_quad_intersection(m2d, n2d, o2d, pd, qd, rd, sd)),
          "TQGeom::tri_quad_intersection() failed." );
  ASSERT(!(TQGeom::tri_quad_intersection(m3d, n3d, o3d, pd, qd, rd, sd)),
          "TQGeom::tri_quad_intersection() failed." );
  ASSERT(!(TQGeom::tri_quad_intersection(m4d, n4d, o4d, pd, qd, rd, sd)),
          "TQGeom::tri_quad_intersection() failed." );
  ASSERT( (TQGeom::tri_quad_intersection(m5d, n5d, o5d, pd, qd, rd, sd)),
          "TQGeom::tri_quad_intersection() failed." );

  /*------------------------------------------------------------------
  | Float
  ------------------------------------------------------------------*/
  const Vec2f pf {-1.0f,-0.5f };
  const Vec2f qf { 0.5f,-0.5f };
  const Vec2f rf { 0.5f, 1.5f };
  const Vec2f sf {-1.0f, 1.5f };

  // Regular intersection
  const Vec2f m1f {-2.0f, 0.5f };
  const Vec2f n1f {-1.5f,-1.5f };
  const Vec2f o1f {-0.5f, 0.0f };

  // No intersection
  const Vec2f m2f { 1.5f,-1.0f };
  const Vec2f n2f { 2.5f,-1.0f };
  const Vec2f o2f { 2.0f,-0.5f };

  // Adjacent edges --> no intersection
  const Vec2f m3f { 0.5f,-0.5f };
  const Vec2f n3f { 1.5f, 0.5f };
  const Vec2f o3f { 0.5f,-1.0f };

  // Adjacent vertex --> no intersection
  const Vec2f m4f { 0.5f, 1.5f };
  const Vec2f n4f { 1.5f, 1.0f };
  const Vec2f o4f { 1.5f, 1.5f };

  // Half edge intersection
  const Vec2f m5f { 0.5f,-0.5f };
  const Vec2f n5f { 1.5f,-0.5f };
  const Vec2f o5f { 0.5f, 0.5f };

  ASSERT( (TQGeom::tri_quad_intersection(m1f, n1f, o1f, pf, qf, rf, sf)),
          "TQGeom::tri_quad_intersection() failed." );
  ASSERT(!(TQGeom::tri_quad_intersection(m2f, n2f, o2f, pf, qf, rf, sf)),
          "TQGeom::tri_quad_intersection() failed." );
  ASSERT(!(TQGeom::tri_quad_intersection(m3f, n3f, o3f, pf, qf, rf, sf)),
          "TQGeom::tri_quad_intersection() failed." );
  ASSERT(!(TQGeom::tri_quad_intersection(m4f, n4f, o4f, pf, qf, rf, sf)),
          "TQGeom::tri_quad_intersection() failed." );
  ASSERT( (TQGeom::tri_quad_intersection(m5f, n5f, o5f, pf, qf, rf, sf)),
          "TQGeom::tri_quad_intersection() failed." );

  /*------------------------------------------------------------------
  | Integer
  ------------------------------------------------------------------*/
  const Vec2i pi {-10,-5 };
  const Vec2i qi { 5,-5 };
  const Vec2i ri { 5, 15 };
  const Vec2i si {-10, 15 };

  // Regular intersection
  const Vec2i m1i {-20, 5 };
  const Vec2i n1i {-15,-15 };
  const Vec2i o1i {-5, 0 };

  // No intersection
  const Vec2i m2i { 15,-10 };
  const Vec2i n2i { 25,-10 };
  const Vec2i o2i { 20,-5 };

  // Adjacent edges --> no intersection
  const Vec2i m3i { 5,-5 };
  const Vec2i n3i { 15, 5 };
  const Vec2i o3i { 5,-10 };

  // Adjacent vertex --> no intersection
  const Vec2i m4i { 5, 15 };
  const Vec2i n4i { 15, 10 };
  const Vec2i o4i { 15, 15 };

  // Half edge intersection
  const Vec2i m5i { 5,-5 };
  const Vec2i n5i { 15,-5 };
  const Vec2i o5i { 5, 5 };

  ASSERT( (TQGeom::tri_quad_intersection(m1i, n1i, o1i, pi, qi, ri, si)),
          "TQGeom::tri_quad_intersection() failed." );
  ASSERT(!(TQGeom::tri_quad_intersection(m2i, n2i, o2i, pi, qi, ri, si)),
          "TQGeom::tri_quad_intersection() failed." );
  ASSERT(!(TQGeom::tri_quad_intersection(m3i, n3i, o3i, pi, qi, ri, si)),
          "TQGeom::tri_quad_intersection() failed." );
  ASSERT(!(TQGeom::tri_quad_intersection(m4i, n4i, o4i, pi, qi, ri, si)),
          "TQGeom::tri_quad_intersection() failed." );
  ASSERT( (TQGeom::tri_quad_intersection(m5i, n5i, o5i, pi, qi, ri, si)),
          "TQGeom::tri_quad_intersection() failed." );

  DBG_MSG("Tests for TQGeom::tri_quad_intersection() succeeded");

} // Test_TQGeom_tri_quad_intersection()


/*********************************************************************
* Test TQGeom::rect_overlap()
*********************************************************************/
void Test_TQGeom_rect_overlap()
{
  /*------------------------------------------------------------------
  | Double
  ------------------------------------------------------------------*/
  const Vec2d al_d { -0.2, -0.2 };
  const Vec2d au_d {  0.2,  0.2 };

  const Vec2d bl_d {  0.0,  0.0 };
  const Vec2d bu_d {  0.2,  0.2 };

  const Vec2d cl_d {  0.1,  0.1 };
  const Vec2d cu_d {  0.1,  0.4 };

  const Vec2d dl_d {  0.2, -0.2 };
  const Vec2d du_d {  0.4,  0.3 };

  const Vec2d el_d {  0.2,  0.2 };
  const Vec2d eu_d {  0.4,  0.4 };

  const Vec2d fl_d {  0.3,  0.3 };
  const Vec2d fu_d {  0.4,  0.4 };

  ASSERT( (TQGeom::rect_overlap(al_d, au_d, bl_d, bu_d)),
          "TQGeom::rect_overlap() failed." );
  ASSERT( (TQGeom::rect_overlap(al_d, au_d, cl_d, cu_d)),
          "TQGeom::rect_overlap() failed." );
  ASSERT( (TQGeom::rect_overlap(al_d, au_d, dl_d, du_d)),
          "TQGeom::rect_overlap() failed." );
  ASSERT( (TQGeom::rect_overlap(al_d, au_d, el_d, eu_d)),
          "TQGeom::rect_overlap() failed." );
  ASSERT(!(TQGeom::rect_overlap(al_d, au_d, fl_d, fu_d)),
          "TQGeom::rect_overlap() failed." );

  DBG_MSG("Tests for TQGeom::rect_overlap() succeeded");

} // Test_TQGeom_rect_overlap()


/*********************************************************************
* Test TQGeom::in_on_rect()
*********************************************************************/
void Test_TQGeom_in_on_rect()
{
  /*------------------------------------------------------------------
  | Double
  ------------------------------------------------------------------*/
  const Vec2d al_d {-0.2, -0.2};
  const Vec2d au_d { 0.2,  0.2};
  
  const Vec2d p_d { 0.0,  0.0};
  const Vec2d q_d { 0.2,  0.0};
  const Vec2d r_d { 0.2,  0.2};
  const Vec2d s_d { 0.3,  0.3};


  ASSERT( (TQGeom::in_on_rect(p_d, al_d, au_d)),
          "TQGeom::in_on_rect() failed." );
  ASSERT( (TQGeom::in_on_rect(q_d, al_d, au_d)),
          "TQGeom::in_on_rect() failed." );
  ASSERT( (TQGeom::in_on_rect(r_d, al_d, au_d)),
          "TQGeom::in_on_rect() failed." );
  ASSERT(!(TQGeom::in_on_rect(s_d, al_d, au_d)),
          "TQGeom::in_on_rect() failed." );


  DBG_MSG("Tests for TQGeom::in_on_rect() succeeded");

} // Test_TQGeom_in_on_rect()


/*********************************************************************
* Test TQGeom::in_rect()
*********************************************************************/
void Test_TQGeom_in_rect()
{
  /*------------------------------------------------------------------
  | Double
  ------------------------------------------------------------------*/
  const Vec2d al_d {-0.2, -0.2};
  const Vec2d au_d { 0.2,  0.2};
  
  const Vec2d p_d { 0.0,  0.0};
  const Vec2d q_d { 0.2,  0.0};
  const Vec2d r_d { 0.2,  0.2};
  const Vec2d s_d { 0.3,  0.3};


  ASSERT( (TQGeom::in_rect(p_d, al_d, au_d)),
          "TQGeom::in_rect() failed." );
  ASSERT(!(TQGeom::in_rect(q_d, al_d, au_d)),
          "TQGeom::in_rect() failed." );
  ASSERT(!(TQGeom::in_rect(r_d, al_d, au_d)),
          "TQGeom::in_rect() failed." );
  ASSERT(!(TQGeom::in_rect(s_d, al_d, au_d)),
          "TQGeom::in_rect() failed." );


  DBG_MSG("Tests for TQGeom::in_rect() succeeded");

} // Test_TQGeom_in_rect()

/*********************************************************************
* Test TQGeom::in_on_triangle()
*********************************************************************/
void Test_TQGeom_in_on_triangle()
{
  /*------------------------------------------------------------------
  | Double
  ------------------------------------------------------------------*/
  const Vec2d t1_d {-0.2, -0.2};
  const Vec2d t2_d { 0.2, -0.2};
  const Vec2d t3_d { 0.0,  0.2};

  const Vec2d a_d { 0.0 , 0.0};
  const Vec2d b_d { 0.0 ,-0.2};
  const Vec2d c_d {-0.2 ,-0.2};
  const Vec2d d_d {-0.2 , 0.2};

  ASSERT( (TQGeom::in_on_triangle(a_d, t1_d, t2_d, t3_d)),
          "TQGeom::in_on_triangle() failed." );
  ASSERT( (TQGeom::in_on_triangle(b_d, t1_d, t2_d, t3_d)),
          "TQGeom::in_on_triangle() failed." );
  ASSERT( (TQGeom::in_on_triangle(c_d, t1_d, t2_d, t3_d)),
          "TQGeom::in_on_triangle() failed." );
  ASSERT(!(TQGeom::in_on_triangle(d_d, t1_d, t2_d, t3_d)),
          "TQGeom::in_on_triangle() failed." );


  DBG_MSG("Tests for TQGeom::in_on_triangle() succeeded");

} // Test_TQGeom_in_on_triangle()

/*********************************************************************
* Test TQGeom::in_triangle()
*********************************************************************/
void Test_TQGeom_in_triangle()
{
  /*------------------------------------------------------------------
  | Double
  ------------------------------------------------------------------*/
  const Vec2d t1_d {-0.2, -0.2};
  const Vec2d t2_d { 0.2, -0.2};
  const Vec2d t3_d { 0.0,  0.2};

  const Vec2d a_d { 0.0 , 0.0};
  const Vec2d b_d { 0.0 ,-0.2};
  const Vec2d c_d {-0.2 ,-0.2};
  const Vec2d d_d {-0.2 , 0.2};

  ASSERT( (TQGeom::in_triangle(a_d, t1_d, t2_d, t3_d)),
          "TQGeom::in_triangle() failed." );
  ASSERT(!(TQGeom::in_triangle(b_d, t1_d, t2_d, t3_d)),
          "TQGeom::in_triangle() failed." );
  ASSERT(!(TQGeom::in_triangle(c_d, t1_d, t2_d, t3_d)),
          "TQGeom::in_triangle() failed." );
  ASSERT(!(TQGeom::in_triangle(d_d, t1_d, t2_d, t3_d)),
          "TQGeom::in_triangle() failed." );


  DBG_MSG("Tests for TQGeom::in_triangle() succeeded");

} // Test_TQGeom_in_triangle()

/*********************************************************************
* Test TQGeom::in_on_quad()
*********************************************************************/
void Test_TQGeom_in_on_quad()
{
  /*------------------------------------------------------------------
  | Double
  ------------------------------------------------------------------*/
  const Vec2d q1_d { -0.2, -0.2 };
  const Vec2d q2_d {  0.2, -0.2 };
  const Vec2d q3_d {  0.2,  0.2 };
  const Vec2d q4_d { -0.2,  0.2 };

  const Vec2d a_d { 0.0, 0.0 };
  const Vec2d b_d {-0.2, 0.0 };
  const Vec2d c_d {-0.2, 0.0 };
  const Vec2d d_d {-0.3, 0.3 };

  ASSERT( (TQGeom::in_on_quad(a_d, q1_d, q2_d, q3_d, q4_d)),
          "TQGeom::in_on_quad() failed." );
  ASSERT( (TQGeom::in_on_quad(b_d, q1_d, q2_d, q3_d, q4_d)),
          "TQGeom::in_on_quad() failed." );
  ASSERT( (TQGeom::in_on_quad(c_d, q1_d, q2_d, q3_d, q4_d)),
          "TQGeom::in_on_quad() failed." );
  ASSERT(!(TQGeom::in_on_quad(d_d, q1_d, q2_d, q3_d, q4_d)),
          "TQGeom::in_on_quad() failed." );

  DBG_MSG("Tests for TQGeom::in_on_quad() succeeded");

} // Test_TQGeom_in_on_quad()

/*********************************************************************
* Test TQGeom::in_quad()
*********************************************************************/
void Test_TQGeom_in_quad()
{
  /*------------------------------------------------------------------
  | Double
  ------------------------------------------------------------------*/
  const Vec2d q1_d { -0.2, -0.2 };
  const Vec2d q2_d {  0.2, -0.2 };
  const Vec2d q3_d {  0.2,  0.2 };
  const Vec2d q4_d { -0.2,  0.2 };

  const Vec2d a_d { 0.0, 0.0 };
  const Vec2d b_d {-0.2, 0.0 };
  const Vec2d c_d {-0.2, 0.0 };
  const Vec2d d_d {-0.3, 0.3 };

  ASSERT( (TQGeom::in_quad(a_d, q1_d, q2_d, q3_d, q4_d)),
          "TQGeom::in_quad() failed." );
  ASSERT(!(TQGeom::in_quad(b_d, q1_d, q2_d, q3_d, q4_d)),
          "TQGeom::in_quad() failed." );
  ASSERT(!(TQGeom::in_quad(c_d, q1_d, q2_d, q3_d, q4_d)),
          "TQGeom::in_quad() failed." );
  ASSERT(!(TQGeom::in_quad(d_d, q1_d, q2_d, q3_d, q4_d)),
          "TQGeom::in_quad() failed." );

  DBG_MSG("Tests for TQGeom::in_quad() succeeded");

} // Test_TQGeom_in_quad()

} // Namespace GeometryTests

/*********************************************************************
* Run geometry tests
*********************************************************************/
void run_geometry_tests()
{
  MSG("\n#===== Geometry tests =====");

  GeometryTests::Test_TQGeom_orientation();
  GeometryTests::Test_TQGeom_is_left();
  GeometryTests::Test_TQGeom_is_lefton();
  GeometryTests::Test_TQGeom_in_segment();
  GeometryTests::Test_TQGeom_in_on_segment();
  GeometryTests::Test_TQGeom_line_line_intersection();
  GeometryTests::Test_TQGeom_line_tri_intersection();
  GeometryTests::Test_TQGeom_line_quad_intersection();
  GeometryTests::Test_TQGeom_tri_tri_intersection();
  GeometryTests::Test_TQGeom_quad_quad_intersection();
  GeometryTests::Test_TQGeom_tri_quad_intersection();
  GeometryTests::Test_TQGeom_rect_overlap();
  GeometryTests::Test_TQGeom_in_on_rect();
  GeometryTests::Test_TQGeom_in_rect();
  GeometryTests::Test_TQGeom_in_on_triangle();
  GeometryTests::Test_TQGeom_in_triangle();
  GeometryTests::Test_TQGeom_in_on_quad();
  GeometryTests::Test_TQGeom_in_quad();

}
