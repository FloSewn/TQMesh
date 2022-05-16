/*
* This source file is part of the tqmesh library.  
* This code was written by Florian Setzwein in 2022, 
* and is covered under the MIT License
* Refer to the accompanying documentation for details
* on usage and license.
*/
#pragma once

#include <vector>         // std::vector
#include <memory>         // std::shared_ptr
#include <utility>        // std::move
#include <array>          // std::array
#include <functional>     // std::function

#include "utils.h"
#include "Boundary.h"

namespace TQMesh {
namespace TQAlgorithm {

using namespace TQUtils;

/*********************************************************************
* A size function definition by the user
*********************************************************************/
using UserSizeFunction = std::function<double(const Vec2d& xy)>;

/*********************************************************************
* This class is a simple container for boundaries
*********************************************************************/
class Domain
{
public:
  using Vector       = std::vector<std::unique_ptr<Boundary>>;
  using size_type    = typename Vector::size_type;
  using value_type   = Boundary;
  using VertexVector = std::vector<Vertex*>;

  using iterator       = typename Vector::iterator;
  using const_iterator = typename Vector::const_iterator;

  iterator begin() { return boundaries_.begin(); }
  iterator end() { return boundaries_.end(); }

  const_iterator begin() const { return boundaries_.begin(); }
  const_iterator end() const { return boundaries_.end(); }

  /*------------------------------------------------------------------
  | Constructor
  ------------------------------------------------------------------*/
  Domain( UserSizeFunction f = [](const Vec2d& p){return 1.0;},
          double qtree_scale = TQ_QTREE_SCALE,
          size_t qtree_items = TQ_QTREE_ITEMS, 
          size_t qtree_depth = TQ_QTREE_DEPTH,
          double min_size    = TQ_DOMAIN_MIN_ELEMSIZE,
          double min_scaling = TQ_DOMAIN_MIN_SCALING )
  : f_ { f }
  , verts_ { qtree_scale, qtree_items, qtree_depth }
  , min_size_ { min_size }
  , min_scaling_ { min_scaling }
  { }

  /*------------------------------------------------------------------
  | Get the current number of all boundaries in the domain
  ------------------------------------------------------------------*/
  size_type size() const { return boundaries_.size(); }

  /*------------------------------------------------------------------
  | Getter
  ------------------------------------------------------------------*/
  const Vertices& vertices() const { return verts_; }
  Vertices& vertices() { return verts_; }

  /*------------------------------------------------------------------
  | Insert any boundary through constructor behind 
  | a specified position
  ------------------------------------------------------------------*/
  template <typename... Args>
  Boundary& insert_boundary( auto pos, Args&&... args )
  {
    std::unique_ptr<Boundary> b_ptr 
      = std::make_unique<Boundary>(args...);

    Boundary* ptr = b_ptr.get();

    boundaries_.insert( pos, std::move(b_ptr) );

    return *ptr;
  }

  /*------------------------------------------------------------------
  | Add a boundary through constructor to the end of the vector
  ------------------------------------------------------------------*/
  template <typename... Args>
  Boundary& add_boundary( Args&&... args )
  { return insert_boundary( boundaries_.end(), args... ); }

  /*------------------------------------------------------------------
  | Add an exterior boundary 
  ------------------------------------------------------------------*/
  Boundary& add_exterior_boundary()
  { return add_boundary( BdryType::EXTERIOR ); }

  /*------------------------------------------------------------------
  | Add an interior boundary 
  ------------------------------------------------------------------*/
  Boundary& add_interior_boundary()
  { return add_boundary( BdryType::INTERIOR ); }

  /*------------------------------------------------------------------
  | Remove a boundary from the domain
  ------------------------------------------------------------------*/
  void remove_boundary(size_t pos) 
  { boundaries_.erase( boundaries_.begin()+pos ); }

  /*------------------------------------------------------------------
  | Access operator
  ------------------------------------------------------------------*/
  Boundary& operator[](size_t i) const 
  { 
    ASSERT( (i < boundaries_.size()),
            "Invalid access to Domain");
    return *(boundaries_[i].get());
  } 

  /*------------------------------------------------------------------
  | Check if an object is inside of all exterior boundaries
  | and outside of all interior boundaries
  ------------------------------------------------------------------*/
  template <typename T>
  bool is_inside(const T& s) const 
  {
    bool inside = true;

    for ( const auto& b_ptr : boundaries_ )
    {
      Boundary* b = b_ptr.get();
      inside &= (  ( b->is_exterior() && b->is_inside(s)   )
                || ( b->is_interior() && !b->is_inside(s) ) );
    }

    return inside;
  }

  /*------------------------------------------------------------------
  | Evaluate the domain's size function at a given point
  ------------------------------------------------------------------*/
  inline double size_function(const Vec2d& xy) const
  {
    double scaling = 1.0;

    const double h = f_(xy);

    // Reduces value by contribution of boundary vertices
    for ( const auto& boundary : (*this) )
    {
      for ( const auto& e : boundary.get()->edges() )
      {
        const Vertex& v1 = e->v1();

        // User defined vertex sizing
        const double s_v = v1.sizing();

        // At small boundary segments, reduce size function locally
        // in order to obtain smooth element distribution
        const double fxy   = scaling * f_( v1.xy() );
        const double ratio = sqrt( e->length() / fxy );
        const double s_e   = std::min(1.0, ratio );

        const double range = v1.range();
        const double l2 = range * range; 
        const double r2 = ( xy - v1.xy() ).length_squared();
        const double fac  = exp( -r2 / l2);
        scaling *= 1.0 + fac * (s_e*s_v - 1.0);
      }
    }

    // Reduces value by contribution of fixed vertices
    for ( auto& v : fixed_verts_ )
    {
      const double s_v = v->sizing();
      const double range = v->range();

      const double l2 = range * range; 
      const double r2 = ( xy - v->xy() ).length_squared();
      const double fac  = exp( -r2 / l2);
      scaling *= 1.0 + fac * (s_v - 1.0);
    }

    // Apply lower threshold for scale factor 
    scaling = std::max(scaling, min_scaling_);

    return ( std::max( h * scaling, min_size_ ) );
  }

  /*------------------------------------------------------------------
  | Add a vertex to the domain
  ------------------------------------------------------------------*/
  template <typename... Args>
  Vertex& add_vertex( Args&&... args )
  {
    Vertex& v_new = verts_.push_back( args... );
    v_new.is_fixed( false );
    v_new.on_front( true );
    v_new.on_boundary( false );

    return v_new;
  } // Domain::add_vertex() 

  /*------------------------------------------------------------------
  | Remove a vertex from the domain
  ------------------------------------------------------------------*/
  void remove_vertex(Vertex& v) 
  { 
    verts_.remove( v );

  } // Domain::remove_vertex()


  /*------------------------------------------------------------------
  | Add a fixed vertex, which will be incorporated in the meshing
  | process
  ------------------------------------------------------------------*/
  template <typename... Args>
  Vertex& add_fixed_vertex( Args&&... args )
  {
    Vertex& v_new = verts_.push_back( args... );
    v_new.is_fixed( true );
    v_new.on_front( true );
    v_new.on_boundary( false );

    fixed_verts_.push_back( &v_new );

    return v_new;

  } // Domain::add_fixed_vertex()

  /*------------------------------------------------------------------
  | Remove a fixed vertex
  ------------------------------------------------------------------*/
  void remove_fixed_vertex(Vertex& v) 
  { 
    auto itr = std::remove_if(
        fixed_verts_.begin(), fixed_verts_.end(), 
        [&](Vertex* a){return *a == v;}
    );

    fixed_verts_.erase( itr, fixed_verts_.end() );

    verts_.remove( v );

  } // Domain::remove_fixed_vertex()



  /*------------------------------------------------------------------
  | This function prints out the size function of the domain onto
  | a cartesian grid 
  ------------------------------------------------------------------*/
  void export_size_function(const Vec2d& xy_min, const Vec2d& xy_max,
                            unsigned int Nx, unsigned int Ny)
  {
    const Vec2d len = xy_max - xy_min;
    const Vec2d dxy = { len.x / static_cast<double>(Nx),
                        len.y / static_cast<double>(Ny) };

    std::vector<double> values {};
    values.reserve( Nx * Ny );

    // Compute size function values at various positions
    for (unsigned int j = 0; j < Ny; ++j)
    {
      for (unsigned int i = 0; i < Nx; ++i)
      {
        const Vec2d xy = { xy_min.x + static_cast<double>(i)*dxy.x,
                           xy_min.y + static_cast<double>(j)*dxy.y };
        const double h = this->size_function( xy );

        values.push_back( h );
      }
    }

    std::cout << "SIZE-FUNCTION " 
      << std::setprecision(5) << std::fixed 
      << xy_min.x << " " << xy_min.y << " "
      << xy_max.x << " " << xy_max.y << " "
      << Nx << " " << Ny << "\n";

    // Print data to the user
    unsigned int nx = 10;
    unsigned int ny = (Nx*Ny) / nx;
    unsigned int index = 0;

    for ( unsigned int j = 0; j < ny; ++j )
    {
      for ( unsigned int i = 0; i < nx; ++i )
      {
        std::cout << std::setprecision(5) << std::fixed 
                  << values[index]  << ( (i==nx-1) ? "" : "," );
        ++index;
      }
      std::cout << "\n";
    }

    for ( ; index < Nx * Ny; ++index )
      std::cout << std::setprecision(5) << std::fixed 
                << values[index]  << ( (index==Nx*Ny-1) ? "" : "," );
    std::cout << "\n";

  } // Domain::export_size_function()


  /*------------------------------------------------------------------
  | This function calculates the area enclosed by all boundaries
  ------------------------------------------------------------------*/
  double area() const
  {
    double area = 0.0;
    
    for ( const auto& boundary : (*this) )
      area += boundary.get()->area();

    return area;

  } // Domain::area()


private:
  Vector           boundaries_;

  UserSizeFunction f_;
  Vertices         verts_;
  double           min_size_    { 0.0 };
  double           min_scaling_ { 0.0 };

  VertexVector     fixed_verts_ {};

}; // Domain

} // namespace TQAlgorithm
} // namespace TQMesh
