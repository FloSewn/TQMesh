/*
* This source file is part of the tqmesh library.  
* This code was written by Florian Setzwein in 2022, 
* and is covered under the MIT License
* Refer to the accompanying documentation for details
* on usage and license.
*/
#pragma once

#include "TQMesh.h"

namespace TQMesh {

using namespace CppUtils;

/*********************************************************************
* A size function definition by the user
*********************************************************************/
using UserSizeFunction = std::function<double(const Vec2d& xy)>;


/*********************************************************************
* This class defines the local mesh size
*********************************************************************/
class SizeFunction
{
public:
  /*------------------------------------------------------------------
  | Constructor
  ------------------------------------------------------------------*/
  SizeFunction(UserSizeFunction f) : f_ { f } {}

  /*------------------------------------------------------------------
  | Getter
  ------------------------------------------------------------------*/
  const UserSizeFunction& user_size_function() const { return f_; }
  UserSizeFunction& user_size_function() { return f_; }

  /*------------------------------------------------------------------
  | Evaluate the domain's size function at a given point
  ------------------------------------------------------------------*/
  template <typename Domain>
  inline double evaluate(const Vec2d& xy, 
                         const Domain& domain) const
  {
    // The underlying size function
    const double h_fun = f_(xy);

    if ( h_fun <= 0.0 )
      TERMINATE("SizeFunction::evaluate(): Encountered invalid value (<=0).");

    double h = h_fun;

    // Gather distance contribution of each boundary vertices
    for ( const auto& boundary : domain )
      for ( const auto& edge : boundary.get()->edges() )
      {
        const Vec2d&  v1_xy = edge->v1().xy();
        const Vec2d&  v2_xy = edge->v2().xy();

        const double el = edge->length();
        const double r = MAX(h_fun/el, el/h_fun);

        const double d1_sqr = (xy - v1_xy).norm_sqr();
        const double d2_sqr = (xy - v2_xy).norm_sqr();

        const double s1 = (edge->v1().size_range() <= 0.0) 
                        ? el : edge->v1().size_range();
        const double s2 = (edge->v2().size_range() <= 0.0) 
                        ? el : edge->v2().size_range();

        const double s1_inv = 1.0 / (r* s1);
        const double s2_inv = 1.0 / (r* s2);

        const double z1 = exp(-d1_sqr * s1_inv * s1_inv);
        const double z2 = exp(-d2_sqr * s2_inv * s2_inv);

        const double h1 = (edge->v1().mesh_size() <= 0.0)
                        ? h_fun : edge->v1().mesh_size();

        const double h2 = (edge->v2().mesh_size() <= 0.0)
                        ? h_fun : edge->v2().mesh_size();

        const double hv1 = z1*MIN(h1,el) + (1.0-z1)*h_fun;
        const double hv2 = z2*MIN(h2,el) + (1.0-z2)*h_fun;

        h = MIN(h, hv1);
        h = MIN(h, hv2);
      }

    // Gather distance contribution of fixed vertices
    for ( auto& vertex : domain.fixed_vertices() )
    {
      if (vertex->mesh_size() <= 0.0) 
        continue;

      const double d_sqr = (xy - vertex->xy()).norm_sqr();
      const double s_inv = (vertex->size_range() <= 0.0) 
                         ? 1.0/h_fun : 1.0/vertex->size_range();

      const double z = exp(-d_sqr * s_inv * s_inv);
      h = MIN(h, z*vertex->mesh_size() + (1.0-z)*h_fun);
    }

    return h;

  } // evaluate()

  /*------------------------------------------------------------------
  | Export the size function to a cartesian grid
  ------------------------------------------------------------------*/
  template<typename Domain>
  void export_size_function(std::ostream& os,
                            const Domain& domain,
                            const Vec2d& xy_min, const Vec2d& xy_max,
                            unsigned int Nx, unsigned int Ny) const
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
        const double h = this->evaluate( xy, domain );

        values.push_back( h );
      }
    }

    os << "SIZE-FUNCTION " 
       << std::setprecision(5) << std::fixed 
       << xy_min.x << " " << xy_min.y << " "
       << xy_max.x << " " << xy_max.y << " "
       << Nx << " " << Ny << "\n";

    // Print data to the output stream
    unsigned int nx = 10;
    unsigned int ny = (Nx*Ny) / nx;
    unsigned int index = 0;

    for ( unsigned int j = 0; j < ny; ++j )
    {
      for ( unsigned int i = 0; i < nx; ++i )
      {
        os << std::setprecision(5) << std::fixed 
           << values[index]  << ( (i==nx-1) ? "" : "," );
        ++index;
      }
      os << "\n";
    }

    for ( ; index < Nx * Ny; ++index )
      os << std::setprecision(5) << std::fixed 
                << values[index]  << ( (index==Nx*Ny-1) ? "" : "," );
    os << "\n";

  } // export_size_function()

private:

  /*------------------------------------------------------------------
  | Attributes
  ------------------------------------------------------------------*/
  UserSizeFunction f_;

}; // SizeFunction


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
          double qtree_scale = ContainerQuadTreeScale,
          size_t qtree_items = ContainerQuadTreeItems, 
          size_t qtree_depth = ContainerQuadTreeDepth )
  : size_fun_ { f }
  , verts_ { qtree_scale, qtree_items, qtree_depth }
  { }

  /*------------------------------------------------------------------
  | Get the current number of all boundaries in the domain
  ------------------------------------------------------------------*/
  size_type size() const { return boundaries_.size(); }

  /*------------------------------------------------------------------
  | Get the maximum x-/y-extent of all vertices of this domain
  ------------------------------------------------------------------*/
  std::pair<Vec2d,Vec2d> extent() const 
  {
    Vec2d x_extent { 0.0, 0.0 };
    Vec2d y_extent { 0.0, 0.0 };

    for ( const auto& v : verts_ )
    {
      x_extent[0] = MIN( x_extent[0], v->x() );
      x_extent[1] = MAX( x_extent[1], v->x() );
      
      y_extent[0] = MIN( y_extent[0], v->y() );
      y_extent[1] = MAX( y_extent[1], v->y() );
    }

    return {x_extent, y_extent};
  }

  /*------------------------------------------------------------------
  | Get edges within a given point and radius
  ------------------------------------------------------------------*/
  std::vector<Edge*> 
  get_edges (const Vec2d& center, double radius) const 
  {
    std::vector<Edge*> found {};

    for ( const auto& boundary : *this )
    {
      auto cur_found = boundary->get_edges(center, radius);
      found.insert(found.end(), cur_found.begin(), cur_found.end());
    }

    return std::move( found );
  }

  /*------------------------------------------------------------------
  | Getter
  ------------------------------------------------------------------*/
  const Vertices& vertices() const { return verts_; }
  Vertices& vertices() { return verts_; }

  const UserSizeFunction& user_size_function() const 
  { return size_fun_.user_size_function(); }
  UserSizeFunction& user_size_function() 
  { return size_fun_.user_size_function(); }

  const VertexVector& fixed_vertices() const { return fixed_verts_; }
  VertexVector& fixed_vertices() { return fixed_verts_; }

  /*------------------------------------------------------------------
  | Setter 
  ------------------------------------------------------------------*/
  void quad_tree_scale(double v) { verts_.quad_tree().scale(v); }
  void quad_tree_max_item(size_t v) { verts_.quad_tree().max_item(v); }
  void quad_tree_max_depth(size_t v) { verts_.quad_tree().max_depth(v); }
  void quad_tree_center(const Vec2d& v) { verts_.quad_tree().center(v); }

  /*------------------------------------------------------------------
  | Evaluate the domain's size function at a given point
  ------------------------------------------------------------------*/
  inline double size_function(const Vec2d& xy) const
  { return size_fun_.evaluate(xy, *this); }

  /*------------------------------------------------------------------
  | Insert any boundary through constructor behind 
  | a specified position
  ------------------------------------------------------------------*/
  template <typename... Args>
  Boundary& insert_boundary( const_iterator pos, Args&&... args )
  {
    std::unique_ptr<Boundary> b_ptr 
      = std::make_unique<Boundary>(this->verts_, args...);

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
  | Count the number of edge overlaps between this and another domain
  ------------------------------------------------------------------*/
  size_t count_edge_overlaps(const Domain& nbr_domain)
  {
    size_t n_overlaps = 0;

    for ( const auto& boundary : *this )
    {
      for ( const auto& e : boundary->edges() )
      {
        const Vec2d v1 = e->v1().xy();
        const Vec2d w1 = e->v2().xy();

        // search in vicinity of current edge for edges of the 
        // neighboring domain. 
        const Vec2d c  = e->xy();
        double radius  = edge_overlap_range_ * e->length();

        std::vector<Edge*> nbr_edges = nbr_domain.get_edges(c, radius);

        for ( Edge* e_nbr : nbr_edges )
        {
          const Vec2d v2 = e_nbr->v1().xy();
          const Vec2d w2 = e_nbr->v2().xy();

          if ( segment_overlap(v1,w1, v2,w2) )
            ++n_overlaps;
        }
      }
    }

    return n_overlaps;

  } // Domain::count_edge_overlaps()

  /*------------------------------------------------------------------
  | Add a vertex to the domain
  ------------------------------------------------------------------*/
  template <typename... Args>
  Vertex& add_vertex( Args&&... args )
  {
    Vertex& v_new = verts_.push_back( args... );

    return v_new;

  } // Domain::add_vertex() 

  /*------------------------------------------------------------------
  | Remove a vertex from the domain
  ------------------------------------------------------------------*/
  void remove_vertex(Vertex& v) 
  { verts_.remove( v ); } 

  /*------------------------------------------------------------------
  | Add a fixed vertex, which will be incorporated in the meshing
  | process
  ------------------------------------------------------------------*/
  template <typename... Args>
  Vertex& add_fixed_vertex( Args&&... args )
  {
    Vertex& v_new = verts_.push_back( args... );
    v_new.add_property(VertexProperty::is_fixed);

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
  void export_size_function(std::ostream& os,
                            const Vec2d& xy_min, const Vec2d& xy_max,
                            unsigned int Nx, unsigned int Ny)
  { size_fun_.export_size_function(os, *this, xy_min, xy_max, Nx, Ny); }

  /*------------------------------------------------------------------
  | This function calculates the area enclosed by all boundaries
  ------------------------------------------------------------------*/
  double area() const
  {
    double area = 0.0;
    for ( const auto& boundary : (*this) )
      area += boundary.get()->area();
    return area;
  } 


private:
  Vector           boundaries_;

  SizeFunction     size_fun_;
  Vertices         verts_;
  VertexVector     fixed_verts_ {};

  double           edge_overlap_range_ { 1.5 };

}; // Domain

} // namespace TQMesh
