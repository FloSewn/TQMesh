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

namespace TQMesh {
namespace TQAlgorithm {

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

} // namespace TQAlgorithm
} // namespace TQMesh
