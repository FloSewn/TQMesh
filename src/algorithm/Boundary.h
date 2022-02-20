/*
* This source file is part of the tqmesh library.  
* This code was written by Florian Setzwein in 2022, 
* and is covered under the MIT License
* Refer to the accompanying documentation for details
* on usage and license.
*/
#pragma once

#include "Vec2.h"
#include "utils.h"
#include "geometry.h"

#include "EdgeList.h"
#include "Edge.h"
#include "Vertex.h"

namespace TQMesh {
namespace TQAlgorithm {

using namespace TQUtils;

/********************************************************************* 
* Boundary type orientation
*********************************************************************/
enum class BdryType
{
  EXTERIOR,
  INTERIOR 
}; 


/*********************************************************************
* A mesh boundary - defined by a list of edges
* > Interior boundaries are defined clockwise (CW)
* > Exterior boundaries are defined counter-clockwise (CCW)
*********************************************************************/
class Boundary : public EdgeList
{
public:

  /*------------------------------------------------------------------
  | Constructor
  ------------------------------------------------------------------*/
  Boundary(BdryType btype)
  : EdgeList( (btype == BdryType::EXTERIOR 
              ? TQGeom::Orientation::CCW : TQGeom::Orientation::CW) ) 
  , btype_ { btype }
  { }

  /*------------------------------------------------------------------
  | Getters
  ------------------------------------------------------------------*/
  bool is_exterior() const 
  { return (btype_ == BdryType::EXTERIOR); }
  bool is_interior() const
  { return (btype_ == BdryType::INTERIOR); }

private:

  BdryType btype_;

}; // Boundary

} // namespace TQAlgorithm
} // namespace TQMesh
