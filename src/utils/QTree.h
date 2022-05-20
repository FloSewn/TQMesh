/*
* This source file is part of the tqmesh library.  
* This code was written by Florian Setzwein in 2022, 
* and is covered under the MIT License
* Refer to the accompanying documentation for details
* on usage and license.
*/
#pragma once

#include <list>      // std::list
#include <vector>    // std::vector
#include <array>     // std::array
#include <stdexcept> // std::runtime_error
#include <iomanip>   // std::setprecision
#include <iostream>  // std::to_string

#include "Vec2.h"
#include "geometry.h"

/*
 *
 * TO DO:
 * ------
 * * Replace "throw" statements with "ASSERT()" 
 * * QTree removal performes bad if size is very large
 *   --> scale should not be constant in container, but 
 *       a variable value that depends on the case
 *
*/  

namespace TQMesh {
namespace TQUtils {

/*********************************************************************
* This class refers to a quad tree for 2D simplices
*
*
*                  scale_
*         |<-------------------->| (upright_)
*         *----------------------*--
*         |                      | ^
*         |                      | |
*         |                      | |
*         |       (center_)      | |
*         |           *          | | scale_
*         |                      | |
*         |                      | |
*         |                      | |
*         |                      | v
*         *----------------------*--
*  (lowleft_)
*
*********************************************************************/
template <typename T, typename V>
class QTree
{
public:
  using List   = std::list<T*>;
  using Array  = std::array<QTree<T,V>*,4>; 
  using Vector = std::vector<T*>;
  using ListIterator = typename List::iterator;

  /*------------------------------------------------------------------ 
  | Constructor
  ------------------------------------------------------------------*/
  QTree(double          scale, 
        size_t          max_item,
        size_t          max_depth, 
        const Vec2<V>&  center={0.0,0.0}, 
        size_t          depth=0,
        QTree<T,V>*     parent=nullptr)
  : scale_      { scale     }
  , max_item_   { max_item  }
  , max_depth_  { max_depth }
  , center_     { center    }
  , depth_      { depth     }
  , parent_     { parent    }
  {  
    halfscale_ = 0.5 * scale_;
    Vec2<V> hsv  = { halfscale_, halfscale_ };
    lowleft_   = center_ - hsv;
    upright_   = center_ + hsv;

  } // QTree()

  /*------------------------------------------------------------------ 
  | Destructor
  ------------------------------------------------------------------*/
  ~QTree() 
  {
    if ( split_ )
    {
      for ( auto child : children_ )
      {
        ASSERT( child, "QTree structure is corrupted." );
        delete  child;
      }
    }

  } // ~QTree() 

  /*------------------------------------------------------------------ 
  | Getters
  ------------------------------------------------------------------*/
  size_t size() const { return n_items_; }
  bool split() const { return split_; }
  const List& items() const { return items_; }
  const Array& children() const { return children_; }
  double scale() const { return scale_; }
  size_t max_items() const { return max_item_; }
  size_t max_depth() const { return max_depth_; }
  const Vec2<V>& center() const { return center_; }
  const Vec2<V>& lowleft() const { return lowleft_; }
  const Vec2<V>& upright() const { return upright_; }

  /*------------------------------------------------------------------
  | Return the total number of qtree leafs
  ------------------------------------------------------------------*/
  int n_leafs( int n = 0 ) const 
  {
    if ( split_ )
    {
      for ( auto child : children_ )
      {
        ASSERT( child, "QTree structure is corrupted." );
        n = child->n_leafs( n );
      }
    }
    else
    {
      ++n;
    }

    return n;

  } /* n_leafs() */

  /*------------------------------------------------------------------ 
  | Get items within bounding box
  ------------------------------------------------------------------*/
  size_t get_items(const Vec2<V>& lowleft, 
                   const Vec2<V>& upright,
                   Vector& found) const
  {
    bool overlap = TQGeom::rect_overlap(lowleft_, upright_,  
                                        lowleft,  upright  );
                                 
    if ( !overlap ) return false;

    size_t n_found = 0;

    if ( split_ )
    {
      for ( auto child : children_ )
      {
        ASSERT( child, "QTree structure is corrupted." );
        n_found += child->get_items(lowleft,upright,found);
      }
    }
    else
    {
      for ( auto item : items_ )
        if ( TQGeom::in_on_rect(item->xy(), lowleft, upright) )
        {
          found.push_back( item );
          ++n_found;
        }
    }

    return n_found;

  } /* get_items() */

  /*------------------------------------------------------------------ 
  | Get items within circle 
  ------------------------------------------------------------------*/
  size_t get_items(const Vec2<V>& center, 
                   const double radius,
                   Vector& found) const
  {
    const double radius_scaled = 1.4142 * radius;
    const Vec2<V> lowleft = center - radius_scaled;
    const Vec2<V> upright = center + radius_scaled;

    bool overlap = TQGeom::rect_overlap( lowleft_, upright_,  
                                         lowleft,  upright  );
                                 
    if ( !overlap ) return false;

    size_t n_found = 0;

    if ( split_ )
    {
      for ( auto child : children_ )
      {
        ASSERT( child, "QTree structure is corrupted." );
        n_found += child->get_items(center, radius, found);
      }
    }
    else
    {
      const double radius_squared = radius * radius;

      for ( auto item : items_ )
      {
        const Vec2<V> dist = item->xy() - center;
        if ( dist.length_squared() < radius_squared )
        {
          found.push_back( item );
          ++n_found;
        }
      }
    }

    return n_found;

  } /* get_items() */

  /*------------------------------------------------------------------ 
  | Add a new item to the tree
  ------------------------------------------------------------------*/
  bool add(T* item)
  {
    // Check if item is located within this quad
    if (!item || !TQGeom::in_on_rect(item->xy(), lowleft_, upright_))
      return false;

    // If quad is splitted, pass item to children
    if ( split_ ) 
      return pass_children( item );

    // Add item to this quad
    items_.push_back( item );

    // Increase total number of items in entire tree
    add_item_number( 1 );

    // Split this quad if maxmium item number is reached
    if ( (items_.size() > max_item_) && (depth_ < max_depth_) )
    {
      try 
      { 
        split_qtree(); 
      } 
      catch ( std::exception const& e ) 
      {
        std::cout << "[QTree-Error]: " << e.what() << "\n"; 
      }
    }

    return true;

  } /* add() */

  /*------------------------------------------------------------------ 
  | Remove an item from the tree
  ------------------------------------------------------------------*/
  bool remove(T* item)
  {
    if (!item || !TQGeom::in_on_rect(item->xy(), lowleft_, upright_) )
      return false;

    if ( split_ )
    {
      bool child_removed = 
         (  children_[0]->remove(item) || children_[1]->remove(item)
         || children_[2]->remove(item) || children_[3]->remove(item) );

      ASSERT( child_removed, "Failed to remove item from QTree children.");
      (void) child_removed;

      if (n_items_ <= max_item_) 
      {
        bool merged = merge();
        ASSERT( merged, "Failed to merge QTree structure.");
        (void) merged;
        return true;
      }

    }
    else
    {
      size_t orig_size = items_.size();
      items_.remove( item );
      size_t new_size = items_.size();

      ASSERT( (new_size < orig_size), 
          "Failed to remove item from QTree." );
      add_item_number( -1 );

      (void) orig_size;
      (void) new_size;
    }

    return true;

  } /* remove() */

protected:

  /*------------------------------------------------------------------ 
  | Add to total number of items in entire tree 
  ------------------------------------------------------------------*/
  inline void add_item_number( int i ) 
  {
    n_items_ += i;

    if ( parent_ ) 
      parent_->add_item_number( i );

  } /* add_item_number() */

private:

  /*------------------------------------------------------------------ 
  | Merge QTree children
  ------------------------------------------------------------------*/
  bool merge()
  {
    ASSERT( children_[0], "QTree structure is corrupted." );
    ASSERT( children_[1], "QTree structure is corrupted." );
    ASSERT( children_[2], "QTree structure is corrupted." );
    ASSERT( children_[3], "QTree structure is corrupted." );

    if (  children_[0]->split() 
       || children_[1]->split() 
       || children_[2]->split() 
       || children_[3]->split() )
      return false;

    for ( auto child : children_ )
    {
      for ( auto item : child->items() )
        items_.push_back( item );
      
      delete child;
      child = nullptr;
    }

    split_ = false;

    return true;

  } /* merge() */

  /*------------------------------------------------------------------ 
  | Pass an item to children
  ------------------------------------------------------------------*/
  bool pass_children(T* item)
  {
    ASSERT( children_[0], "QTree structure is corrupted." );
    ASSERT( children_[1], "QTree structure is corrupted." );
    ASSERT( children_[2], "QTree structure is corrupted." );
    ASSERT( children_[3], "QTree structure is corrupted." );

    if (  children_[0]->add( item ) 
       || children_[1]->add( item ) 
       || children_[2]->add( item ) 
       || children_[3]->add( item ) )
      return true;

    return false;
  }

  /*------------------------------------------------------------------ 
  | Split this quad
  ------------------------------------------------------------------*/
  void split_qtree() 
  { 
    size_t child_depth = depth_ + 1;

    double h = 0.5 * halfscale_;

    // Centroids of children
    Vec2<V> c0 = { center_.x + h,  center_.y + h };
    Vec2<V> c1 = { center_.x - h,  center_.y + h };
    Vec2<V> c2 = { center_.x - h,  center_.y - h };
    Vec2<V> c3 = { center_.x + h,  center_.y - h };

    // Child quad: NORTH-EAST (NE)
    children_[0] = new QTree<T,V> 
    { halfscale_, max_item_, max_depth_, c0, child_depth, this};
                               
    // Child quad: NORTH-WEST (NW)
    children_[1] = new QTree<T,V> 
    { halfscale_, max_item_, max_depth_, c1, child_depth, this};

    // Child quad: SOUTH-WEST (SW)
    children_[2] = new QTree<T,V> 
    { halfscale_, max_item_, max_depth_, c2, child_depth, this};

    // Child quad: SOUTH-EAST (SE)
    children_[3] = new QTree<T,V> 
    { halfscale_, max_item_, max_depth_, c3, child_depth, this};

    // Distribute items among children
    try 
    { 
      distribute_items(); 
    }
    catch ( std::exception const& e ) 
    { 
      throw e; 
    }

    // Check that quad is empty
    if ( items_.size() > 0 )
      throw std::runtime_error(
          "Failed to split QTree. "
          "Data structure might be corrupted.");

    // Mark this quad as splitted
    split_ = true;
    
  } /* split_qtree() */

  /*------------------------------------------------------------------ 
  | Distribute items to children qtrees
  ------------------------------------------------------------------*/
  void distribute_items()
  {
    while ( items_.size() > 0 )
    {
      auto item = items_.back();

      ASSERT( children_[0], "QTree structure is corrupted." );
      ASSERT( children_[1], "QTree structure is corrupted." );
      ASSERT( children_[2], "QTree structure is corrupted." );
      ASSERT( children_[3], "QTree structure is corrupted." );
      
      if (  children_[0]->add(item)
         || children_[1]->add(item)
         || children_[2]->add(item)
         || children_[3]->add(item) )
        add_item_number( -1 );
      else
        throw std::runtime_error(
            "Failed to distribute items. "
            "Data structure might be corrupted.");

      items_.pop_back();
    }

    ASSERT( (items_.size() == 0),
            "Failed to empty list.");

  } /* distribute_items() */

  /*------------------------------------------------------------------
  | Attributes
  ------------------------------------------------------------------*/
  double       scale_     { 0.0 };
  double       halfscale_ { 0.0 };
  size_t       max_item_  { 0 };
  size_t       max_depth_ { 0 };

  Vec2<V>      center_   { 0.0, 0.0 };
  Vec2<V>      lowleft_  { 0.0, 0.0 };
  Vec2<V>      upright_  { 0.0, 0.0 };

  bool         split_     { false };
  size_t       depth_     { 0 };
  size_t       n_items_   { 0 };

  List         items_;
  Array        children_  { nullptr };
  QTree<T,V>*  parent_    { nullptr };



}; // QTree

/*********************************************************************
* Stream to std::cout
*********************************************************************/
template<typename T, typename V>
std::ostream& operator<<(std::ostream& os, 
                         const QTree<T,V>& qt)
{
  if ( qt.split() )
  {
    for ( auto child : qt.children() )
      if ( child ) os << *child;
  }
  else
  {
    os << std::setprecision(5) << std::fixed
       << qt.center().x << ","
       << qt.center().y << ","
       << qt.scale()    << ","
       << qt.size()     << "\n";
  }

  return os;
  
}


} // namespace TQUtils
} // namespace TQMesh 
