/*
* This file is part of the CppUtils library.  
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

#include "VecND.h"
#include "Geometry.h"
#include "Helpers.h"
#include "Log.h"

namespace CppUtils {


/*********************************************************************
* Default query function for QuadTree item search within rectangle
*********************************************************************/
template <typename T, typename V>
static inline bool quadtree_rect_query_fun(T* item, 
                                           const Vec2<V>& lowleft, 
                                           const Vec2<V>& upright)
{ return in_on_rect(item->xy(), lowleft, upright); }

/*********************************************************************
* Default query function for QuadTree item search within circle
*********************************************************************/
template <typename T, typename V>
static inline bool quadtree_circ_query_fun(T* item, 
                                           const Vec2<V>& c, 
                                           const V r_sqr)
{ return ( (item->xy() - c).norm_sqr() < r_sqr ); }

/*********************************************************************
* Default query function for QuadTree nearest neighbor search 
*********************************************************************/
template <typename T, typename V>
static inline bool quadtree_nearest_query_fun(T* item, 
                                              const Vec2<V>& c, 
                                              V& min_d_sqr)
{
  const V d_sqr = (item->xy() - c).norm_sqr();
  if ( d_sqr < min_d_sqr )
  {
    min_d_sqr = d_sqr;
    return true;
  }
  return false;
}

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
class QuadTree
{
public:
  using List   = std::list<T*>;
  using Array  = std::array<QuadTree<T,V>*,4>; 
  using Vector = std::vector<T*>;
  using ListIterator = typename List::iterator;

  /*------------------------------------------------------------------ 
  | Constructor
  ------------------------------------------------------------------*/
  QuadTree(V               scale, 
           size_t          max_item,
           size_t          max_depth, 
           const Vec2<V>&  center={0.0,0.0}, 
           size_t          depth=0,
           QuadTree<T,V>*  parent=nullptr)
  { 
    update_attributes( scale, max_item, max_depth, center ); 
    depth_  = depth;
    parent_ = parent;
  }

  /*------------------------------------------------------------------ 
  | Destructor
  ------------------------------------------------------------------*/
  ~QuadTree() 
  {
    if ( split_ )
    {
      for ( auto child : children_ )
      {
        ASSERT( child, "QuadTree structure is corrupted." );
        delete  child;
      }
    }
  } // ~QuadTree() 

  /*------------------------------------------------------------------ 
  | Getters
  ------------------------------------------------------------------*/
  size_t size() const { return n_items_; }
  bool split() const { return split_; }
  const List& items() const { return items_; }
  const Array& children() const { return children_; }
  V scale() const { return scale_; }
  size_t max_items() const { return max_item_; }
  size_t max_depth() const { return max_depth_; }
  const Vec2<V>& center() const { return center_; }
  const Vec2<V>& lowleft() const { return lowleft_; }
  const Vec2<V>& upright() const { return upright_; }

  /*------------------------------------------------------------------ 
  | Setters
  ------------------------------------------------------------------*/
  void scale(double v)  
  { update_attributes(v, max_item_, max_depth_, center_); }

  void max_item(size_t v)  
  { update_attributes(scale_, v, max_depth_, center_); }

  void max_depth(size_t v)  
  { update_attributes(scale_, max_item_, v, center_); }

  void center(const Vec2<V>& v)  
  { update_attributes(scale_, max_item_, max_depth_, v); }

  /*------------------------------------------------------------------ 
  | Query function pointers
  ------------------------------------------------------------------*/
  typedef bool (*RectQuery)(T* item, 
                            const Vec2<V>& lowleft, 
                            const Vec2<V>& upright);
  typedef bool (*CircQuery)(T* item, 
                            const Vec2<V>& center, 
                            const V radius_squared);
  typedef bool (*NearestQuery)(T* item, 
                               const Vec2<V>& center, 
                               V& min_dist_squared);

  /*------------------------------------------------------------------
  | Return the total number of qtree leafs
  ------------------------------------------------------------------*/
  int n_leafs( int n = 0 ) const 
  {
    if ( split_ )
    {
      for ( auto child : children_ )
      {
        ASSERT( child, "QuadTree structure is corrupted." );
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
  | Get the single nearest item to a given query point
  ------------------------------------------------------------------*/
  T* get_nearest(const Vec2<V>& query,
                 NearestQuery qfun = &(quadtree_nearest_query_fun)) const
  {
    auto quad = get_leaf(query);

    V min_dist_sqr = 4 * quad->scale() * quad->scale();
    T* winner      = nullptr;

    for ( auto item : quad->items() )
      if ( qfun(item, query, min_dist_sqr ) )
        winner = item;

    Vector found {};
    
    get_items(query, 4 * sqrt(min_dist_sqr), found);

    for ( auto item : found )
      if ( qfun(item, query, min_dist_sqr ) )
        winner = item;

    return winner;
  }

  /*------------------------------------------------------------------ 
  | Get leaf quad that encloses a given query point
  ------------------------------------------------------------------*/
  const QuadTree* get_leaf(const Vec2<V>& query) const
  {
    if ( !in_on_rect(query, lowleft_, upright_) )
      return nullptr;

    if ( split_ )
    {
      for ( auto child : children_ )
      {
        const QuadTree* quad = child->get_leaf(query);
        if (quad)
          return quad;
      }
    }
    return this;
  }



  /*------------------------------------------------------------------ 
  | Get items within bounding box
  ------------------------------------------------------------------*/
  size_t get_items(const Vec2<V>& lowleft, 
                   const Vec2<V>& upright,
                   Vector& found,
                   RectQuery qfun = &(quadtree_rect_query_fun)) const
  {
    bool overlap = rect_overlap(lowleft_, upright_,  
                                lowleft,  upright  );
                                 
    if ( !overlap ) return false;

    size_t n_found = 0;

    if ( split_ )
    {
      for ( auto child : children_ )
      {
        ASSERT( child, "QuadTree structure is corrupted." );
        n_found += child->get_items(lowleft,upright,found);
      }
    }
    else
    {
      for ( auto item : items_ )
        if ( qfun(item, lowleft, upright) ) 
        {
          found.push_back( item );
          ++n_found;
        }
    }

    return n_found;

  } // get_items()

  /*------------------------------------------------------------------ 
  | Get items within circle 
  ------------------------------------------------------------------*/
  size_t get_items(const Vec2<V>& center, 
                   V radius,
                   Vector& found,
                   CircQuery qfun = &(quadtree_circ_query_fun)) const
  {
    const V radius_scaled = 1.4142 * radius;
    const Vec2<V> lowleft = center - radius_scaled;
    const Vec2<V> upright = center + radius_scaled;

    bool overlap = rect_overlap( lowleft_, upright_,  
                                 lowleft,  upright  );
                                 
    if ( !overlap ) return false;

    size_t n_found = 0;

    if ( split_ )
    {
      for ( auto child : children_ )
      {
        ASSERT( child, "QuadTree structure is corrupted." );
        n_found += child->get_items(center, radius, found);
      }
    }
    else
    {
      const V radius_squared = radius * radius;

      for ( auto item : items_ )
      {
        if ( qfun(item, center, radius_squared) )
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
    if (!item || !in_on_rect(item->xy(), lowleft_, upright_))
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
        LOG(ERROR) << e.what(); 
      }
    }

    return true;

  } /* add() */

  /*------------------------------------------------------------------ 
  | Remove an item from the tree
  ------------------------------------------------------------------*/
  bool remove(T* item)
  {
    if (!item || !in_on_rect(item->xy(), lowleft_, upright_) )
      return false;

    if ( split_ )
    {
      bool child_removed = 
         (  children_[0]->remove(item) || children_[1]->remove(item)
         || children_[2]->remove(item) || children_[3]->remove(item) );

      ASSERT( child_removed, "Failed to remove item from QuadTree children.");
      (void) child_removed;

      if (n_items_ <= max_item_) 
      {
        bool merged = merge();
        ASSERT( merged, "Failed to merge QuadTree structure.");
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
          "Failed to remove item from QuadTree." );
      add_item_number( -1 );

      (void) orig_size;
      (void) new_size;
    }

    return true;

  } /* remove() */

protected:

  /*------------------------------------------------------------------ 
  | Initialize the quadtree geometry
  ------------------------------------------------------------------*/
  void update_attributes(double scale, size_t max_item, 
                         size_t max_depth, const Vec2<V>& center)
  {
    if ( split_ || parent_ || n_items_ > 0 )
      TERMINATE( "QuadTree::update_attributes(): "
        "Failed to update splitted quadtree.");

    ASSERT( scale > 0.0, "Invalid quad tree scale" );

    scale_       = scale;
    max_item_    = max_item;
    max_depth_   = max_depth;
    center_      = center;

    halfscale_   = scale_ / 2;
    Vec2<V> hsv  = { halfscale_, halfscale_ };
    lowleft_     = center_ - hsv;
    upright_     = center_ + hsv;
  }

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
  | Merge QuadTree children
  ------------------------------------------------------------------*/
  bool merge()
  {
    ASSERT( children_[0], "QuadTree structure is corrupted." );
    ASSERT( children_[1], "QuadTree structure is corrupted." );
    ASSERT( children_[2], "QuadTree structure is corrupted." );
    ASSERT( children_[3], "QuadTree structure is corrupted." );

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
    ASSERT( children_[0], "QuadTree structure is corrupted." );
    ASSERT( children_[1], "QuadTree structure is corrupted." );
    ASSERT( children_[2], "QuadTree structure is corrupted." );
    ASSERT( children_[3], "QuadTree structure is corrupted." );

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

    V h = halfscale_ / 2;

    // Centroids of children
    Vec2<V> c0 = { center_.x + h,  center_.y + h };
    Vec2<V> c1 = { center_.x - h,  center_.y + h };
    Vec2<V> c2 = { center_.x - h,  center_.y - h };
    Vec2<V> c3 = { center_.x + h,  center_.y - h };

    // Child quad: NORTH-EAST (NE)
    children_[0] = new QuadTree<T,V> 
    { halfscale_, max_item_, max_depth_, c0, child_depth, this};
                               
    // Child quad: NORTH-WEST (NW)
    children_[1] = new QuadTree<T,V> 
    { halfscale_, max_item_, max_depth_, c1, child_depth, this};

    // Child quad: SOUTH-WEST (SW)
    children_[2] = new QuadTree<T,V> 
    { halfscale_, max_item_, max_depth_, c2, child_depth, this};

    // Child quad: SOUTH-EAST (SE)
    children_[3] = new QuadTree<T,V> 
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
          "Failed to split QuadTree. "
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

      ASSERT( children_[0], "QuadTree structure is corrupted." );
      ASSERT( children_[1], "QuadTree structure is corrupted." );
      ASSERT( children_[2], "QuadTree structure is corrupted." );
      ASSERT( children_[3], "QuadTree structure is corrupted." );
      
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
  V              scale_     { 0.0 };
  V              halfscale_ { 0.0 };
  size_t         max_item_  { 0 };
  size_t         max_depth_ { 0 };

  Vec2<V>        center_    { 0.0, 0.0 };
  Vec2<V>        lowleft_   { 0.0, 0.0 };
  Vec2<V>        upright_   { 0.0, 0.0 };

  bool           split_     { false };
  size_t         depth_     { 0 };
  size_t         n_items_   { 0 };

  List           items_;
  Array          children_  { nullptr };
  QuadTree<T,V>* parent_    { nullptr };

}; // QuadTree


/*********************************************************************
* 
*********************************************************************/
template <typename T, typename V>
class QuadTreeBuilder {

public:

  // Enforce singleton
  QuadTreeBuilder<T,V>(const QuadTreeBuilder<T,V>&) = delete;
  QuadTreeBuilder<T,V>& operator=(const QuadTreeBuilder<T,V>&) = delete;

  // Set quadtree properties
  QuadTreeBuilder<T,V>& set_scale(V s) { scale_ = s; return *this; }
  QuadTreeBuilder<T,V>& set_max_items(size_t n) { max_items_ = n; return *this; }
  QuadTreeBuilder<T,V>& set_max_depth(size_t n){ max_depth_ = n; return *this; }
  QuadTreeBuilder<T,V>& set_center(const Vec2<V>& c) { center_ = c; return *this; }

  // Get singleton instance
  static QuadTreeBuilder<T,V>& get_instance() 
  {
    static QuadTreeBuilder<T,V> instance;
    return instance;
  }

  // Built a new quad tree
  QuadTree<T,V> build() const
  {
    return QuadTree<T,V> { scale_, max_items_, max_depth_, center_ };
  }


private:

  // Enforce singleton
  QuadTreeBuilder<T,V>() {}

  V       scale_     { V{} };
  size_t  max_items_ { 0 };
  size_t  max_depth_ { 0 };
  Vec2<V> center_    { V{}, V{} };

}; // QuadTreeBuilder


/*********************************************************************
* Stream to std::cout
*********************************************************************/
template<typename T, typename V>
std::ostream& operator<<(std::ostream& os, 
                         const QuadTree<T,V>& qt)
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


} // namespace CppUtils
