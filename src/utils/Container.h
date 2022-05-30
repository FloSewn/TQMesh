/*
* This file is part of the CppUtils library.  
* This code was written by Florian Setzwein in 2022, 
* and is covered under the MIT License
* Refer to the accompanying documentation for details
* on usage and license.
*/
#pragma once

#include <list>           // std::list
#include <vector>         // std::vector
#include <memory>         // std::unique_ptr
#include <utility>        // std::move

#include "QuadTree.h"
#include "Vec2.h"
#include "Helpers.h"

namespace CppUtils {

constexpr double ContainerQuadTreeScale = 10000.0;
constexpr size_t ContainerQuadTreeItems = 100;
constexpr size_t ContainerQuadTreeDepth = 25;
static SimpleLogger ContainerLogger { std::clog, "[Container]: " };

/*********************************************************************
* This class is a container for two-dimensional objects that also 
* keeps track of its objects using a quadtree
*********************************************************************/
template <typename T>
class Container
{
public:
  using List           = std::list<std::unique_ptr<T>>;
  using WasteVector    = std::vector<std::unique_ptr<T>>;
  using size_type      = typename List::size_type;
  using value_type     = std::unique_ptr<T>;
  using Vector         = std::vector<T*>;
  using iterator       = typename List::iterator;
  using const_iterator = typename List::const_iterator;

  iterator begin() { return items_.begin(); }
  iterator end() { return items_.end(); }

  const_iterator begin() const { return items_.begin(); }
  const_iterator end() const { return items_.end(); }

  /*----------------------------------------------------------------------------
  | Constructor
  ----------------------------------------------------------------------------*/
  Container(double        qtree_scale = ContainerQuadTreeScale,
            size_t        qtree_items = ContainerQuadTreeItems, 
            size_t        qtree_depth = ContainerQuadTreeDepth,
            SimpleLogger& logger      = ContainerLogger)
  : qtree_ { qtree_scale, qtree_items, qtree_depth }  
  , logger_ { &logger }
  {}


  /*------------------------------------------------------------------
  | Copy constructor
  ------------------------------------------------------------------*/
  Container(const Container<T>& c) 
  : qtree_ { c.qtree().scale(),
             c.qtree().max_items(),
             c.qtree().max_depth() }  
  , logger_ { c.logger_ }
  {
    // Copy the list items
    for (auto& obj : c)
    {
      const Vec2d& xy = obj->xy();
      push_back( xy.x, xy.y, obj->sizing() );
    }

    // Copy the waste items
    for (auto& obj : c.waste())
      waste_.push_back( obj );
  }

  /*------------------------------------------------------------------
  | Move constructor
  ------------------------------------------------------------------*/
  Container(Container<T>&& c) 
  : qtree_ { c.qtree().scale(),
             c.qtree().max_items(),
             c.qtree().max_depth() }  
  , logger_ { c.logger_ }
  {
    items_ = std::move(c.items_);
    qtree_ = std::move(c.qtree_);
    waste_ = std::move(c.waste_);
  }

  /*------------------------------------------------------------------
  | Get a reference to the waste elements
  ------------------------------------------------------------------*/
  List& waste() { return waste_; }
  const List& waste() const { return waste_; }

  /*------------------------------------------------------------------
  | Get the current size of the list
  ------------------------------------------------------------------*/
  size_type size() const { return items_.size(); }

  /*------------------------------------------------------------------
  | Get reference to the  container qtreee
  ------------------------------------------------------------------*/
  const QuadTree<T,double>& qtree() const { return qtree_; }

  /*------------------------------------------------------------------
  | Get all items in a specified rectangle
  ------------------------------------------------------------------*/
  Vector get_items(const Vec2d& lowleft, 
                   const Vec2d& upright) const
  { 
    Vector found; 
    qtree_.get_items( lowleft, upright, found ); 
    return std::move( found );
  }

  /*------------------------------------------------------------------
  | Get all items in a specified radius
  ------------------------------------------------------------------*/
  Vector get_items(const Vec2d& center, 
                   const double radius) const
  { 
    Vector found; 
    qtree_.get_items( center, radius, found ); 
    return std::move( found );
  }

  /*------------------------------------------------------------------
  | Insert any simplex through constructor before 
  | a specified position
  ------------------------------------------------------------------*/
  template <typename... Args>
  T& insert( auto pos, Args&&... args )
  {
    std::unique_ptr<T> u_ptr = std::make_unique<T>(args...);
    T* ptr = u_ptr.get();
    iterator iter = items_.insert( pos, std::move(u_ptr) );
    ptr->pos_          = iter;
    ptr->in_container_ = true;
    bool in_qtree = qtree_.add( ptr );

    // Failed to add element to qtree -> cleanup
    if (!in_qtree)
    {
      *logger_ << "Failed to add element to the quad tree. "
                "Maybe the element is outside of the defined domain."
               << std::endl;
    }
    
    return *ptr;
  }

  /*------------------------------------------------------------------
  | Add a simplex through constructor to the end of 
  | the list
  ------------------------------------------------------------------*/
  template <typename... Args>
  T& push_back( Args&&... args )
  { return insert( items_.end(), args... ); }

  /*------------------------------------------------------------------
  | Function to remove an object from the container item list
  ------------------------------------------------------------------*/
  bool remove(T& item) 
  { 
    if ( qtree_.remove( &item ) )
    {
      // Move item to waste list
      waste_.splice(waste_.end(), items_, item.pos_);
      // Call the item container destructor
      item.container_destructor();
      // Mark element to be in the waste list
      item.in_container_ = false;

      return true;
    }
    return false;
  }

  /*------------------------------------------------------------------
  | Function to finally remove all items in the container garbage 
  | collector 
  ------------------------------------------------------------------*/
  void clear_waste() { waste_.clear(); }

  /*------------------------------------------------------------------
  | Access operator
  ------------------------------------------------------------------*/
  const T& operator[](size_t i) const
  {
    ASSERT( (i >= 0 && i < items_.size()),
            "Invalid access to Container" );
    auto iter = items_.begin();
    std::advance( iter, i );

    T& out = *(iter->get());
    return out;
  }

  T& operator[](size_t i)
  {
    ASSERT( (i >= 0 && i < items_.size()),
            "Invalid access to Container" );
    auto iter = items_.begin();
    std::advance( iter, i );
    T& out = *(iter->get());
    return out;
  }

  /*------------------------------------------------------------------
  | Return a reference to the last element in container
  | Calling std::list.back() function on an empty container 
  | causes undefined behavior.
  ------------------------------------------------------------------*/
  T& back() 
  { 
    ASSERT( (items_.size() > 0), 
        "Invalid access to container." );
    return *items_.back(); 
  } 
  const T& back() const 
  { 
    ASSERT( (items_.size() > 0), 
        "Invalid access to container." );
    return *items_.back(); 
  }

  /*------------------------------------------------------------------
  | Return a reference to the first element in container
  ------------------------------------------------------------------*/
  T& front() 
  { 
    ASSERT( (items_.size() > 0), 
        "Invalid access to container." );
    return *items_.front(); 
  } 
  const T& front() const 
  { 
    ASSERT( (items_.size() > 0), 
        "Invalid access to container." );
    return *items_.front(); 
  }

  /*------------------------------------------------------------------
  | Sort the elements in the container
  ------------------------------------------------------------------*/
  template <class Compare>
  void sort(Compare comp)
  { items_.sort( comp ); }



private:

  /*------------------------------------------------------------------
  | Container attributes
  ------------------------------------------------------------------*/
  List               items_;
  QuadTree<T,double> qtree_;
  List               waste_;
  SimpleLogger*      logger_;


}; // Container


} // namespace CppUtils
