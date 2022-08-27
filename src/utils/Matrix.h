/*
* This file is part of the CppUtils library.  
* This code was written by Florian Setzwein in 2022, 
* and is covered under the MIT License
* Refer to the accompanying documentation for details
* on usage and license.
*/
#pragma once

#include <vector>

namespace CppUtils {


/*********************************************************************
* Inspired by:
* -----------
*  https://stackoverflow.com/questions/20873768/more-efficient-way-\
*  with-multi-dimensional-arrays-matrices-using-c
*********************************************************************/
template<typename T, typename Allocator = std::allocator<T>>
class Matrix
{
public:
  /*------------------------------------------------------------------
  | Constructor
  ------------------------------------------------------------------*/
  Matrix() 
  : rows_ { 0 }
  , cols_ { 0 }
  {}

  Matrix(int r, int c) 
  : rows_ { r }
  , cols_ { c }
  , data_ (r*c, 0) 
  { }

  Matrix(T* data, int r, int c)
  : rows_ { r }
  , cols_ { c }
  , data_ (r*c, 0) 
  {
    std::copy(&data[0], &data[0] + r*c, const_cast<T*>(data_.data()));
  }

  Matrix(T** data, int r, int c)
  : rows_ { r }
  , cols_ { c }
  , data_ (r*c, 0) 
  {
    std::copy(data[0], data[0] + r*c, const_cast<T*>(data_.data()));
  }

  virtual ~Matrix() {}

  /*------------------------------------------------------------------
  | Copy / Move
  ------------------------------------------------------------------*/
  Matrix(const Matrix& m)
  : rows_ { m.rows_ }
  , cols_ { m.cols_ }
  , data_ { m.data_ }
  { }

  Matrix(Matrix&& m)
  : rows_ { m.rows_ }
  , cols_ { m.cols_ }
  , data_ { std::move(m.data_) }
  {}

  /*------------------------------------------------------------------
  | Operators
  ------------------------------------------------------------------*/
  inline Matrix& operator = (const Matrix& m)
  {
    rows_ = m.rows_;
    cols_ = m.cols_;
    data_.assign( m.data_.begin(), m.data_.end() );
    return *this;
  }

  inline T* operator[](const int i)
  { return data_.data() + i*cols_; }

  inline const T* operator [](const int i) const
  { return data_.data() + i*cols_; }

  /*------------------------------------------------------------------
  | Swap data 
  ------------------------------------------------------------------*/
  inline Matrix& swap(Matrix& m)
  {
    rows_ = m.rows_;
    cols_ = m.cols_;
    data_.swap( m.data_ );
    return *this;
  }

  /*------------------------------------------------------------------
  | Resize data
  ------------------------------------------------------------------*/
  inline void resize(int r, int c)
  { 
    rows_ = r;
    cols_ = c;
    data_.resize( r * c ); 
  }

  /*------------------------------------------------------------------
  | Getters
  ------------------------------------------------------------------*/
  inline int rows() { return rows_; }
  inline int rows() const { return rows_; }

  inline int columns() { return cols_; }
  inline int columns() const { return cols_; }

  inline std::size_t size() { return data_.size(); }
  inline std::size_t size() const { return data_.size(); }

private:
  /*------------------------------------------------------------------
  | Attributes
  ------------------------------------------------------------------*/
  int rows_ { 0 };
  int cols_ { 0 };

  std::vector<T, Allocator> data_;

};

} // namespace CppUtils
