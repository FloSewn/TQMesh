/*
* This source file is part of the tqmesh library.  
* This code was written by Florian Setzwein in 2022, 
* and is covered under the MIT License
* Refer to the accompanying documentation for details
* on usage and license.
*/
#pragma once

#include <iostream>
#include <string>

#include "MathUtility.h"

namespace CppUtils {

/*--------------------------------------------------------------------
| A simple command line progress bar
--------------------------------------------------------------------*/
class ProgressBar
{
public:

  ProgressBar(int w = 70)
  : width_ { w }
  {}

  /*------------------------------------------------------------------
  | 
  ------------------------------------------------------------------*/
  void update(int prog) 
  { 
    int p = MIN(100, MAX(0, prog));

    if ( p <= progress_ ) 
    {
      update_ = false;
      return;
    }

    old_progress_ = progress_;
    progress_     = p;
    update_       = true;
  }

  /*------------------------------------------------------------------
  | 
  ------------------------------------------------------------------*/
  void show(std::ostream& os) const
  {
    if ( !update_ ) return; 

    os << "# [";
    int pos = width_ * progress_ / 100;

    for (int i = 0; i < width_; ++i) {
        if (i < pos) os << "=";
        else if (i == pos) os << ">";
        else os << " ";
    }

    os << "] " << progress_ << " %\r";
    os.flush();
  }

private:
  int   width_;
  int   progress_     {0};
  int   old_progress_ {0};
  bool  update_       {false};

}; // ProgressBar


/*********************************************************************
* 
*********************************************************************/
static inline std::ostream& operator<<(std::ostream& os, 
                                       const ProgressBar& pbar)
{
  pbar.show( os );
  return os;
}

} // namespace CppUtils
