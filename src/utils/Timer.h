/*
* This file is part of the CppUtils library.  
* This code was written by Florian Setzwein in 2022, 
* and is covered under the MIT License
* Refer to the accompanying documentation for details
* on usage and license.
*/
#pragma once

#include <chrono>
#include <vector>

namespace CppUtils {

/*********************************************************************
* A simple class to measure benchmark times
*********************************************************************/
class Timer
{
public:
  using Clock = std::chrono::high_resolution_clock;
  using Second = std::chrono::duration<double, std::ratio<1> >;
  using Timepoint = std::chrono::time_point<Clock>;
  using Timevector = std::vector<Timepoint>;
  using Msgvector = std::vector<std::string>;

  Timer() = default;

  // Getter
  const Timevector& times() const { return tv_; }
  const Msgvector& messages() const { return msg_; }
  const std::string& message(std::size_t i) const { return msg_[i]; }
  std::size_t size() const { return tv_.size(); }

  // Measure time
  void count(const std::string& msg="") 
  { 
    tv_.push_back( Clock::now() ); 
    msg_.push_back( msg );
  }

  // Get delta between timepoints in seconds
  double delta(int i) const
  {
    int sz = tv_.size();
    if ( i > (sz-2) || i < 0 )
      throw std::runtime_error("Invalid timer access");

    auto t1 = tv_[i+1];
    auto t0 = tv_[i];

    return std::chrono::duration_cast<Second>(t1-t0).count();
  }

private:

  Timevector tv_;
  Msgvector  msg_;


}; // Timer


} // namespace CppUtils
