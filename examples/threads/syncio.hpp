// syncoed ip  for HPC course
// (c) 2012 Matthias Troyer, ETH Zurich

#ifndef HPC12_SYNCIO_HPP
#define HPC12_SYNCIO_HPP

#include <thread>
#include <mutex>

std::mutex io_mutex;  // global

struct sync
{
  sync( std::ostream& os )
  : os(os)
  , lock(io_mutex) {}
  
  template <class T>
  std::ostream& operator<<(T const& x)
  {
    return os << x;
  }
  
private:
  std::ostream& os;
  std::lock_guard<std::mutex> lock;
};

#endif