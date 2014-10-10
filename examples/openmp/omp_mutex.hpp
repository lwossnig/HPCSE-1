// Example codes for HPC course
// (c) 2012 Matthias Troyer, ETH Zurich

#include <omp.h>

class omp_mutex
{
public:
  omp_mutex() { omp_init_lock(&_mutex); }
  ~omp_mutex() { omp_destroy_lock(&_mutex); }
  void lock() { omp_set_lock(&_mutex); }
  void unlock() { omp_unset_lock(&_mutex); }
  omp_mutex(omp_mutex const&) = delete;
  omp_mutex& operator=(omp_mutex const&) = delete;
private:
  omp_lock_t _mutex;
};


class omp_recursive_mutex
{
public:
  omp_recursive_mutex() { omp_init_nest_lock(&_mutex); }
  ~omp_recursive_mutex() { omp_destroy_nest_lock(&_mutex); }
  void lock() { omp_set_nest_lock(&_mutex); }
  void unlock() { omp_unset_nest_lock(&_mutex); }
  omp_recursive_mutex(omp_recursive_mutex const&) = delete;
  omp_recursive_mutex& operator=(omp_recursive_mutex const&) = delete;
private:
  omp_nest_lock_t _mutex;
};


