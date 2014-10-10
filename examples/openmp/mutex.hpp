#include <omp.h>

class mutex {
public:
  mutex() { omp_init_lock(&mutex_);}
  ~mutex() { omp_destroy_lock(&mutex_);}
  void lock() { omp_set_lock(&mutex_);}
  void unlock() { omp_unset_lock(&mutex_);}
private:
  omp_lock_t mutex_;

};

template <class Mutex>
class lock_guard {
public:
  lock_guard(Mutex& m) 
  : mutex_(m)
  { 
    mutex_.lock();
  }

  ~lock_guard()
  {
    mutex_.unlock();
  }

private:
  Mutex& mutex_;
};
