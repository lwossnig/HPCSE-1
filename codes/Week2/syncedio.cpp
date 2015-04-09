#include<iostream>
#include<vector>
#include<thread>

std::mutex io_mutex;

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
