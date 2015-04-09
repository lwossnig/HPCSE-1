/// sync_io.cpp
#include <iostream>
#include <thread>
#include <vector>

std::mutex io_mutex; // global mutex
struct sync
{
	sync(std::ostream& os)
		: os(os),
		lock(io_mutex) {}

	template <class T>
		std::ostream& operator<< (T const& x)
		{
			return os << x;
		}
	private:
	std::ostream& os;
	std::lock_guard<std::mutex> lock;
};

/// usage:
/// sync(std::cout) << ... 

