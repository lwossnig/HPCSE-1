/// simpson.hpp
#include <iostream>
#include <functional>
#include <thread>
#include <mutex>
#include <cmath>
#include <vector>

#ifndef SIMPSON_HPP
#define SIMPSON_HPP

void simpson(const double, const double, const unsigned, std::function<double(double)>, std::pair<double, std::mutex>&);

#endif
