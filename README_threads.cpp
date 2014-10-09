# To compile code using pthreads (c++11), requirements: g++-4.7 or higher
g++ -std=c++11 -pthread ...

# or alternatively use clang++ (v: 4.0)
clang++ -std=c++11 -stdlib=libc++

# launch thread by using the thread constructor
std::thread t (foo, arg1);

# or using c++11 lambda functions:
auto f = [] () {std::cout << "Hello world!\n";};

# call the funtion f in a thread
std::thread t(f);

# threads are joined by using:
t.join();

# How to move/ use threads
# can store them in vectors
std::vector<std::thread> v(10);

# Can pass to/return from functions
std::thread t = make_thread();
do_something(make_thread());

# using detach() the thread runs without refering to it
# but be CAREFUL! It can become a security hole or bug!
# Boost silently detaches a joinable (still running) thread
# C++11 terminates joinable threads





##########################################################
#		           C++11	  		 #
##########################################################

# C++11 random number generators and generator engines
std::mt19937 mt; // create Mersenne twister generator engine
# seed the generator
mt.seed(<arbitrary_number>); // or time(0)

# build a vector for the see sequence
int N = ...;
std::vector<int> seeds(N);
std::seed_seq seq(seeds.begin(), seeds.end());	//create sequence and fill the vector
mt.seed(seq);					//feed the generator with the sequence

# std::function  - C++11 functionals
#include<functional>				// Header for functions
std::function<Result(Arg1, Arg2, Arg3)>		// declaration of the result and argument types

# example:
double simpson(std::function<double(double)> f,  double a, double b, unsigned int N);

# new auto type
auto x = 3.14159 + 5;

###############################
# EXAMPLE FOR THE USAGE:
################################
#include<iostream>
#include<functional>

int f(int x){return x+1;}
int main()
{
    // function pointer
    int (*p1)(int) = f;

    // easier function pointer with auto
    auto p2 =f;

    // or here we could just have used std::function
    std::function<int(int)> p3=f;

    std::cout << (*p1)(42) << std::endl;
    std::cout << (*p2)(42) << std::endl;
    std::cout << p3(42)	   << std::endl;
}

# How to integrate a function with two parameters with a lambda function using simpson rule?
#include "simpson.hpp"
#include <iostream>
#include <cmath>

int main()
{
    double a = 3.4;

    //create a lambda function
    // [=] indicates that the variable a should be used inside the lambda
    std::cout << simpson([=] (double x) { return std::exp(a*x)},0.,1.,100) << std::endl;
    return 0;

}

# How to use the lambda function?
# [] - capture nothing
# [&] - capture any reference variable by reference
# [=] - capture any reference variable by making a copy
# [=, &foo] - capture any reference variable by making a copy, but capture variable foo by reference
# [bar] - capture bar by making a copy; don't copy anything else
# [this] - Capture the this pointer of the enclosing class


    


