#include <iostream>

int main()
{
  // create a function and store a pointer to it in f
  auto f = []() {std::cout << "Hello world!\n";};
  
  // call the function
  f();
}