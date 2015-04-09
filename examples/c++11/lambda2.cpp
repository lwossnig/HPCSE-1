#include <iostream>
#include <thread>

int main()
{
  // create a function and store a pointer to it in f
  auto f = []() {std::cout << "Hello world!\n";};
  
  // call the function in a thread
  std::thread t(f);
  t.join();
}