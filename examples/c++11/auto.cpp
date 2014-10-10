#include <iostream>
#include <functional>

int f(int x) { return x+1;}
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
  std::cout << p3(42)    << std::endl;
}
