
#include <iostream>
#include "dvector.hpp"

int main() {

  dvector<int> v;

  int x = 1;
  std::cout << "add " << x << " to position " << v.add(x) << std::endl;
  std::cout << v << std::endl;

  x++;
  std::cout << "add " << x << " to position " << v.add(x) << std::endl;
  std::cout << v << std::endl;

  x++;  
  std::cout << "add " << x << " to position " << v.add(x) << std::endl;
  std::cout << v << std::endl;
  
  x++;
  std::cout << "add " << x << " to position " << v.add(x) << std::endl;
  std::cout << v << std::endl;

  std::cout << "remove at index " << 2 << std::endl;
  v.remove(2);
  
  x++;
  std::cout << "add " << x << " to position " << v.add(x) << std::endl;
  std::cout << v << std::endl;
  
  
  
}

