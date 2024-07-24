#include <iostream>

#if __has_include(<generator>) && __cplusplus >= 202302
#define GENERATOR_AVAILABLE 1
#include <generator>
#else
#define GENERATOR_AVAILABLE 0
#endif



int main() {
#if GENERATOR_AVAILABLE
  std::cout << "generator available" << std::endl;
#else
  std::cout << "generator not available" << std::endl;
#endif
}

