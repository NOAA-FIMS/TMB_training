#include <iostream>

int main() {
  //initiate a variable
  float a = 3.1459;
  //initiate a new variable with the same value as a
  float b = a;
  //initiate a variable c that points to the same address as a
  float* c = &a;
  *c = 100;
  std::cout << a << "=" << b << " " << c << std::endl;
  c = &b;
  std::cout << a << "=" << b << " " << c;
  return 0;
}

