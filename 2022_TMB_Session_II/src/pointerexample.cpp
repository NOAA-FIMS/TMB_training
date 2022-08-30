#include <iostream>
#include <cmath>

template <class T>
  T my_add(T *x, T y){
    *x += y;
    return 0;
}

template <class T>
  T my_exp(T *x){
    *x = std::exp(*x);
    return 0;
}

int main() {
  //initiate a variable
  float par = 0;
  //initiate a new constant used to update par
  float b = 0.6931472;
  //intiate a pointer to par
  float *c = &par;
  
  //update par using functions
  my_add(c, b);
  std::cout << "par = " << par << std::endl;
  my_exp(c);
  std::cout << "par = " << par << std::endl;

  return 0;
}

