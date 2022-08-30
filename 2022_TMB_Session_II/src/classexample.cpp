#include <iostream>
#include <cmath>

template <class T>
class MyClass{
    public:
    T a;
    T b;

    MyClass(){}

    T evaluate(){
        return a + b;
    }

};
int main() {
  //initiate class
  MyClass<double> *m = new MyClass<double>();
  m -> a = 1.2;
  m -> b = 3.4;

  std::cout << m -> evaluate() << std::endl;  

  //intiate class
  MyClass<double> myclass;
  myclass.a = 2;
  myclass.b = 3;

  std::cout << myclass.evaluate() << std::endl;

  return 0;
}

