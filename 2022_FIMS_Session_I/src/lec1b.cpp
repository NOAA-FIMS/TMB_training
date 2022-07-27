#include <iostream>
#include <string>

int main() {
    double x = 1.1;
    int y = x;
    std::cout << "x = " << x << "; y = " << y << std::endl;
    std::string a = "a";
    std::string b = std::to_string(y);
    std::cout << a + b;
    return 0;
}
