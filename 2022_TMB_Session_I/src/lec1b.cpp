#include <iostream>
#include <string>

int main() {
    int x = 1;
    std::string a = "a";
    std::string b = std::to_string(x);
    std::cout << a + b;
    return 0;
}
