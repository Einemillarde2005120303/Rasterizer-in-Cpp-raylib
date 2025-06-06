#include <iostream>
#include <thread>

void function1(char symbol) {
    for (size_t i = 0; i < 500; i++) std::cout << symbol;
}

void function2(char symbol) {
    for (size_t i = 0; i < 500; i++) std::cout << symbol;
}

int main() {
    std::thread t1(function1, '#');
    std::thread t2(function2, '-');

    t1.join(); // Wait for t1 to finish
    t2.join(); // Wait for t2 to finish

    return 0;
}