#include <numcpp/functions.hpp>
#include "utils.hpp"
#include <iostream>

void testHypergeometric2F1() {

    assert(isClose(numcpp::functions::hypergeometric2F1(2.0,3.0,4.0,0.5),2.728935333122,1e-10));
}

int main() {

    testHypergeometric2F1();
    std::cout << "All tests for the module functions have been passed successfully!" << std::endl;
    return 0;
}