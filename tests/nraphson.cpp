#include <numcpp/optim/newtonraphson.hpp>
#include <iostream>
#include <cassert>

double testFunctionNR(double x) {
    return x * x - 4;  
}

double testFunctionDerivativeNR(double x) {
    return 2 * x;
}

void testNewtonRaphson() {

    numcpp::optim::NewtonRaphsonResult result = numcpp::optim::newtonRaphson(1.0, testFunctionNR, testFunctionDerivativeNR);
    assert(std::abs(result.x - 2.0) < 1e-4 || std::abs(result.x+ 2.0) < 1e-4);

}

int main() {

    testNewtonRaphson(); 
    std::cout << "All tests successfully passed for the newton raphson optimizer!" << std::endl;
}

