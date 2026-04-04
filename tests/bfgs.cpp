#include <iostream>
#include <cassert>
#include <numcpp/objects/base.hpp>
#include <numcpp/optim/bfgs.hpp>

double f(const numcpp::objects::Vector& v) {
    double x = v(0), y = v(1);
    return std::pow(1.0 - x, 2) + 100.0 * std::pow(y - x * x, 2);
};

numcpp::objects::Vector g(const numcpp::objects::Vector& v) {
    double x = v(0), y = v(1);
    numcpp::objects::Vector grad(2);
    grad(0) = -2.0 * (1.0 - x) - 400.0 * x * (y - x * x);
    grad(1) =  200.0 * (y - x * x);
    return grad;
};

void testBFGSRosenbrock() {

    numcpp::objects::Vector x0(2);
    x0 << -1.5, 0.5;
 
    auto result = numcpp::optim::bfgs(f, g, x0);
    assert(std::abs(result.x(0)-1.0)<=1e-6);
    assert(std::abs(result.x(1)-1.0)<=1e-6);
    assert(std::abs(result.f)<=1e-10);


    auto result2 = numcpp::optim::bfgsForwardDifference(f, x0);
    assert(std::abs(result2.x(0)-1.0)<=1e-5);
    assert(std::abs(result2.x(1)-1.0)<=1e-4);
    assert(std::abs(result2.f)<=1e-10);

}

int main() {

    testBFGSRosenbrock(); 
    std::cout << "All tests for the BFGS optimization method have been passed!" << std::endl;
    return 0;

}