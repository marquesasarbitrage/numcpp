#include <numcpp/optim/neldermead.hpp>
#include <iostream>
#include <cassert>

double testFunctionNM(const std::vector<double>& x) {
    double a = 1.0, b = 100.0;
    return pow(a - x[0], 2) + b * pow(x[1] - x[0] * x[0], 2);
}

double himmelblauFunction(const std::vector<double>& x) {
    double a = x[0] * x[0] + x[1] - 11;
    double b = x[0] + x[1] * x[1] - 7;
    return a * a + b * b;
}

void testNelderMead() {

    std::vector<double> x0 = {-1.2, 1.0};  

    numcpp::optim::NelderMeadResult result = numcpp::optim::nelderMead(x0, testFunctionNM,500,1e-6,{.05,1.0,2.0,0.5,0.75});

    assert(std::abs(result.x[0] - 1.0) < 1e-3);
    assert(std::abs(result.x[1] - 1.0) < 1e-3);
}

void testNelderMeadHimmelblau() {

    std::vector<double> x0 = {-5.0, 5.0};  

    numcpp::optim::NelderMeadResult result = numcpp::optim::nelderMead(x0, himmelblauFunction,1000,1e-6,{.05,1.0,2.0,0.5,0.5});

    bool near_minimum = (
        (std::abs(result.x[0] - 3.0) < 1e-3 && std::abs(result.x[1] - 2.0) < 1e-3) ||
        (std::abs(result.x[0] + 2.805) < 1e-3 && std::abs(result.x[1] - 3.131) < 1e-3) ||
        (std::abs(result.x[0] + 3.779) < 1e-3 && std::abs(result.x[1] + 3.283) < 1e-3) ||
        (std::abs(result.x[0] - 3.584) < 1e-3 && std::abs(result.x[1] + 1.848) < 1e-3)
    );

    assert(near_minimum);
}

int main()
{
    testNelderMead();
    testNelderMeadHimmelblau();
    std::cout << "All tests successfully passed for the nelder mead optimizer!" << std::endl;
    return 0;
}