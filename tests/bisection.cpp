#include <numcpp/optim/bisection.hpp>
#include <numcpp/optim/newtonraphson.hpp>
#include <iostream>


double f(double x) { return x * x - 2.0; };
double df(double x) { return 2.0 * x; };

void testBisectionVersusNewton() {

    numcpp::optim::NewtonRaphsonResult nr = numcpp::optim::newtonRaphson(1.0,f,df); 
    numcpp::optim::BisectionResult bi = numcpp::optim::bisection(0.0, 2.0, f); 

    std::cout << "Newton result:\n";
    std::cout << "x = " << nr.x 
              << ", f(x) = " << nr.f
              << ", iter = " << nr.iter << "\n\n";

    std::cout << "Bisection result:\n";
    std::cout << "x = " << bi.x 
              << ", f(x) = " << bi.f
              << ", iter = " << bi.iter << "\n";


}

int main() {

    testBisectionVersusNewton();
    return 0; 
}