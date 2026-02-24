#include "numcpp/probability.hpp"
#include <cassert>
#include <iostream>

bool isClose(double a, double b, double eps) { return std::abs(a-b)<eps;}

void testNormalDistributions() {

    using namespace numcpp::probability;
    Normal n{}; 

    assert(std::abs(cdf(n,0.0) - 0.5) < 1e-6); 
    assert(std::abs(cdf(n,1.0) - 0.841344746) < 1e-6);
    assert(std::abs(cdf(n,-1.0) - 0.158655254) < 1e-6);
    assert(std::abs(cdf(n,2.0) - 0.977249868) < 1e-6);
    assert(std::abs(cdf(n,-2.0) - 0.022750132) < 1e-6);
    assert(std::abs(cdf(n,3.0) - 0.998650102) < 1e-6);

    assert(std::abs(invCdf(n, 0.5) - 0.0) < 1e-6); 
    assert(std::abs(invCdf(n, .841344746) - 1.0) < 1e-6);
    assert(std::abs(invCdf(n, .158655254) + 1.0) < 1e-6);
    assert(std::abs(invCdf(n, 0.977249868) - 2.0) < 1e-6);
    assert(std::abs(invCdf(n, 0.022750132) + 2.0) < 1e-6);
    assert(std::abs(invCdf(n, 0.998650102) - 3.0) < 1e-6);

    assert(std::abs(cf(n,0.0).real() - 1.0) < 1e-12);
    assert(std::abs(cf(n,0.0).imag()) < 1e-12);

    assert(std::abs(cf(n,0.1).real() - .995012479) < 1e-9);
    assert(std::abs(cf(n,0.1).imag()) < 1e-12);

    assert(std::abs(cf(n,1.0).real() - .6065306597) < 1e-9);
    assert(std::abs(cf(n,1.0).imag()) < 1e-12);

    assert(std::abs(pdf(n,0.0) - numcpp::constants::ONE_OVER_SQRT_TWO_PI) < 1e-6); // pdf(0) for standard normal

    n = {2.0,3.0};

    assert(std::abs(cdf(n,2.0) - 0.5) < 1e-6);
    assert(std::abs(cdf(n,5.0) - 0.841344746) < 1e-6);
    assert(std::abs(cdf(n,-1.0) - 0.158655254) < 1e-6);

    assert(std::abs(invCdf(n,.5) - 2.0) < 1e-6);
    assert(std::abs(invCdf(n,.841344746) - 5.0) < 1e-6);
    assert(std::abs(invCdf(n,.158655254) + 1.0) < 1e-6);
   
    assert(std::abs(cf(n,0.0).real() - 1.0) < 1e-12);
    assert(std::abs(cf(n,0.0).imag()) < 1e-12);

    assert(std::abs(cf(n,0.1).real() - std::exp(-0.5 * 9.0 * 0.1 * 0.1) * std::cos(2.0 * 0.1)) < 1e-9);
    assert(std::abs(cf(n,0.1).imag() - std::exp(-0.5 * 9.0 * 0.1 * 0.1) * std::sin(2.0 * 0.1)) < 1e-9);

    assert(std::abs(pdf(n,2.0) - 1.0 / (3.0 * std::sqrt(2 * M_PI))) < 1e-6);
}

void testGammaDistribution() {

    using namespace numcpp::probability;
    Gamma g{2.0,2.0};
    
    assert(isClose(pdf(g,2.0),0.1839397206, 1e-9));
    assert(isClose(pdf(g,7.0),0.05284542099, 1e-9));
    assert(isClose(pdf(g,14.0),0.003191586879, 1e-9));

    assert(isClose(cdf(g,2.0),0.2642411177, 1e-9));
    assert(isClose(cdf(g,7.0),0.8641117746, 1e-9));
    assert(isClose(cdf(g,14.0),0.9927049443, 1e-9));


    
}

int main() {

    testNormalDistributions();
    testGammaDistribution();
    std::cout << "All tests for probability module have been passed!" << std::endl;
    return 0;
}