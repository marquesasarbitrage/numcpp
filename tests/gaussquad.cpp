#include <iostream>
#include <chrono>
#include <numcpp/gaussquad.hpp>
#include <cassert>

void testGaussLaguerre() {
    // Test the computation of roots and weights
    int points = 20;
    numcpp::gaussquad::GaussLaguerreQuad glq = numcpp::gaussquad::getGaussLaguerreQuadrature(points);

    // Verify the number of roots and weights
    std::cout << glq.roots.size() << std::endl;
    assert(glq.roots.size() == points);
    assert(glq.weights.size() == points);

    // Test the integration method
    auto f = [](double x) { return x; }; // Integral of x * exp(-x) from 0 to infinity is 1
    double result = numcpp::gaussquad::gaussLaguerreQuadIntegrate(f,glq);

    // Verify the integration result
    assert(std::abs(result - 1.0) < 1e-6);
}

void testGaussLaguerreHighDimension() {
    // Test the computation of roots and weights
    int points = 300;
    numcpp::gaussquad::GaussLaguerreQuad glq = numcpp::gaussquad::getGaussLaguerreQuadrature(points);

    // Verify the number of roots and weights
    assert(glq.roots.size() == points);
    assert(glq.weights.size() == points);

    // Test the integration method
    auto f = [](double x) { return x; }; // Integral of x * exp(-x) from 0 to infinity is 1
    double result = numcpp::gaussquad::gaussLaguerreQuadIntegrate(f,glq);

    // Verify the integration result
    assert(std::abs(result - 1.0) < 1e-6);
}

void testGaussLaguerreTime()
{
    // Test the timing of Gauss-Laguerre quadrature
    int points = 300;
    auto start = std::chrono::high_resolution_clock::now();
    numcpp::gaussquad::GaussLaguerreQuad glq = numcpp::gaussquad::getGaussLaguerreQuadrature(points);
    auto end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed = end - start;
    std::cout << "Time taken to compute Gauss-Laguerre quadrature with " << points << " points: " << elapsed.count() << " seconds" << std::endl;
}

int main() {
    testGaussLaguerre();
    testGaussLaguerreHighDimension();
    testGaussLaguerreTime();
    std::cout << "All tests for the gauss quad module have been passed successfully!" << std::endl;
    return 0;
}