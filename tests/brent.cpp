#include "numcpp/probability/uniform.hpp"
#include <numcpp/optim/brent.hpp>
#include <iostream>
#include <cassert>


numcpp::optim::BrentResult bruteMinimum(double a, double b, const std::function<double(double)>& f, int N) {

    double minVal=1e100, minX;

    for (int i = 0; i <= N; ++i) {
        double x = a + (b - a) * i / N;
        if (f(x)<minVal) {minVal = f(x); minX = x;}
    }
    return {minX, f(minX), 100};
}

std::function<double(double)> randomQuadratic(std::mt19937 & gen) {

    double a = numcpp::probability::sample(numcpp::probability::Uniform{-5,5},gen);

    auto f = [&](double x) {
        return (x - a)*(x - a);
    };

    return f;
}

void testBrentRandomQuadratic() {

    std::random_device rd;      
    std::mt19937 gen(rd());
    for (int i = 0; i<1000; i++) {
       
        std::function<double(double)> f = randomQuadratic(gen); 
        numcpp::optim::BrentResult res = numcpp::optim::brent(f, -10.0, 10.0,50,1e-10); 
        assert(std::abs(res.f) < 1e-6);
    }
}

void testOscilllaryFunction() {

    std::function<double(double)> f = [](double x) {
        return std::sin(x) + 0.1*x;
    };

    numcpp::optim::BrentResult brut = bruteMinimum(-5.0,5.0,f,10000);
    numcpp::optim::BrentResult brent = numcpp::optim::brent(f,-5.0,5.0,50,1e-15);

    assert(std::abs(brut.x-brent.x)<1e-4);
    assert(std::abs(brut.f-brent.f)<1e-5);
};

void testScaledQuadratic() {

    auto f = [](double x) {
        return 3.0*(x + 1.5)*(x + 1.5) + 5.0;
    };

    numcpp::optim::BrentResult brut = bruteMinimum(-5.0,5.0,f,10000);
    numcpp::optim::BrentResult brent = numcpp::optim::brent(f,-5.0,5.0,50,1e-15);

    assert(std::abs(brut.x-brent.x)<1e-4);
    assert(std::abs(brut.f-brent.f)<1e-5);
};

int main() {

    testBrentRandomQuadratic();
    testOscilllaryFunction();
    testScaledQuadratic();
    std::cout << "All tests for the brent optimizer have been passed successfully!" << std::endl;

}