#include <cassert>
#include <cmath>
#include <iostream>
#include "numcpp/polysolver.hpp"

using namespace numcpp::objects;

bool isClose(double a, double b, double tol = 1e-9) {
    return std::abs(a - b) <= tol;
}

void testQuadraticPolynomialSolver() {

    double extremum1 = numcpp::polysolver::getQuadraticExtremum(Polynomial({{2,1.0}, {1,-4.0}, {0,5.0}}));
    assert(isClose(extremum1,2.0));

    double extremum2 = numcpp::polysolver::getQuadraticExtremum(Polynomial({{2,-1.0}, {1,4.0}, {0,-1.0}}));
    assert(isClose(extremum2,2.0));

    std::vector<double> roots1 = numcpp::polysolver::getQuadraticRoots(Polynomial({{2,2.0}, {1,20.0}, {0,25.0}}));
    assert(roots1.size() ==2);
    assert(isClose(roots1[0], -8.53553, 1e-4));
    assert(isClose(roots1[1], -1.46447, 1e-4));

    std::vector<double> roots2 = numcpp::polysolver::getQuadraticRoots(Polynomial({{2,1.0}, {1,1.0}, {0,1.0}}));
    assert(roots2.size() ==0);

    std::vector<double> roots3 = numcpp::polysolver::getQuadraticRoots(Polynomial({{2,1.0}, {1,2.0}, {0,1.0}}));
    double extremum3 = numcpp::polysolver::getQuadraticExtremum(Polynomial({{2,1.0}, {1,2.0}, {0,1.0}}));
    assert(roots3.size() ==1);
    assert(isClose(roots3[0], -1.0));
    assert(isClose(roots3[0], extremum3));
}

void testQuadraticCompanionPolynomialSolver() {

    Polynomial p({{2,-1.0}, {1,4.0}, {0,-1.0}});
    std::vector<double> roots = numcpp::polysolver::getCompanionPolynomialRoots(p); 
    assert(roots.size()==2); 
    assert(isClose(roots[0], 2-std::sqrt(3), 1e-6));
    assert(isClose(roots[1], 2+std::sqrt(3), 1e-6));
    std::vector<double> extremum = numcpp::polysolver::getCompanionPolynomialExtremums(p); 
    assert(extremum.size()==1);
    assert(isClose(extremum[0], 2.0, 1e-6));

}

void testQuarticCompanionPolynomialSolver() {

    Polynomial p({{4,1.0}, {3,-10.0}, {2,35.0}, {1,-50.0}, {0,24.0}});
    std::vector<double> roots = numcpp::polysolver::getCompanionPolynomialRoots(p); 
    assert(roots.size()==4); 
    assert(isClose(roots[0], 1.0, 1e-6));
    assert(isClose(roots[1], 2.0, 1e-6));
    assert(isClose(roots[2], 3.0, 1e-6));
    assert(isClose(roots[3], 4.0, 1e-6));
    std::vector<double> extremum = numcpp::polysolver::getCompanionPolynomialExtremums(p); 
    assert(extremum.size()==3);
    assert(isClose(extremum[2] , 2.5 + std::sqrt(5)/2, 1e-6));
    assert(isClose(extremum[0], 2.5 - std::sqrt(5)/2, 1e-6));
    assert(isClose(extremum[1] , 2.5, 1e-6));
   
}

int main() {

    testQuadraticPolynomialSolver();
    testQuadraticCompanionPolynomialSolver();
    testQuarticCompanionPolynomialSolver();
    std::cout << "All tests for the Polynomial object and PolynomialSolver functions have been passed successfully!" << std::endl;
    return 0;
}