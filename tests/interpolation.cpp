#include <cassert>
#include <iostream>
#include "numcpp/interpolation.hpp"
#include "numcpp/errors.hpp"

void testInterpolation()
{
    // Sample data for interpolation
    std::map<double, double> testData = {
        {1.0, 2.0},
        {2.0, 3.0},
        {3.0, 5.0},
        {4.0, 7.0},
        {5.0, 11.0}
    };
    
    // Testing Linear Interpolation
    numcpp::interpolation::LinearInterpolation linear(testData);
    
    assert(linear.evaluate(1.0) == 2.0);
    assert(linear.evaluate(3.0) == 5.0);
    assert(linear.evaluate(5.0) == 11.0);

    
    // Testing Cubic Spline Interpolation
    numcpp::interpolation::CubicSpline cubic(testData);
    
    assert(cubic.evaluate(1.0) == 2.0);
    assert(cubic.evaluate(3.0) == 5.0);
    assert(cubic.evaluate(5.0) == 11.0);

    //testing cubic spline natural boundaries 
    assert(cubic.evaluateSecondDerivative(1.0) == 0.0);
    assert(cubic.evaluateSecondDerivative(5.0) == 0.0);

    // Testing out-of-range handling
    try {
        linear.evaluate(0.5); // Should throw EvaluationOutOfRange
        assert(false); // Should not reach this line
    } catch (const numcpp::errors::Error& e) {
        std::cout << e.what() << std::endl;
    }
    
    try {
        cubic.evaluate(6.0); // Should throw EvaluationOutOfRange
        assert(false);
    } catch (const numcpp::errors::Error& e) {
        std::cout << e.what() << std::endl;
    }




}

int main()
{
    testInterpolation();
    std::cout << "All tests for the interpolation module have been passed successfully!" << std::endl;
    return 0;
}