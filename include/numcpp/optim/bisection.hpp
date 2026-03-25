#pragma once 
#include <cmath>
#include <functional>
#include <numcpp/constants.hpp>

namespace numcpp::optim {

    struct BisectionResult {
        double x;
        double f;
        int iter;
    };

    inline BisectionResult bisection(
        double a,
        double b,
        const std::function<double(double)>& f,
        int maxIter = 100,
        double toleranceThreshold = 1e-12) {

        double fa = f(a);
        double fb = f(b);

        if (fa * fb > 0) {
            return {constants::DOUBLE_NAN, constants::DOUBLE_NAN, 0};
        }

        double x = a;
        double fx = fa;
        int numberIterations = 0;

        for (int i = 1; i <= maxIter; ++i) {
            numberIterations++;

            x = 0.5 * (a + b);
            fx = f(x);

            if (std::abs(fx) < toleranceThreshold) break;
            if (std::abs(b - a) < toleranceThreshold) break;

            if (fa * fx < 0) {
                b = x;
                fb = fx;
            } else {
                a = x;
                fa = fx;
            }
        }

        return {x, fx, numberIterations};
    }
}