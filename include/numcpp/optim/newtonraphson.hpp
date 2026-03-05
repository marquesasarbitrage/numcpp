#pragma once
#include <functional>
#include "numcpp/constants.hpp"

namespace numcpp {

    namespace optim {

        struct NewtonRaphsonResult {double x; double f; double fd; int iter;};

        inline NewtonRaphsonResult newtonRaphson(double x0, const std::function<double(double)>& f, const std::function<double(double)>& fDeriv, int maxIter = 100, double toleranceThreshold = 1e-12, double zeroThreshold = 1e-20) {

            double x = x0;
            double fX; 
            double dfX; 
            double xStep = 0;
            int numberIterations = 0;
            for (int i = 1; i <= maxIter; ++i) {

                numberIterations += 1;
                fX = f(x);
                dfX = fDeriv(x);
                if (std::abs(fX) < toleranceThreshold){
                    break;
                }
                
                if (std::abs(dfX) < zeroThreshold) {x=constants::DOUBLE_NAN,fX=constants::DOUBLE_NAN,dfX=constants::DOUBLE_NAN;break;}
                
                double xNew = x - fX / dfX;
                xStep = std::abs(xNew - x);
                x = xNew;

                if (xStep < toleranceThreshold) {break;}
                if (i==maxIter) break;
                
            }
            return {x,fX,dfX, numberIterations};
        }
    }
}