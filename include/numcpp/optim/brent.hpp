#pragma once
#include <functional>

namespace numcpp::optim {

    struct BrentResult {double x; double f; int iter;};

    inline BrentResult brent(const std::function<double(double)>& f, double ax, double bx, int maxIter = 100, double toleranceThreshold = 1e-12) {

        const double phi = 0.5 * (3.0 - std::sqrt(5.0));

        double a = ax;
        double b = bx;

        double x = a + phi * (b - a);
        double w = x;
        double v = x;

        double fx = f(x);
        double fw = fx;
        double fv = fx;

        double d = 0.0;  
        double e = 0.0;  
        int i =0; 
        for (int iter = 0; iter < maxIter; ++iter) {

            double m = 0.5 * (a + b);
            double tol1 = toleranceThreshold * std::abs(x) + 1e-12;
            double tol2 = 2.0 * tol1;

            if (std::abs(x - m) <= tol2 - 0.5 * (b - a)) break;

            bool parabolicStepAccepted = false;

            if (std::abs(e) > tol1) {
                
                double r = (x - w) * (fx - fv);
                double q = (x - v) * (fx - fw);
                double p = (x - v) * q - (x - w) * r;
                q = 2.0 * (q - r);

                if (q > 0.0) p = -p;
                q = std::abs(q);

                if (std::abs(p) < std::abs(0.5 * q * e) &&
                    p > q * (a - x) &&
                    p < q * (b - x)) {

                    d = p / q;
                    parabolicStepAccepted = true;
                }
            }

            if (!parabolicStepAccepted) {
                e = (x < m) ? (b - x) : (a - x);
                d = phi * e;
            }

            double u = x + ((std::abs(d) >= tol1) ? d : (d > 0 ? tol1 : -tol1));
            double fu = f(u);

            if (fu <= fx) {
                if (u < x) b = x;
                else       a = x;

                v = w; fv = fw;
                w = x; fw = fx;
                x = u; fx = fu;
            } else {
                if (u < x) a = u;
                else       b = u;

                if (fu <= fw || w == x) {
                    v = w; fv = fw;
                    w = u; fw = fu;
                } else if (fu <= fv || v == x || v == w) {
                    v = u; fv = fu;
                }
            }

            e = d;
            i++;
        }

        return {x,fx,i}; // return minimum value
    }

}