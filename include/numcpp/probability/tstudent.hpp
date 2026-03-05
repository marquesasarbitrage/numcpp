#pragma once 
#include "numcpp/objects/base.hpp"
#include <random>
#include <numcpp/probability/chisquared.hpp>
#include <numcpp/optim/newtonraphson.hpp>

namespace numcpp {

    namespace probability {

        struct Student{double nu;};

        inline double logbeta(double a,double b) {return std::lgamma(a + b)- std::lgamma(a)- std::lgamma(b);}

        inline double incompleteBetaContinuedFraction(double a, double b, double x) {

            double qab = a + b;
            double qap = a + 1.0;
            double qam = a - 1.0;

            double c = 1.0;
            double d = 1.0 - qab * x / qap;
            if (std::fabs(d) < std::numeric_limits<double>::min()) d = std::numeric_limits<double>::min();
            d = 1.0 / d;
            double h = d;

            for (int m = 1; m <= 200; ++m)
            {
                int m2 = 2 * m;

                // even step
                double aa = m * (b - m) * x /
                            ((qam + m2) * (a + m2));

                d = 1.0 + aa * d;
                if (std::fabs(d) < std::numeric_limits<double>::min()) d = std::numeric_limits<double>::min();
                c = 1.0 + aa / c;
                if (std::fabs(c) < std::numeric_limits<double>::min()) c = std::numeric_limits<double>::min();
                d = 1.0 / d;
                h *= d * c;

                // odd step
                aa = -(a + m) * (qab + m) * x /
                    ((a + m2) * (qap + m2));

                d = 1.0 + aa * d;
                if (std::fabs(d) < std::numeric_limits<double>::min()) d = std::numeric_limits<double>::min();
                c = 1.0 + aa / c;
                if (std::fabs(c) < std::numeric_limits<double>::min()) c = std::numeric_limits<double>::min();
                d = 1.0 / d;
                double del = d * c;
                h *= del;

                if (std::fabs(del - 1.0) < 1e-14)
                    break;
            }

            return h;
        }

        inline double incompleteBeta(double a, double b, double x) {
            
            if (x == 0.0) return 0.0;
            if (x == 1.0) return 1.0;

            bool flip = false;
            if (x > (a + 1.0) / (a + b + 2.0))
            {
                flip = true;
                std::swap(a, b);
                x = 1.0 - x;
            }

            double front = std::exp(
                a * std::log(x) +
                b * std::log(1.0 - x) +
                logbeta(a,b)
            ) / a;

            double cf = incompleteBetaContinuedFraction(a, b, x);
            double result = front * cf;

            return flip ? 1.0 - result : result;
        }
        
        inline double pdf(Student params, double x){ return std::pow(1.0+x*x/params.nu,-.5*(params.nu+1.0))*std::tgamma(.5*(params.nu+1.0))/(std::sqrt(params.nu*constants::PI)*std::tgamma(.5*params.nu));}
        
        inline double cdf(Student params, double x) {
            return x>=0.0 ? 1.0-.5*incompleteBeta(params.nu/2.0,.5,params.nu/(params.nu+x*x)) : .5*incompleteBeta(params.nu/2.0,.5,params.nu/(params.nu+x*x));
        }

        inline double cornishFisherStudentInverseCdfExpansion(double p, double nu) {

            double z = acklamStandardGaussianInverseCdf(p);

            double z2 = z*z;
            double z3 = z2*z;
            double z5 = z3*z2;

            return z + (z3 + z)/(4.0*nu)
                        + (5.0*z5 + 16.0*z3 + 3.0*z)/(96.0*nu*nu);
        }
        
        inline objects::TheoricalMoments theoricalMoments(Student params){
            return {
                0.0,
                params.nu>2 ? params.nu/(params.nu-2) : constants::DOUBLE_NAN,
                params.nu>3 ? 0.0: constants::DOUBLE_NAN,
                params.nu>4 ? 6.0/(params.nu-4.0) : constants::DOUBLE_NAN
            };
        }

        inline double sample(Student params, std::mt19937& gen) {

            return sample(Normal{},gen)/std::sqrt(sample(ChiSquared{params.nu},gen)/params.nu);
        }

        inline objects::Vector sample(Student params, std::mt19937& gen, int n) {
            objects::Vector vec(n);
            for (int i = 0; i < n; ++i) {vec[i] = sample(params,gen);}
            return vec;
        }

        inline double invCdf(Student params,double p) {

            double z0 = cornishFisherStudentInverseCdfExpansion(p, params.nu);
            std::function<double(double)> target = [params,p](double x) { return p-cdf(params,x); };
            std::function<double(double)> targetDeriv = [params,p](double x) { return -pdf(params,x); };
            optim::NewtonRaphsonResult result = optim::newtonRaphson(z0, target, targetDeriv,20,1e-12);
            return result.x;
        }

    }
}