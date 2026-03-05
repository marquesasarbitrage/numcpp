#pragma once
#include <cmath>
#include <random>
#include "numcpp/probability/normal.hpp"

namespace numcpp {

    namespace probability {

        struct Gamma {double k, theta;}; 

        inline double incompleteGammaSeries(double a, double x, int maxIter, double eps) {
            double sum = 1.0 / a;
            double term = sum;
            for (int n = 1; n < maxIter; ++n) {
                term *= x / (a + n);
                sum += term;
                if (std::fabs(term) < eps * std::fabs(sum)) break;
            }
            return sum * std::exp(-x + a * std::log(x) - std::lgamma(a));
        }

        inline double incompleteGammaContinuedFraction(double a, double x, int maxIter, double eps) {
            double b = x + 1.0 - a;
            double c = 1.0 / std::numeric_limits<double>::min();  
            double d = 1.0 / b;
            double h = d;
            for (int i = 1; i < maxIter; ++i) {
                double an = -i * (i - a);
                b += 2.0;
                d = an * d + b;
                if (std::fabs(d) < std::numeric_limits<double>::min()) d = std::numeric_limits<double>::min();
                c = b + an / c;
                if (std::fabs(c) < std::numeric_limits<double>::min()) c = std::numeric_limits<double>::min();
                d = 1.0 / d;
                double delta = d * c;
                h *= delta;
                if (std::fabs(delta - 1.0) < eps) break;
            }
            return 1.0 - std::exp(-x + a * std::log(x) - std::lgamma(a)) * h;
        }

        inline double incompleteGamma(double a, double x) {
            if (x < 0.0 || a <= 0.0) return 0.0;
            if (x == 0.0) return 0.0;
            if (x < a + 1.0) return incompleteGammaSeries(a, x,10000,1e-14);   
            else return incompleteGammaContinuedFraction(a, x,10000,1e-14);           
        }
    
        inline objects::TheoricalMoments theoricalMoments(Gamma params) {return {params.k*params.theta,params.k*params.theta*params.theta,2/std::sqrt(params.k),6/params.k};}
        
        inline double pdf(Gamma params, double x) {return std::pow(x,params.k-1.0)*std::exp(-x/params.theta)/(std::tgamma(params.k)*std::pow(params.theta,params.k));}

        inline double cdf(Gamma params, double x) {return incompleteGamma(params.k, x/params.theta)/std::tgamma(params.k);}

        inline double sample(Gamma params, std::mt19937& gen) {

            if (params.k < 1.0) {
                double u = sample(Uniform{},gen);
                return sample(Gamma{params.k + 1.0, params.theta},gen) * std::pow(u, 1.0/params.k);
            }

            double d = params.k - 1.0/3.0;
            double c = 1.0 / std::sqrt(9.0*d);

            while (true) {
                double z = sample(Normal{},gen);
                double v = 1.0 + c*z;

                if (v <= 0.0) continue;

                v = v*v*v;

                double u = sample(Uniform{},gen);
                u = std::max(u, std::numeric_limits<double>::min());

                if (u < 1.0 - 0.0331*(z*z)*(z*z))
                    return params.theta*d*v;

                if (std::log(u) < 0.5*z*z + d*(1 - v + std::log(v)))
                    return params.theta*d*v;
            }
        };

        inline objects::Vector sample(Gamma params, std::mt19937& gen, int n) {

            objects::Vector vec(n);
            for (int i = 0; i < n; ++i) {vec[i] = sample(params,gen);}
            return vec;
        }
        
    }
}