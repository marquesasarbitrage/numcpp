#pragma once 
#include "numcpp/objects/base.hpp"
#include <random>

namespace numcpp {

    namespace probability {

        struct Uniform {double a=0.0, b=1.0;};

        inline double pdf(Uniform params, double x) {return (x>=params.a and x<=params.b) ? 1.0/(params.b-params.a) : 0.0;}

        inline double cdf(Uniform params, double x) {return (x<params.a) ? 0.0 : (x>params.b) ? 1.0 : (x-params.a)/(params.b-params.a);}

        inline double invCdf(Uniform params, double p) {return params.a + p*(params.b-params.a);}

        inline objects::TheoricalMoments theoricalMoments(Uniform params) {return {.5*(params.a+params.b),(params.b-params.a)*(params.b-params.a)/12.0,0.0,-6.0/5.0};}
        
        inline objects::Vector sample(Uniform params, std::mt19937& gen, int n) {
            objects::Vector vec(n);
            std::uniform_real_distribution<> dist(params.a, params.b);
            for (int i = 0; i < n; ++i) {vec[i] = dist(gen);}
            return vec;
        }

        inline double sample(Uniform params, std::mt19937& gen) {
            std::uniform_real_distribution<> dist(params.a, params.b);
            return dist(gen);
        }
    
    }
}