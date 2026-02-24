#pragma once 
#include "numcpp/functions.hpp"
#include "numcpp/objects.hpp"
#include <cmath>
#include <complex>
#include <random>

namespace numcpp {

    namespace probability {

        template<typename T>
        struct probabilityDistributionCategory;

        struct theoricalMoments {double mean,variance,skewness,excessKurtosis;};

        struct continuousProbabilityTag {};
        struct discreteProbabilityTag {};
        struct multivariateProbabilityTag {};

        struct Uniform {double a=0.0, b=1.0;};
        template<>
        struct probabilityDistributionCategory<Uniform> {using type = continuousProbabilityTag;};
        inline double pdf(const Uniform& params, double x) {return (x>=params.a and x<=params.b) ? 1.0/(params.b-params.a) : 0.0;}
        inline double cdf(const Uniform& params, double x) {return (x<params.a) ? 0.0 : (x>params.b) ? 1.0 : (x-params.a)/(params.b-params.a);}
        inline double invCdf(const Uniform& params, double p) {return params.a + p*(params.b-params.a);}
        inline theoricalMoments getTheoricalMoments(const Uniform& params) {return {.5*(params.a+params.b),(params.b-params.a)*(params.b-params.a)/12.0,0.0,-6.0/5.0};}
        inline objects::Vector sample(const Uniform& params, int n, std::mt19937& gen) {
            objects::Vector vec(n);
            std::uniform_real_distribution<> dist(params.a, params.b);
            for (int i = 0; i < n; ++i) {vec[i] = dist(gen);}
            return vec;
        }

        // Gaussian Normal 
        struct Normal {double mean=0.0, sigma=1.0;};
        template<>
        struct probabilityDistributionCategory<Normal> {using type = continuousProbabilityTag;};

        inline double pdf(const Normal& params, double x) {double z =(x -params.mean)/params.sigma;return constants::ONE_OVER_SQRT_TWO_PI*std::exp(-.5*z*z)/params.sigma;}
        inline double cdf(const Normal& params, double x) {return .5 * (1.0 + std::erf(constants::ONE_OVER_SQRT_TWO*(x-params.mean)/params.sigma));}
        inline double invCdf(const Normal& params, double p) {return params.mean+params.sigma*functions::acklamStandardGaussianInverseCdf(p);}
        inline std::complex<double> cf(const Normal& params, double t) {std::complex<double> i(0.0, 1.0); return std::exp(i * t * params.mean - 0.5 * params.sigma * params.sigma * t * t);}
        inline theoricalMoments getTheoricalMoments(const Normal& params) {return {params.mean,params.sigma*params.sigma,0.0,0.0};}
        inline objects::Vector sample(const Normal& params, int n, std::mt19937& gen) {
            objects::Vector u = sample(Uniform{},n,gen); 
            return u.unaryExpr([params](double x) { return invCdf(params,x); });
        }
     

        struct Gamma {double k, theta;}; 
        template<>
        struct probabilityDistributionCategory<Gamma> {using type = continuousProbabilityTag;};

        inline theoricalMoments getTheoricalMoments(const Gamma& params) {return {params.k*params.theta,params.k*params.theta*params.theta,2/std::sqrt(params.k),6/params.k};}
        inline double pdf(const Gamma& params, double x) {return std::pow(x,params.k-1.0)*std::exp(-x/params.theta)/(std::tgamma(params.k)*std::pow(params.theta,params.k));}
        inline double cdf(const Gamma& params, double x) {return functions::incompleteGamma(params.k, x/params.theta)/std::tgamma(params.k);}
        inline objects::Vector sample(const Gamma& params,int n);

    }
}