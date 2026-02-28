#pragma once 
#include "numcpp/constants.hpp"
#include "numcpp/functions.hpp"
#include "numcpp/objects.hpp"
#include "numcpp/optim.hpp"
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
        inline double pdf(Uniform params, double x) {return (x>=params.a and x<=params.b) ? 1.0/(params.b-params.a) : 0.0;}
        inline double cdf(Uniform params, double x) {return (x<params.a) ? 0.0 : (x>params.b) ? 1.0 : (x-params.a)/(params.b-params.a);}
        inline double invCdf(Uniform params, double p) {return params.a + p*(params.b-params.a);}
        inline theoricalMoments getTheoricalMoments(Uniform params) {return {.5*(params.a+params.b),(params.b-params.a)*(params.b-params.a)/12.0,0.0,-6.0/5.0};}
        inline objects::Vector sample(Uniform params, int n, std::mt19937& gen) {
            objects::Vector vec(n);
            std::uniform_real_distribution<> dist(params.a, params.b);
            for (int i = 0; i < n; ++i) {vec[i] = dist(gen);}
            return vec;
        }
        inline double sample(Uniform params, std::mt19937& gen) {std::uniform_real_distribution<> dist(params.a, params.b);return dist(gen);}

        struct Normal {double mean=0.0, sigma=1.0;};
        template<>
        struct probabilityDistributionCategory<Normal> {using type = continuousProbabilityTag;};

        inline double pdf(Normal params, double x) {double z =(x -params.mean)/params.sigma;return constants::ONE_OVER_SQRT_TWO_PI*std::exp(-.5*z*z)/params.sigma;}
        inline double cdf(Normal params, double x) {return .5 * (1.0 + std::erf(constants::ONE_OVER_SQRT_TWO*(x-params.mean)/params.sigma));}
        inline double invCdf(Normal params, double p) {return params.mean+params.sigma*functions::acklamStandardGaussianInverseCdf(p);}
        inline std::complex<double> cf(Normal params, double t) {std::complex<double> i(0.0, 1.0); return std::exp(i * t * params.mean - 0.5 * params.sigma * params.sigma * t * t);}
        inline theoricalMoments getTheoricalMoments(Normal params) {return {params.mean,params.sigma*params.sigma,0.0,0.0};}
        inline objects::Vector sample(Normal params, int n, std::mt19937& gen) {
            objects::Vector u = sample(Uniform{},n,gen); 
            return u.unaryExpr([params](double x) { return invCdf(params,std::max(x, std::numeric_limits<double>::min())); });
        }
        inline double sample(Normal params, std::mt19937& gen) {return invCdf(params,sample(Uniform{},gen));}
     

        struct Gamma {double k, theta;}; 
        template<>
        struct probabilityDistributionCategory<Gamma> {using type = continuousProbabilityTag;};

        inline theoricalMoments getTheoricalMoments(Gamma params) {return {params.k*params.theta,params.k*params.theta*params.theta,2/std::sqrt(params.k),6/params.k};}
        inline double pdf(Gamma params, double x) {return std::pow(x,params.k-1.0)*std::exp(-x/params.theta)/(std::tgamma(params.k)*std::pow(params.theta,params.k));}
        inline double cdf(Gamma params, double x) {return functions::incompleteGamma(params.k, x/params.theta)/std::tgamma(params.k);}
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
        inline objects::Vector sample(Gamma params, int n, std::mt19937& gen) {

            objects::Vector vec(n);
            for (int i = 0; i < n; ++i) {vec[i] = sample(params,gen);}
            return vec;
        }

        struct ChiSquared{double nu; Gamma gamma() const {return Gamma{nu/2.0,2.0};}};
        template<>
        struct probabilityDistributionCategory<ChiSquared> {using type = continuousProbabilityTag;};
        inline double pdf(ChiSquared params, double x) {return pdf(params.gamma(),x);}
        inline double cdf(ChiSquared params, double x) {return cdf(params.gamma(),x);}
        inline double sample(ChiSquared params, std::mt19937& gen) { return sample(params.gamma(),gen); }
        inline theoricalMoments getTheoricalMoments(ChiSquared params) {return getTheoricalMoments(params.gamma());}
        inline objects::Vector sample(ChiSquared params, int n, std::mt19937& gen) { return sample(params.gamma(),n,gen); };
        
        struct Student{double nu;};
        template<>
        struct probabilityDistributionCategory<Student> {using type = continuousProbabilityTag;};
        inline double pdf(Student params, double x){ return std::pow(1.0+x*x/params.nu,-.5*(params.nu+1.0))*std::tgamma(.5*(params.nu+1.0))/(std::sqrt(params.nu*constants::PI)*std::tgamma(.5*params.nu));}
        inline double cdf(Student params, double x) {
            return x>=0.0 ? 1.0-.5*functions::incompleteBeta(params.nu/2.0,.5,params.nu/(params.nu+x*x)) : .5*functions::incompleteBeta(params.nu/2.0,.5,params.nu/(params.nu+x*x));
        }
        
        inline theoricalMoments getTheoricalMoments(Student params){
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

        inline objects::Vector sample(Student params, int n, std::mt19937& gen) {
            objects::Vector vec(n);
            for (int i = 0; i < n; ++i) {vec[i] = sample(params,gen);}
            return vec;
        }

        inline double invCdf(Student params,double p) {

            double z0 = functions::cornishFisherStudentInverseCdfExpansion(p, params.nu);
            std::function<double(double)> target = [params,p](double x) { return p-cdf(params,x); };
            std::function<double(double)> targetDeriv = [params,p](double x) { return -pdf(params,x); };
            optim::NewtonRaphsonResult result = optim::newtonRaphson(z0, target, targetDeriv,20,1e-12);
            return result.x;
        }


    }
}