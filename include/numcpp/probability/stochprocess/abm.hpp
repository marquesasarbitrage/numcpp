#pragma once 
#include <map>
#include <numcpp/stats.hpp>
#include <numcpp/objects/base.hpp>
#include <numcpp/probability/normal.hpp>

namespace numcpp {

    namespace probability {

        namespace stochprocess {

            struct ArithmeticBrownianMotion {double x0, mu, sigma;}; 

            inline double mean(const ArithmeticBrownianMotion& abm, double T){return abm.x0 + abm.mu*T;}

            inline double variance(const ArithmeticBrownianMotion& abm, double T){return abm.sigma*abm.sigma*T;}
            

            inline double sampleStep(const ArithmeticBrownianMotion& abm, double x, double dt, double z) {

                return x + abm.mu * dt + abm.sigma * std::sqrt(dt) * z;
            }

            inline double sampleStep(const ArithmeticBrownianMotion& abm, double x, double dt,std::mt19937& gen) {

                return sampleStep(abm,x,dt,sample(Normal{},gen));
            }

            inline numcpp::objects::Vector samplePath(const ArithmeticBrownianMotion& abm, double T, const numcpp::objects::Vector& Z) {

                double dt = T/double(Z.size());
                numcpp::objects::Vector x(Z.size()); 
                x(0) = abm.x0; 
                for (size_t i = 1; i<Z.size();i++) { x(i) = x(i-1) + abm.mu * dt + abm.sigma * std::sqrt(dt) * Z[i];}
                return x; 
            }

            inline numcpp::objects::Vector samplePath(const ArithmeticBrownianMotion& abm, double T, int numberSteps, std::mt19937& gen) {

                return samplePath(abm,T,sample(Normal{},gen,numberSteps));
            }

            inline numcpp::objects::Matrix samplePath(const ArithmeticBrownianMotion& abm, double T, const numcpp::objects::Matrix& Z) {

                size_t steps = Z.cols(),paths = Z.rows();
                double dt = T / double(steps);
                numcpp::objects::Matrix X(paths, steps);
                for (size_t j = 0; j < paths; ++j) X(j, 0) = abm.x0;
                for (size_t i = 1; i < steps; ++i) {
                    for (size_t j = 0; j < paths; ++j) {
                        X(j,i) = X(j,i-1) + abm.mu * dt + abm.sigma * std::sqrt(dt) * Z(j,i);
                    }
                }
                return X;
            }

            inline numcpp::objects::Matrix samplePath(const ArithmeticBrownianMotion& abm, double T, int numberSteps, int size, std::mt19937& gen) {

                return samplePath(abm,T,sample(Normal{},gen,size,numberSteps));
            }

            inline ArithmeticBrownianMotion fitArithmeticBrownianMotion(const std::map<double,double>& data) {

                numcpp::objects::Vector X(data.size()-1); 
                double mu = (data.rbegin()->second - data.begin()->second)/(data.rbegin()->first - data.begin()->first); 
                double x0 = data.rbegin()->second; 
                numcpp::objects::Vector residuals(data.size()-1); 
                size_t i = 0;
                for (auto it = std::next(data.begin()); it != data.end(); ++it) {
                    auto prev = std::prev(it);
                    double sqdt = std::sqrt(it->first - prev->first);
                    residuals(i) = (it->second - prev->second)/sqdt-sqdt*mu;
                    ++i;
                }
                return {x0,mu,numcpp::stats::populationStandardDeviation(residuals)};
            }
        }
    }
}
