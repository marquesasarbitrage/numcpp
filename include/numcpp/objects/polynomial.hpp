#pragma once
#include <map> 
#include <numcpp/objects/base.hpp>
#include "numcpp/errors.hpp"

namespace numcpp {

    namespace objects {

        struct Polynomial {

            std::map<int,double> polynomialMap_;

            Polynomial(std::map<int,double> polynomialMap) {

                if (polynomialMap.size()==0) errors::throwInvalidInputr("Polynomial map cannot be empty");
                for (auto const& [key, val] : polynomialMap) {

                    if (key>=0) { polynomialMap_[key] = val; }
                    else errors::throwInvalidInputr("Polynomial power cannot be negative");
                }

                while (std::abs(polynomialMap_.rbegin()->second)<=1e-12) {

                    polynomialMap_.erase(polynomialMap_.rbegin()->first);
                    if (polynomialMap_.size()==0) errors::throwInvalidInputr("Polynomial map cannot be empty");
                }
            }

            Polynomial getFirstDerivative() const {

                std::map<int, double> mapDerivative; 
                for (auto const& [key, val] : polynomialMap_) { if ( (key-1)>=0 ) mapDerivative[key-1] = key*val; }
                return {mapDerivative};
            }
            double getValue(double x) const {

                double result = 0.0; 
                for (auto const& [key, val] : polynomialMap_) { result += val * std::pow(x,key);}
                return result;
            }

            int getDegree() const { return polynomialMap_.rbegin()->first; }

            double getLeadingCoefficient() const { return polynomialMap_.rbegin()->second; }

            double getCoefficient(int x) const {

                auto it = polynomialMap_.find(x);
                return (it != polynomialMap_.end()) ? it->second : 0.0;
            }

            Eigen::MatrixXd getCompanionMatrix() const {

                int degree = getDegree();
                Eigen::MatrixXd C = Eigen::MatrixXd::Zero(degree, degree);

                for (int i = 1; i < degree; ++i)
                    C(i, i - 1) = 1.0;

                for (int power = 0; power < degree; ++power) {
                    C(power, degree - 1) = -getCoefficient(power) / getLeadingCoefficient();
                }

                return C;
            }
        };
    }
}