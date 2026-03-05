#pragma once 
#include "numcpp/objects/polynomial.hpp"
#include "numcpp/errors.hpp"

namespace numcpp {

    namespace polysolver {
        
        inline double getQuadraticExtremum(const objects::Polynomial& quadraticPolynomial) {

            if (quadraticPolynomial.getDegree() != 2) errors::throwInvalidInputr("Polynomial is not quadratic");
            double a = quadraticPolynomial.getCoefficient(2), b = quadraticPolynomial.getCoefficient(1);
            return -b/(2.0*a);
        }

        inline std::vector<double> getQuadraticRoots(const objects::Polynomial& quadraticPolynomial) {

            if (quadraticPolynomial.getDegree() != 2) errors::throwInvalidInputr("Polynomial is not quadratic");
            double a = quadraticPolynomial.getCoefficient(2), b = quadraticPolynomial.getCoefficient(1), c = quadraticPolynomial.getCoefficient(0);
            double delta = b * b - 4 * a * c;
            if (delta < 0.0) {

                return {};
                
            } else if (std::abs(delta) < 1e-15) {

                return {std::vector<double>({-b/(2.0*a)})};

            } else {
                
                double q = -0.5 * (b + std::copysign(std::sqrt(delta), b));
                return {std::vector<double>({q / a, c / q})};

            }
        }

        inline std::vector<double> getCompanionPolynomialRoots(const objects::Polynomial& polynomial) {

            Eigen::EigenSolver<Eigen::MatrixXd> solver(polynomial.getCompanionMatrix());
            Eigen::VectorXcd eigenvalues = solver.eigenvalues();
            std::vector<double> realRoots;
            realRoots.reserve(eigenvalues.size());

            constexpr double EPS = 1e-12;

            for (int i = 0; i < eigenvalues.size(); ++i) {
                if (std::abs(eigenvalues[i].imag()) < EPS) {
                    realRoots.push_back(eigenvalues[i].real());
                }
            }

            return {realRoots};
        }

        inline std::vector<double> getCompanionPolynomialExtremums(const objects::Polynomial& polynomial) {

            objects::Polynomial firstDerivative = polynomial.getFirstDerivative();
            std::vector<double> firstDerivativeRoots = getCompanionPolynomialRoots(firstDerivative);
            std::vector<double> output;
            objects::Polynomial secondDerivative = firstDerivative.getFirstDerivative();
            double value;
            for (const double x: firstDerivativeRoots) {
                value = secondDerivative.getValue(x);
                
                output.push_back(x);
            }
            return output; 
        }

    }
}
