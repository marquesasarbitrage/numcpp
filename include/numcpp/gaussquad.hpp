#pragma once
#include "numcpp/objects.hpp"

namespace numcpp {

    namespace gaussquad {

        struct GaussLaguerreQuad {objects::Vector roots; objects::Vector weights;};

        inline GaussLaguerreQuad getGaussLaguerreQuadrature(int n) {

            int points = std::max(n,5);
            Eigen::MatrixXd J = Eigen::MatrixXd::Zero(points, points);

            for (int k = 0; k < points; k++)
            {
                J(k, k) = 2.0 * k + 1.0;
            }

            for (int k = 0; k < points - 1; k++)
            {
                double b = k + 1.0;
                J(k, k + 1) = b;
                J(k + 1, k) = b;
            }

            Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> solver(J);
            Eigen::VectorXd roots = solver.eigenvalues();
            Eigen::VectorXd weights(points);
            Eigen::MatrixXd V = solver.eigenvectors();
            for (int i = 0; i < points; i++)
            {
                weights[i] = std::pow(V(0, i), 2);
            }
            return {roots,weights};
        }

        inline double gaussLaguerreQuadIntegrate(const std::function<double(double)>& f,const GaussLaguerreQuad& quad) {

            double integral = 0.0;
            for (size_t i = 0; i < quad.roots.size(); i++){integral += quad.weights[i] * f(quad.roots[i]);}
            return integral;
        }
    }
}