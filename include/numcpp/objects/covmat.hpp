#pragma once 
#include <numcpp/objects/base.hpp>
#include "numcpp/errors.hpp"

namespace numcpp {

    namespace objects {

        struct CovarianceMatrix;

        struct CorrelationMatrix {

            Matrix matrix_; 
            CorrelationMatrix(): matrix_(Eigen::MatrixXd::Identity(1, 1)) {}

            CorrelationMatrix(double coeff) {

                if (abs(coeff) > 1.0) errors::throwInvalidInputr("Correlation coefficient must be between -1 and 1");
                matrix_ = Eigen::MatrixXd::Zero(2,2);
                matrix_ << 1.0, coeff, 
                            coeff, 1.0;
            }

            CorrelationMatrix(const Vector& coeffs) {

                size_t c = coeffs.size();
                int p = static_cast<int>(1.0 + std::sqrt(1.0 + 8.0 * coeffs.size())) / 2.0;
                if (p * (p - 1) / 2 != coeffs.size()) errors::throwInvalidInputr("The number of correlation coeficients is invalid");

                matrix_ = Eigen::MatrixXd::Identity(p, p);

                int idx = 0;
                for (int i = 0; i < p; ++i) {
                    for (int j = i + 1; j < p; ++j) {
                        if (abs(coeffs[idx]) > 1.0) errors::throwInvalidInputr("Correlation coefficient must be between -1 and 1");
                        matrix_(i, j) = coeffs[idx];
                        matrix_(j, i) = coeffs[idx]; 
                        ++idx;
                    }
                }

                if (matrix_.rows() != matrix_.cols()) errors::throwInvalidInputr("Correlation matrix must be square");
                if (!((matrix_ - matrix_.transpose()).norm() < 1e-10)) errors::throwInvalidInputr("Correlation matrix must be symmetric");

                Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> solver(matrix_);
                if (!((solver.eigenvalues().array() > 0).all())) errors::throwInvalidInputr("Correlation matrix must be positive semi definite");

            }

            CorrelationMatrix(const Matrix& matrix) {

                if (matrix.rows() != matrix.cols()) errors::throwInvalidInputr("Correlation matrix must be square");
                if (!((matrix - matrix.transpose()).norm() < 1e-10)) errors::throwInvalidInputr("Correlation matrix must be symmetric");
                Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> solver(matrix);
                if (!((solver.eigenvalues().array() > 0).all())) errors::throwInvalidInputr("Correlation matrix must be positive semi definite");
                Vector variances = matrix.diagonal(); 
                for (double v: variances) if (v!=1.0) errors::throwInvalidInputr("Correlation matrix must have its variances standardized to 1");
                for (size_t i = 0; i < matrix.rows(); ++i) {
                    for (size_t j = i + 1; j < matrix.rows(); ++j) {
                        if (std::abs(matrix(i, j)) > 1.0) errors::throwInvalidInputr("Correlation coefficient must be between -1 and 1");
                    }
                }
                matrix_ = matrix;

            }

            double getCorrelation(int i, int j) const { return matrix_(i,j); }
            Matrix getEigenVectors() const {Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> solver(matrix_);return solver.eigenvectors(); }
            Vector getEigenValues() const {Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> solver(matrix_);return solver.eigenvalues(); }
            Matrix getCholeskyDecomposition() const {Eigen::LLT<Eigen::MatrixXd> llt(matrix_);return llt.matrixL();}
        };

        struct CovarianceMatrix {
 
            CorrelationMatrix correlationMatrix_; 
            Vector variances_;

            CovarianceMatrix(const Matrix& matrix) {

                variances_ = matrix.diagonal();
                for (double v:variances_) { if (v<=0.0) errors::throwInvalidInputr("Variance value must be positive");}
                int n = variances_.size();
                Eigen::MatrixXd corrMatrix = Eigen::MatrixXd::Identity(n,n); 
                double corrCoeff;
                for (int i = 0; i < n; ++i) {
                    for (int j = i + 1; j < n; ++j) {
                        corrCoeff = matrix(i,j) / ( std::sqrt(variances_[i])*std::sqrt(variances_[j]));
                        corrMatrix(i, j) = corrCoeff; 
                        corrMatrix(j, i) = corrCoeff; 
                    }
                }

                correlationMatrix_ = CorrelationMatrix(corrMatrix);
            }
            CovarianceMatrix(const CorrelationMatrix& corrMatrix, const Vector& variances) {

                if (corrMatrix.matrix_.rows() != variances.size())  errors::throwInvalidInputr("The length of the vector of variances must equal the size of the colleration matrix");
                for (double v:variances) { if (v<=0.0) errors::throwInvalidInputr("Variance value must be positive");}
                variances_ = variances;
                correlationMatrix_ = corrMatrix;
            }

            Matrix get() const {

                int n = correlationMatrix_.matrix_.rows();
                Eigen::MatrixXd matrix = variances_.asDiagonal();
                double cov; 
                for (int i = 0; i < n; ++i) {
                    for (int j = i + 1; j < n; ++j) {
                        cov = getCovariance(i,j);
                        matrix(i, j) = cov; 
                        matrix(j, i) = cov; 
                    }
                }

                return matrix;
            }
            double getCovariance(int i, int j) const {return correlationMatrix_.getCorrelation(i,j)*getStandardDeviation(i)*getStandardDeviation(j);}
            double getVariance(int i) const {return variances_[i]; }
            double getStandardDeviation(int i) const{return std::sqrt(getVariance(i));}

        };
        

    }
}