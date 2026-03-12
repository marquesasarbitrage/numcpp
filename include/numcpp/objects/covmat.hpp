#pragma once 
#include <cassert>
#include <numcpp/objects/base.hpp>

namespace numcpp {

    namespace objects {

        struct CovarianceMatrix;

        struct CorrelationMatrix {

            Matrix matrix_; 
            CorrelationMatrix(): matrix_(Eigen::MatrixXd::Identity(1, 1)) {}

            CorrelationMatrix(double coeff) {

                Matrix mat = Eigen::MatrixXd::Zero(2,2);
                mat << 1.0, coeff, 
                       coeff, 1.0;
                assert(isValid(mat)); 
                matrix_ = mat;
            }

            CorrelationMatrix(const Vector& coeffs) {

                assert(isNumberCoefficientsValid(coeffs.size()));
                size_t c = coeffs.size();
                int p = static_cast<int>(1.0 + std::sqrt(1.0 + 8.0 * coeffs.size())) / 2.0;

                Matrix mat = Eigen::MatrixXd::Identity(p, p);

                int idx = 0;
                for (int i = 0; i < p; ++i) {
                    for (int j = i + 1; j < p; ++j) {
                        mat(i, j) = coeffs[idx];
                        mat(j, i) = coeffs[idx]; 
                        ++idx;
                    }
                }

                assert(isValid(mat)); 
                matrix_ = mat;

            }

            CorrelationMatrix(const Matrix& matrix) {assert(isValid(matrix));matrix_ = matrix;}

            double getCorrelation(int i, int j) const { return matrix_(i,j); }

            Matrix getEigenVectors() const {Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> solver(matrix_);return solver.eigenvectors(); }

            Vector getEigenValues() const {Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> solver(matrix_);return solver.eigenvalues(); }

            Matrix getCholeskyDecomposition() const {Eigen::LLT<Eigen::MatrixXd> llt(matrix_);return llt.matrixL();}

            static bool isValid(const Matrix& matrix) {

                if (matrix.rows() != matrix.cols()) return false;
                if (!((matrix - matrix.transpose()).norm() < 1e-10)) return false;
                Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> solver(matrix);
                if (!((solver.eigenvalues().array() > 0).all())) return false;
                for (double v: matrix.diagonal()) if (std::abs(v-1.0)>1e-10) return false;
                for (size_t i = 0; i < matrix.rows(); ++i) {
                    for (size_t j = i + 1; j < matrix.rows(); ++j) {
                        if (std::abs(matrix(i, j))-1.0 > 1e-20) return false;
                    }
                }
                return true;
            }

            static bool isNumberCoefficientsValid(int n) {
                int p = static_cast<int>(1.0 + std::sqrt(1.0 + 8.0 * n)) / 2.0;
                return !(p * (p - 1) / 2 != n);
            }

        };

        struct CovarianceMatrix {
 
            CorrelationMatrix correlationMatrix_; 
            Vector variances_;

            CovarianceMatrix(const Matrix& matrix) {

                assert(isVarianceValid(matrix.diagonal()));
                variances_ = matrix.diagonal();
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
                CorrelationMatrix::isValid(corrMatrix);
                correlationMatrix_ = CorrelationMatrix(corrMatrix);
            }

            CovarianceMatrix(const CorrelationMatrix& corrMatrix, const Vector& variances) {

                assert(isCorrMatrixVarianceVectorMatch(corrMatrix.matrix_,variances));
                assert(isVarianceValid(variances));
                variances_ = variances;
                correlationMatrix_ = corrMatrix;
            }

            static bool isVarianceValid(const Vector& variances) {for (double v:variances) { if (v<=1e-20) {return false;}}return true;}

            static bool isCorrMatrixVarianceVectorMatch(const Matrix& corrMatrix, const Vector& variances) {return corrMatrix.rows() == variances.size();}

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

            double getStandardDeviation(int i) const{return std::sqrt(variances_[i]);}

        };
    }
}