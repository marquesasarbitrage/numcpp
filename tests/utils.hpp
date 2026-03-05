#pragma once 
#include <cmath>
#include <Eigen/Dense>
#include <cassert>

inline bool isClose(double a, double b, double eps) { return std::abs(a-b)<eps;}

inline void assetIsMatrixesSimilar(const Eigen::MatrixXd& A, const Eigen::MatrixXd& B, double tol = 1e-9) {

    assert(A.rows() == B.rows() && A.cols() == B.cols());
    assert(A.isApprox(B, tol));
}
