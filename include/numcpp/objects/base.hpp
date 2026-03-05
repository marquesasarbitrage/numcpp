#pragma once 
#include <Eigen/Dense>

namespace numcpp {

    namespace objects {

        using Vector = Eigen::VectorXd;
        using Matrix = Eigen::MatrixXd;

        struct TheoricalMoments {double mean,variance,skewness,excessKurtosis;};

    }
}