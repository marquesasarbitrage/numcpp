#pragma once 
#include "numcpp/objects.hpp"

namespace numcpp {

    namespace multstats {

        inline objects::CovarianceMatrix sampleCovarianceMatrix(const objects::Matrix& X) {

            Eigen::MatrixXd centered = X.colwise() - X.rowwise().mean();
            return objects::CovarianceMatrix((centered * centered.transpose()) / (X.cols() - 1.0));
        }

        inline objects::CovarianceMatrix populationCovarianceMatrix(const objects::Matrix& X) {

            Eigen::MatrixXd centered = X.colwise() - X.rowwise().mean();
            return objects::CovarianceMatrix((centered * centered.transpose()) / X.cols() );
        }

        inline objects::CorrelationMatrix correlationMatrix(const objects::Matrix& X) {

            return sampleCovarianceMatrix(X).correlationMatrix_;
        }
    }
}