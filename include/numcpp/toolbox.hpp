#pragma once 
#include <numcpp/objects/base.hpp>
#include <vector>
#include <numcpp/errors.hpp>

namespace numcpp {

    namespace toolbox {

        inline objects::Vector getVectorObject(const std::vector<double>& data) {

            return Eigen::Map<const Eigen::VectorXd>(data.data(), data.size());
        }

        inline objects::Matrix getMatrixObject(const std::vector<std::vector<double>>& data) {

            int rows = data.size();
            int cols = data[0].size();

            if (rows == 0) numcpp::errors::throwInvalidInputr("The data sizing is invalid to construct a matrix, number of rows must be non-null");
            if (cols == 0) numcpp::errors::throwInvalidInputr("The data sizing is invalid to construct a matrix, number of columns must be non-null");

            // Create Eigen matrix
            Eigen::MatrixXd mat(rows, cols);

            // Fill Eigen matrix
            for (int i = 0; i < rows; ++i) {
                if (data[i].size() != cols) numcpp::errors::throwInvalidInputr("The data sizing is invalid to construct a matrix, The number of columns must be the same for each row");
                for (int j = 0; j < cols; ++j) {
                    mat(i, j) = data[i][j];
                }
            }
            return mat;
        }
    }
}
