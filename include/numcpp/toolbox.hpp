#pragma once 
#include <cassert>
#include <numcpp/objects/base.hpp>
#include <vector>

namespace numcpp {

    namespace toolbox {

        inline objects::Vector getVectorObject(const std::vector<double>& data) {

            return Eigen::Map<const Eigen::VectorXd>(data.data(), data.size());
        }

        inline objects::Matrix getMatrixObject(const std::vector<std::vector<double>>& data) {

            int rows = data.size();
            int cols = data[0].size();

            assert(rows != 0);
            assert(cols != 0);

            Eigen::MatrixXd mat(rows, cols);

            for (int i = 0; i < rows; ++i) {
                assert(data[i].size() == cols);
                for (int j = 0; j < cols; ++j) {mat(i, j) = data[i][j];}
            }
            return mat;
        }
    }
}
