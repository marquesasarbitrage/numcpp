#pragma once 
#include <vector>
#include <numcpp/objects/base.hpp>

namespace numccp {

    namespace regression {

        namespace olstools {

            inline numcpp::objects::Matrix addIntercept(const numcpp::objects::Matrix& X, bool fitIntercept) {

                if (!fitIntercept) return X;
                Eigen::MatrixXd XX(X.rows(), X.cols() + 1);
                XX.col(0) = Eigen::VectorXd::Ones(X.rows()); // first column = ones
                XX.block(0, 1, X.rows(), X.cols()) = X;
                return XX;
            }

            inline numcpp::objects::Matrix projectionMatrix(const numcpp::objects::Matrix& X) {

                Eigen::MatrixXd XtX = X.transpose() * X;
                Eigen::LDLT<Eigen::MatrixXd> ldlt(XtX);
                return X * ldlt.solve(X.transpose());
            }

        }

        struct OrdinaryLeastSquare {
            bool useIntercept; 
            std::vector<double> betas; 
            numcpp::objects::Vector residuals; 
            double rSquared; 
            std::vector<double> standardError; 
        
            std::vector<double> confidenceInterval(bool upper);
        };

        

        //inline OrdinaryLeastSquare ols(const numcpp::objects::Vector& Y, const numcpp::objects::Vector& X, bool fitIntercept);
    }
}