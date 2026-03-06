#pragma once 
#include "Eigen/Core"
#include <numcpp/objects/base.hpp>
#include <numcpp/errors.hpp>
#include <numcpp/probability/tstudent.hpp>

namespace numcpp {

    namespace regression {

        namespace olstools {

            inline void checkInput(const numcpp::objects::Vector& Y,const numcpp::objects::Matrix& X) {

                if(Y.size() != X.rows()) errors::throwInvalidInputr("The size of the vector Y must be equal to the number of rows from matrix X in OLS");
            }

            inline numcpp::objects::Matrix addIntercept(const numcpp::objects::Matrix& X, bool fitIntercept) {

                if (!fitIntercept) return X;
                Eigen::MatrixXd XX(X.rows(), X.cols() + 1);
                XX.col(0) = Eigen::VectorXd::Ones(X.rows());
                XX.block(0, 1, X.rows(), X.cols()) = X;
                return XX;
            }

            inline numcpp::objects::Matrix inverseXtX(const numcpp::objects::Matrix& X) {

                Eigen::MatrixXd XtX = X.transpose() * X;
                Eigen::LDLT<Eigen::MatrixXd> ldlt(XtX);
                return ldlt.solve(Eigen::MatrixXd::Identity(XtX.rows(), XtX.rows()));
            }

            inline numcpp::objects::Matrix projectionMatrix(const numcpp::objects::Matrix& X) {

                Eigen::MatrixXd XtX = X.transpose() * X;
                Eigen::LDLT<Eigen::MatrixXd> ldlt(XtX);
                return X * ldlt.solve(X.transpose());
            }


        }

        struct OLS {
            bool useIntercept; 
            numcpp::objects::Vector betas; 
            numcpp::objects::Vector residuals; 
            numcpp::objects::Matrix inverseXtX;
            double sampleVariance;
            double rSquared;
        
            inline int numberParameters() const {return betas.size();}

            inline int numberObservations() const {return residuals.size();}

            inline int degreesFreedom() const {return numberObservations()-numberParameters();}

            inline Eigen::VectorXd standardErrors() const {return (sampleVariance * inverseXtX.diagonal()).array().sqrt();}

            inline Eigen::VectorXd confidenceInterval(double quantile, bool upper) {

                double tcrit = probability::invCdf(probability::Student{double(degreesFreedom())},0.5 + quantile / 2.0);
                Eigen::VectorXd ci(numberParameters());
                Eigen::VectorXd se = standardErrors();
                for(int j = 0; j < numberParameters(); ++j) {
                    ci(j) = betas(j) + (upper ? tcrit * se(j) : -tcrit * se(j));
                }

                return ci;
            }

            inline Eigen::VectorXd tstats() {

                Eigen::VectorXd ts(numberParameters()); 
                Eigen::VectorXd se = standardErrors();
                for(int j = 0; j < numberParameters(); ++j) {
                    ts(j) = betas(j)/se(j);
                }
                return ts;
            }

            inline double adjustedRSquared() const {return 1.0 - (1.0 - rSquared) * (numberObservations() - 1.0) / double(degreesFreedom());}

            inline Eigen::VectorXd pValues() {

                Eigen::VectorXd tstats_ = tstats();
                Eigen::VectorXd p(tstats_.size());

                for(int i = 0; i < tstats_.size(); ++i) {
                    p(i) = 2.0 * (1.0 - cdf(probability::Student{double(degreesFreedom())}, std::abs(tstats_(i))));
                }

                return p;
            }
        };


        inline OLS ols(const numcpp::objects::Vector& Y,const numcpp::objects::Matrix& X, bool fitIntercept) {

            numcpp::objects::Matrix X_ = olstools::addIntercept(X, fitIntercept); 
            numcpp::objects::Matrix XtXinv = olstools::inverseXtX(X_);
            numcpp::objects::Vector betas = XtXinv*X_.transpose()*Y;
            numcpp::objects::Vector residuals = Y - X_*betas;
            double rss = residuals.dot(residuals);
            double sampleVariance = rss / double(Y.size() - betas.size());
            double ymean = Y.mean();
            double tss = 0.0;
            for (size_t i = 0; i < Y.size(); ++i){double d = Y[i] - ymean;tss += d * d;}
            double rsquared = 1.0 - rss / tss;
            return {fitIntercept,betas,residuals,XtXinv,sampleVariance,rsquared};
        }


        

        //inline OrdinaryLeastSquare ols(const numcpp::objects::Vector& Y, const numcpp::objects::Vector& X, bool fitIntercept);
    }
}