#pragma once 
#include <vector>
#include "numcpp/objects.hpp"

namespace numccp {

    namespace regression {

        struct OrdinaryLeastSquare {
            bool useIntercept; 
            std::vector<double> betas; 
            numcpp::objects::Vector residuals; 
            double rSquared; 
            std::vector<double> standardError; 
        
            std::vector<double> confidenceInterval(bool upper);
        };

        inline OrdinaryLeastSquare ols(numcpp::objects::Vector Y, numcpp::objects::Vector X, bool fitIntercept);
    }
}