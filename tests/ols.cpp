#include <numcpp/regression/ols.hpp>
#include <iostream>
#include "numcpp/objects/base.hpp"
#include "utils.hpp"


void testHeightsWeightsExampleWiki() {

    numcpp::objects::Vector y(15); 
    y << 52.21,53.12,54.48,55.84,57.20,58.57,59.93,61.29,63.11,64.47,66.28,68.10,69.92,72.19,74.46;

    numcpp::objects::Matrix X(15,2);
    double heights[] = {1.47,1.50,1.52,1.55,1.57,1.6,1.63,1.65,1.68,1.70,1.73,1.75,1.78,1.80,1.83};

    for(int i=0;i<15;i++){
        X(i,0) = heights[i];
        X(i,1) = heights[i]*heights[i];
    }

    numcpp::regression::OLS ols = numcpp::regression::ols(y, X, true);

    
    assert(isClose(ols.betas(0), 128.8128, 1e-4));
    assert(isClose(ols.betas(1), -143.1620, 1e-4));
    assert(isClose(ols.betas(2), 61.9603, 1e-4));

    numcpp::objects::Vector se = ols.standardErrors();
    assert(isClose(se(0), 16.3083, 1e-4));
    assert(isClose(se(1), 19.8332, 1e-4));
    assert(isClose(se(2), 6.0084, 1e-4));

    numcpp::objects::Vector ts = ols.tstats();
    assert(isClose(ts(0), 7.8986, 1e-4));
    assert(isClose(ts(1), -7.2183, 1e-4));
    assert(isClose(ts(2), 10.3122, 1e-4));

    numcpp::objects::Vector pv = ols.pValues();
    assert(isClose(pv(0), 0.0, 1e-4));
    assert(isClose(pv(1), 0.0, 1e-4));
    assert(isClose(pv(2), 0.0, 1e-4));

    assert(isClose(ols.rSquared,0.9989,1e-4));
    assert(isClose(ols.adjustedRSquared(),0.9987,1e-4));

}


int main() {

    testHeightsWeightsExampleWiki(); 
    std::cout << "All tests for the OLS object have been passed successfully!" << std::endl;
    return 0;
}