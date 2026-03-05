#include <numcpp/stats.hpp>
#include <numcpp/toolbox.hpp>
#include <cassert>
#include <iostream>
#include "utils.hpp"
#include <cassert>

Eigen::VectorXd getRawData() {

    std::vector<double> data =  {16, 14, 64, 88, 53, 3, 54, 96, 23, 88, 52,
        19, 92, 36, 94, 64, 69, 19, 92, 23, 96, 47,
        63, 17, 80, 59, 99, 7, 75, 23, 62, 96, 11, 78,
        70, 19, 32, 51, 44, 96, 90, 58, 64, 75, 20, 53,
        78, 95, 4, 82};

    return numcpp::toolbox::getVectorObject(data);

}

void testStats() {

    Eigen::VectorXd data = getRawData();

    assert(numcpp::stats::mean(data)==56.06);
    assert(isClose(numcpp::stats::sampleVariance(data),925.8127,1e-4));
    assert(isClose(numcpp::stats::populationVariance(data),907.2964,1e-4));
    assert(isClose(numcpp::stats::sampleStandardDeviation(data),30.4272,1e-4));
    assert(isClose(numcpp::stats::populationStandardDeviation(data),30.1214,1e-4));
    assert(isClose(numcpp::stats::sampleSkewness(data),-0.2365,1e-4));
    assert(isClose(numcpp::stats::populationSkewness(data),-0.2294,1e-4));
    assert(isClose(numcpp::stats::sampleExcessKurtosis(data),-1.28317,1e-4));
    assert(isClose(numcpp::stats::populationExcessKurtosis(data),-1.2760,1e-4));

}

int main() {

    testStats(); 
    std::cout << "All tests for the stat module have been passed successfully!" << std::endl;

    return 0;
}

