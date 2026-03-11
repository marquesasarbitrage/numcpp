#include <numcpp/probability/stochprocess/abm.hpp>
#include "numcpp/objects/base.hpp"
#include "utils.hpp"
#include <iostream>

std::map<double,double> getDataForABM() {

    return {
        { 0.0,  3.75 / 100.0 },
        { 0.008333333333,  3.70 / 100.0 },
        { 0.01111111111,  3.66 / 100.0 },
        { 0.01388888889,  3.65 / 100.0 },
        { 0.01666666667,  3.64 / 100.0 },
        { 0.01944444444,  3.64 / 100.0 },
        { 0.02777777778, 3.64 / 100.0 },
        { 0.03055555556, 3.65 / 100.0 },
        { 0.03333333333, 3.64 / 100.0 },
        { 0.03611111111, 3.66 / 100.0 },
        { 0.03888888889, 3.65 / 100.0 },
        { 0.05, 3.64 / 100.0 },
        { 0.05277777778, 3.63 / 100.0 },
        { 0.05555555556, 3.64 / 100.0 },
        { 0.05833333333, 3.65 / 100.0 },
        { 0.06666666667, 3.66 / 100.0 },
        { 0.06944444444, 3.66 / 100.0 },
        { 0.07222222222, 3.64 / 100.0 },
        { 0.075, 3.65 / 100.0 },
        { 0.07777777778, 3.68 / 100.0 },
        { 0.08611111111,  3.69 / 100.0 }
    };      
}

void testTemplateABMSampling(const numcpp::probability::stochprocess::ArithmeticBrownianMotion& abm, int size, int steps, double T) {

    std::mt19937 gen(42);
    numcpp::objects::Matrix sample = numcpp::probability::stochprocess::samplePath(abm,T,steps,size,gen);
    double dt = T/steps;
    for (size_t i = 0; i<size; ++i) {
        double t = 0.0;
        std::map<double, double> genData;
        for (size_t j = 0; j<steps; ++j) {
            genData[t] = sample(i,j);
            t+=dt; 
        }

        numcpp::probability::stochprocess::ArithmeticBrownianMotion abmFit = numcpp::probability::stochprocess::fitArithmeticBrownianMotion(genData);
        assert(isClose(abm.mu, abmFit.mu, 1e-2));
        assert(isClose(abm.sigma, abmFit.sigma, 1e-2));
    } 
}

void testABM() {

    numcpp::probability::stochprocess::ArithmeticBrownianMotion abm = numcpp::probability::stochprocess::fitArithmeticBrownianMotion(getDataForABM());
    assert(abm.x0==0.0369); 
    assert(isClose(abm.mu, -0.006967741935, 1e-9));
    assert(isClose(abm.sigma, 0.003002421255, 1e-9));

    testTemplateABMSampling(abm,100,1000,2.0);

}

int main() {

    testABM(); 
    std::cout << "All tests for the stochastic processes have been passed successfully!" << std::endl;
}