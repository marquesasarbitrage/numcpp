#include <numcpp/probability.hpp>
#include <numcpp/objects.hpp>
#include <numcpp/teststats.hpp>
#include <cassert>
#include <iostream>
#include "utils.hpp"
#include <cassert>

void templateTestSampler(const std::function<numcpp::objects::Vector(int,std::mt19937&)>& sampler, const std::function<double(double)> cdf) {

    std::random_device rd;      
    std::mt19937 gen(rd());
    int failures = 0;
    int runs = 1000;

    for(int i = 0; i < runs; ++i)
    {
        numcpp::objects::Vector sample = sampler(1000,gen);
        double p = numcpp::teststats::oneSampleKolmogorovSmirnovProbabilityValue(sample,cdf);
        if(p < 0.05)failures++;
    }

    assert(failures < runs * 0.1);
}

void testNormalDistributions() {

    using namespace numcpp::probability;
    Normal n{}; 

    assert(std::abs(cdf(n,0.0) - 0.5) < 1e-6); 
    assert(std::abs(cdf(n,1.0) - 0.841344746) < 1e-6);
    assert(std::abs(cdf(n,-1.0) - 0.158655254) < 1e-6);
    assert(std::abs(cdf(n,2.0) - 0.977249868) < 1e-6);
    assert(std::abs(cdf(n,-2.0) - 0.022750132) < 1e-6);
    assert(std::abs(cdf(n,3.0) - 0.998650102) < 1e-6);

    assert(std::abs(invCdf(n, 0.5) - 0.0) < 1e-6); 
    assert(std::abs(invCdf(n, .841344746) - 1.0) < 1e-6);
    assert(std::abs(invCdf(n, .158655254) + 1.0) < 1e-6);
    assert(std::abs(invCdf(n, 0.977249868) - 2.0) < 1e-6);
    assert(std::abs(invCdf(n, 0.022750132) + 2.0) < 1e-6);
    assert(std::abs(invCdf(n, 0.998650102) - 3.0) < 1e-6);

    assert(std::abs(cf(n,0.0).real() - 1.0) < 1e-12);
    assert(std::abs(cf(n,0.0).imag()) < 1e-12);

    assert(std::abs(cf(n,0.1).real() - .995012479) < 1e-9);
    assert(std::abs(cf(n,0.1).imag()) < 1e-12);

    assert(std::abs(cf(n,1.0).real() - .6065306597) < 1e-9);
    assert(std::abs(cf(n,1.0).imag()) < 1e-12);

    assert(std::abs(pdf(n,0.0) - numcpp::constants::ONE_OVER_SQRT_TWO_PI) < 1e-6); // pdf(0) for standard normal

    n = {2.0,3.0};

    assert(std::abs(cdf(n,2.0) - 0.5) < 1e-6);
    assert(std::abs(cdf(n,5.0) - 0.841344746) < 1e-6);
    assert(std::abs(cdf(n,-1.0) - 0.158655254) < 1e-6);

    assert(std::abs(invCdf(n,.5) - 2.0) < 1e-6);
    assert(std::abs(invCdf(n,.841344746) - 5.0) < 1e-6);
    assert(std::abs(invCdf(n,.158655254) + 1.0) < 1e-6);
   
    assert(std::abs(cf(n,0.0).real() - 1.0) < 1e-12);
    assert(std::abs(cf(n,0.0).imag()) < 1e-12);

    assert(std::abs(cf(n,0.1).real() - std::exp(-0.5 * 9.0 * 0.1 * 0.1) * std::cos(2.0 * 0.1)) < 1e-9);
    assert(std::abs(cf(n,0.1).imag() - std::exp(-0.5 * 9.0 * 0.1 * 0.1) * std::sin(2.0 * 0.1)) < 1e-9);

    assert(std::abs(pdf(n,2.0) - 1.0 / (3.0 * std::sqrt(2 * M_PI))) < 1e-6);

    templateTestSampler([n](int v ,std::mt19937& gen){return numcpp::probability::sample(n,v,gen);}, [n](double x) { return cdf(n,x); });
}

void testGammaDistribution() {

    using namespace numcpp::probability;
    Gamma g{2.0,2.0};
    
    assert(isClose(pdf(g,2.0),0.1839397206, 1e-9));
    assert(isClose(pdf(g,7.0),0.05284542099, 1e-9));
    assert(isClose(pdf(g,14.0),0.003191586879, 1e-9));

    assert(isClose(cdf(g,2.0),0.2642411177, 1e-9));
    assert(isClose(cdf(g,7.0),0.8641117746, 1e-9));
    assert(isClose(cdf(g,14.0),0.9927049443, 1e-9));

    templateTestSampler([g](int n ,std::mt19937& gen){return numcpp::probability::sample(g,n,gen);}, [g](double x) { return cdf(g,x); });
    
}

void testChiSquareDistribution() {

    using namespace numcpp::probability;
    ChiSquared g{2.0};
    
    assert(isClose(pdf(g,2.0),0.1839397206, 1e-9));
    assert(isClose(pdf(g,7.0),0.01509869171, 1e-9));
    assert(isClose(pdf(g,14.0),0.0004559409828, 1e-9));

    assert(isClose(cdf(g,2.0),0.6321205588, 1e-9));
    assert(isClose(cdf(g,7.0),0.9698026166, 1e-9));
    assert(isClose(cdf(g,14.0),0.999088118, 1e-9));

    templateTestSampler([g](int n ,std::mt19937& gen){return numcpp::probability::sample(g,n,gen);}, [g](double x) { return cdf(g,x); });

}

void testStudent() {

    using namespace numcpp::probability;
    Student g{2.0};
    
    assert(isClose(pdf(g,-8.0),0.001865022591, 1e-9));
    assert(isClose(pdf(g,-2.0),0.06804138174, 1e-9));
    assert(isClose(pdf(g,.0),0.3535533906, 1e-9));
    assert(isClose(pdf(g,3.0),0.02741012223, 1e-9));
    assert(isClose(pdf(g,9.0),0.001322460964, 1e-9));

    assert(isClose(cdf(g,-8.0),0.007634036083, 1e-9));
    assert(isClose(cdf(g,-2.0),0.09175170954, 1e-9));
    assert(isClose(cdf(g,.0),0.5, 1e-9));
    assert(isClose(cdf(g,3.0),0.9522670169, 1e-9));
    assert(isClose(cdf(g,9.0),0.99393917, 1e-9));

    assert(isClose(invCdf(g,0.0000001),-2236.067642, 1e-3));
    assert(isClose(invCdf(g,.01),-6.964556734, 1e-9));
    assert(isClose(invCdf(g,.2),-1.060660172, 1e-9));
    assert(isClose(invCdf(g,0.9999),70.70007107, 1e-4));
    assert(isClose(invCdf(g,0.9999999),2236.067642, 1e-3));

    templateTestSampler([g](int n ,std::mt19937& gen){return numcpp::probability::sample(g,n,gen);}, [g](double x) { return cdf(g,x); });
}

int main() {

    testNormalDistributions();
    testGammaDistribution();
    testChiSquareDistribution();
    testStudent();
    std::cout << "All tests for probability module have been passed!" << std::endl;
    return 0;
}