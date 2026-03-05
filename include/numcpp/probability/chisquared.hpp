#pragma once 
#include "numcpp/objects/base.hpp"
#include <random>
#include <numcpp/probability/gamma.hpp>

namespace numcpp {

    namespace probability {


        struct ChiSquared{double nu; Gamma gamma() const {return Gamma{nu/2.0,2.0};}};

        inline double pdf(ChiSquared params, double x) {return pdf(params.gamma(),x);}

        inline double cdf(ChiSquared params, double x) {return cdf(params.gamma(),x);}

        inline double sample(ChiSquared params, std::mt19937& gen) { return sample(params.gamma(),gen); }
        
        inline objects::TheoricalMoments theoricalMoments(ChiSquared params) {return theoricalMoments(params.gamma());}

        inline objects::Vector sample(ChiSquared params, std::mt19937& gen, int n) { return sample(params.gamma(),gen,n); };
        

    }
}