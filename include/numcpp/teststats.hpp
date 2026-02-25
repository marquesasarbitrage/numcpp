#pragma once 
#include "numcpp/objects.hpp"

namespace numcpp {

    namespace teststats {

        inline double oneSampleKolmogorovSmirnovProbabilityValue(const objects::Vector& sample, const std::function<double(double)>& cdfFunction, int sumIter = 100) {

            int n = sample.size();
            objects::Vector sorted = sample;
            std::sort(sorted.begin(), sorted.end());

            double D = 0.0;

            for (int i = 0; i < n; ++i) {
                double F = cdfFunction(sorted[i]);

                double Fn_upper = (i + 1) / (double)n;
                double Fn_lower = i / (double)n;

                double diff1 = std::abs(Fn_upper - F);
                double diff2 = std::abs(F - Fn_lower);

                D = std::max(D, std::max(diff1, diff2));
            }

            double sum = 0.0;
            for (int k = 1; k <= sumIter; ++k) {
                double term = std::pow(-1.0, k - 1) * std::exp(-2.0 * k * k * D * D * n);
                sum += term;
            }
            return 2.0 * sum;
        };

       
    }
}