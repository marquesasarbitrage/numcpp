#pragma once
#include "numcpp/objects.hpp"
#include "numcpp/constants.hpp"

namespace numcpp {

    namespace stats {

        inline double mean(const objects::Vector& data) { return data.size()>0 ? data.mean() : constants::DOUBLE_NAN;}

        inline double sampleVariance(const objects::Vector& data) {

            const Eigen::Index n = data.size();
            if (n == 0 ||  n < 2) return constants::DOUBLE_NAN;

            const double mean = data.mean();
            const double sq = (data.array() - mean).square().sum();

            return sq / (n - 1.0);
        }

        inline double populationVariance(const objects::Vector& data) {

            const Eigen::Index n = data.size();
            if (n == 0) return constants::DOUBLE_NAN;

            const double mean = data.mean();
            const double sq = (data.array() - mean).square().sum();

            return sq / n;
        }

        inline double sampleStandardDeviation(const objects::Vector& data) {return std::sqrt(sampleVariance(data));}

        inline double populationStandardDeviation(const objects::Vector& data) {return std::sqrt(populationVariance(data));}

        inline double populationSkewness(const objects::Vector& data) {
            const Eigen::Index n = data.size();
            if (n == 0) return constants::DOUBLE_NAN;

            const Eigen::ArrayXd centered = data.array() - mean(data);

            const double m2 = centered.square().mean();
            if (m2 == 0.0) return 0.0;

            return centered.pow(3).mean() / std::pow(m2, 1.5); 
        }

        inline double sampleSkewness(const objects::Vector& data) {

            const Eigen::Index n = data.size();
            if (n == 0 ||  n < 3) return constants::DOUBLE_NAN;
            return std::sqrt(n * (n - 1.0)) / (n - 2.0)*populationSkewness(data);
        }

        inline double populationExcessKurtosis(const objects::Vector& data) {

            const Eigen::Index n = data.size();
            if (n == 0) return constants::DOUBLE_NAN;

            const double mean = data.mean();
            const Eigen::ArrayXd centered = data.array() - mean;

            const double m2 = centered.square().mean();
            if (m2 == 0.0) return -3.0;
            return centered.pow(4).mean() / (m2 * m2) - 3.0;  
        }

        inline double sampleExcessKurtosis(const objects::Vector& data) {

            const Eigen::Index n = data.size();
            if (n == 0 ||  n < 4) return constants::DOUBLE_NAN;
            const double n1 = (n - 1.0);
            const double n2 = (n - 2.0);
            const double n3 = (n - 3.0);

            return ((n - 1.0) / ((n - 2.0) * (n - 3.0))) *
                ((n + 1.0) * populationExcessKurtosis(data) + 6.0);
        }
    }
}