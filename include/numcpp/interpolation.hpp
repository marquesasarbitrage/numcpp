#pragma once 
#include <map>
#include "numcpp/errors.hpp"


namespace numcpp {

    namespace interpolation {

        class Interpolation {
            public:
                Interpolation(const std::map<double, double>& data) {

                    if (data.size() < 2) errors::throwInvalidInputr("Interpolation object must be constrcuted from more than 2 points");
                    for(auto const& imap: data){xVector_.push_back(imap.first);yVector_.push_back(imap.second);}    
                }

                virtual ~Interpolation() = default;

                double evaluate(double x) const {

                    if (x < getLowerBoundX() || x > getUpperBoundX()) errors::throwInvalidInputr("Interpolation object cannot extrapolate");
                    return _evaluate(x);
                }

                double evaluateFirstDerivative(double x) const{

                    if (x < getLowerBoundX() || x > getUpperBoundX()) errors::throwInvalidInputr("Interpolation object cannot extrapolate");
                    return _evaluateFirstDerivative(x);
                }

                double evaluateSecondDerivative(double x) const{

                    if (x < getLowerBoundX() || x > getUpperBoundX()) errors::throwInvalidInputr("Interpolation object cannot extrapolate");
                    return _evaluateSecondDerivative(x);
                }

                double getLowerBoundX() const { return xVector_[0];}
                double getUpperBoundX() const { return xVector_[xVector_.size()-1];}
                const std::vector<double>& getVectorX() const { return xVector_; }
                const std::vector<double>& getVectorY() const  { return yVector_; }

            protected: 
                virtual double _evaluate(double x) const = 0;
                virtual double _evaluateFirstDerivative(double x) const = 0; 
                virtual double _evaluateSecondDerivative(double x) const = 0; 

                int findIndex(double x) const {

                    if (x==getLowerBoundX()) return 1; if (x==getUpperBoundX()) return xVector_.size();
                    return std::distance(xVector_.begin(), std::lower_bound(xVector_.begin(), xVector_.end(), x));
                }
                double getX(int i) const { return xVector_[i];}
                double getY(int i) const { return yVector_[i];}

            private:
                std::vector<double> xVector_; 
                std::vector<double> yVector_; 

        };

        class LinearInterpolation final: public Interpolation {
            public: 
                LinearInterpolation(const std::map<double, double>& data): Interpolation(data){}
                ~LinearInterpolation() = default;
            
            protected: 
                double _evaluate(double x) const override {
                    int i = findIndex(x);
                    double x0 = getX(i - 1), x1 = getX(i);
                    double y0 = getY(i - 1), y1 = getY(i);
                    double t = (x - x0) / (x1 - x0);
                    return y0 + t * (y1 - y0);
                }

                double _evaluateFirstDerivative(double x) const override {

                    int i = findIndex(x);
                    double x0 = getX(i - 1), x1 = getX(i);
                    double y0 = getY(i - 1), y1 = getY(i);
                    return (y1 - y0) / (x1 - x0);
                }

                double _evaluateSecondDerivative(double x) const override { return 0.0;}
        };

        class CubicSpline final: public Interpolation {
            public: 
                CubicSpline(const std::map<double, double>& data): Interpolation(data) {

                    x_ = getVectorX();
                    int n = x_.size() - 1; 
                    a_ = getVectorY();
                    
                    std::vector<double> h(n);
                    std::vector<double> alpha(n, 0.0);

                    for (int i = 0; i < n; ++i) {
                        h[i] = x_[i + 1] - x_[i];
                    }

                    for (int i = 1; i < n; ++i) {
                        alpha[i] = (3.0 / h[i]) * (a_[i + 1] - a_[i]) - (3.0 / h[i - 1]) * (a_[i] - a_[i - 1]);
                    }

                    std::vector<double> l(n + 1, 0.0);
                    std::vector<double> mu(n, 0.0);
                    std::vector<double> z(n + 1, 0.0);
                    
                    l[0] = 1.0;  
                    mu[0] = 0.0;
                    z[0] = 0.0;

                    for (int i = 1; i < n; ++i) {
                        l[i] = 2.0 * (x_[i + 1] - x_[i - 1]) - h[i - 1] * mu[i - 1];
                        mu[i] = h[i] / l[i];
                        z[i] = (alpha[i] - h[i - 1] * z[i - 1]) / l[i];
                    }

                    l[n] = 1.0;
                    z[n] = 0.0;
                    
                    c_.resize(n + 1, 0.0);
                    b_.resize(n, 0.0);
                    d_.resize(n, 0.0);
                    
                    for (int j = n - 1; j >= 0; --j) {
                        c_[j] = z[j] - mu[j] * c_[j + 1];
                        b_[j] = (a_[j + 1] - a_[j]) / h[j] - h[j] * (c_[j + 1] + 2.0 * c_[j]) / 3.0;
                        d_[j] = (c_[j + 1] - c_[j]) / (3.0 * h[j]);
                    }
                }
                ~CubicSpline() = default;
                
            protected: 
                double _evaluate(double x) const override {

                    int i = findIndex(x);
                    double dx = x - getX(i-1);
                    return a_[i-1] + b_[i-1] * dx + c_[i-1] * dx * dx + d_[i-1] * dx * dx * dx;
                }

                double _evaluateFirstDerivative(double x) const override {

                    int i = findIndex(x);
                    double dx = x - getX(i-1);
                    return b_[i-1] + 2.0*c_[i-1] * dx + 3.0*d_[i-1] * dx * dx;
                }

                double _evaluateSecondDerivative(double x) const override {

                    int i = findIndex(x);
                    double dx = x - getX(i-1);
                    return 2.0*c_[i-1] + 6.0*d_[i-1] * dx;
                }

            private:
                std::vector<double> a_, b_, c_, d_, x_;
            
        };
    
    }
}