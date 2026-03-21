#pragma once 
#include <numcpp/objects/base.hpp>

namespace numcpp::optim::bfgs {

    inline void bfgsUpdate(objects::Matrix& H, const objects::Vector& s, const objects::Vector& y) {
        double ys = y.dot(s);
        if (ys <= 0.0) return;   // curvature condition violated — skip update
    
        double rho = 1.0 / ys;
        int    n   = s.size();
        objects::Matrix    I   = objects::Matrix::Identity(n, n);
    
        objects::Matrix A = I - rho * s * y.transpose();
        H     = A * H * A.transpose() + rho * s * s.transpose();
    }

    struct LineSearchResult {
        double alpha;
        double f;
        objects::Vector g;
        bool ok;
    };
    
    inline LineSearchResult lineSearchStrongWolf(
        const std::function<double(const objects::Vector&)>& f,
        const std::function<objects::Vector(const objects::Vector&)>& grad,
        const objects::Vector& x,
        const objects::Vector& p,
        double f0,
        const objects::Vector&  g0,
        double c1       = 1e-4,
        double c2       = 0.9,
        double alpha_max = 1.0,
        int max_iter  = 40) {

        double dphi0 = g0.dot(p);
        assert(dphi0 < 0);
    
        auto phi  = [&](double a) { return f(x + a * p); };
        auto dphi = [&](double a, objects::Vector& g_a){ g_a = grad(x + a * p); return g_a.dot(p); };
    
        auto zoom = [&](double a_lo, double phi_lo,
                        double a_hi, double /*phi_hi*/) -> LineSearchResult {
            objects::Vector g_trial;
            for (int i = 0; i < max_iter; ++i) {
                double a_t   = 0.5 * (a_lo + a_hi);   // bisection (could use cubic interp)
                double phi_t = phi(a_t);
    
                if (phi_t > f0 + c1 * a_t * dphi0 || phi_t >= phi_lo) {
                    a_hi = a_t;
                } else {
                    double dphi_t = dphi(a_t, g_trial);
                    if (std::abs(dphi_t) <= -c2 * dphi0)
                        return {a_t, phi_t, g_trial, true};
                    if (dphi_t * (a_hi - a_lo) >= 0) a_hi = a_lo;
                    a_lo   = a_t;
                    phi_lo = phi_t;
                }
            }
            // return best found so far
            objects::Vector g_a; dphi(a_lo, g_a);
            return {a_lo, phi(a_lo), g_a, false};
        };
    
        double a_prev = 0.0, phi_prev = f0;
        double a_curr = alpha_max;
    
        for (int i = 1; i <= max_iter; ++i) {
            double phi_curr = phi(a_curr);
    
            if (phi_curr > f0 + c1 * a_curr * dphi0 || (i > 1 && phi_curr >= phi_prev))
                return zoom(a_prev, phi_prev, a_curr, phi_curr);
    
            objects::Vector g_curr;
            double dphi_curr = dphi(a_curr, g_curr);
    
            if (std::abs(dphi_curr) <= -c2 * dphi0)
                return {a_curr, phi_curr, g_curr, true};
    
            if (dphi_curr >= 0)
                return zoom(a_curr, phi_curr, a_prev, phi_prev);
    
            a_prev   = a_curr;
            phi_prev = phi_curr;
            a_curr   = std::min(2.0 * a_curr, alpha_max);
        }
    
        objects::Vector g_a; dphi(a_curr, g_a);
        return {a_curr, phi(a_curr), g_a, false};
    }

    struct BFGSResult {
        objects::Vector x;
        double f;
        objects::Vector g;
        int iter;
        bool converged;
    };
    
    inline BFGSResult bfgs(
        const std::function<double(const objects::Vector&)>& f,
        const std::function<objects::Vector(const objects::Vector&)>& grad,
        objects::Vector x, 
        double epsilon = 1e-6,
        int maxIter = 1000) {

        int n = x.size();
        objects::Vector g = grad(x);
        double fx = f(x);
        objects::Matrix H = objects::Matrix::Identity(n, n);  
    
        for (int k = 0; k < maxIter; ++k) {
            double gnorm = g.norm();
    
            if (gnorm < epsilon)
                return {x, fx, g, k, true};
    
            objects::Vector p = -(H * g);
    
            auto ls = lineSearchStrongWolf(f, grad, x, p, fx, g);
            if (!ls.ok) {
                H  = objects::Matrix::Identity(n, n);
                p  = -g;
                ls = lineSearchStrongWolf(f, grad, x, p, fx, g);
            }
    
            objects::Vector s = ls.alpha * p;   // step
            objects::Vector y = ls.g - g;       // gradient change
    
            bfgsUpdate(H, s, y);
    
            x  = x + s;
            g  = ls.g;
            fx = ls.f;
        }
    
        return {x, fx, g, maxIter, false};
    }
    
        
 
    
}