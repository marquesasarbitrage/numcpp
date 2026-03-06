#pragma once
#include <vector>
#include <cmath>
#include <functional>

namespace numcpp {

    namespace optim {

        namespace neldermeadtools {

            enum class SimplexInitializationMethod {BASIC, SCALED, SYMMETRIC};

            struct Vertex {
                std::vector<double> points_; 
                double value_;             
            }; 

            struct Simplex {    

                std::vector<Vertex> vertices_; 

                Simplex(const std::vector<Vertex>& vertices): vertices_(vertices) {}

                void sortVertices() {
                    std::sort(vertices_.begin(), vertices_.end(), [](const Vertex& a, const Vertex& b) {return a.value_ < b.value_;});
                }
                void setVertex(const Vertex& vertex_, int i) { vertices_[i] = vertex_; sortVertices(); }
                const Vertex& best() const { return vertices_[0]; }
                const Vertex& worst() const { return vertices_.back(); }
                const Vertex& secondWorst() const { return vertices_[vertices_.size()-2]; }

            };

            inline Simplex initialSimplex(const std::vector<double>& x0, const std::function<double(const std::vector<double>&)>& f, double perturbationFactor, SimplexInitializationMethod initMethod) {

                std::vector<Vertex> vertices;
                double a = perturbationFactor/(2*std::sqrt(double(x0.size()))); 
                vertices.push_back({x0,f(x0)});
                for (size_t i = 0; i < x0.size(); ++i) {
                    std::vector<double> x = x0; 
                    switch (initMethod)
                    {
                    case SimplexInitializationMethod::BASIC: x[i] += perturbationFactor; break;
                    case SimplexInitializationMethod::SCALED: x[i] += perturbationFactor * (1 + std::abs(x0[i])); break;
                    case SimplexInitializationMethod::SYMMETRIC: 
                        for (size_t j = 0; j < x0.size(); ++j) {
                            x[j] += (i == j) ? a : (-a / (x0.size() - 1));
                        };
                        break;
                    default: break;
                    }
                    vertices.push_back({x,f(x)});
                };
                return {vertices};
            }

            inline Vertex centroidVertex(const Simplex& simplex_, const std::function<double(const std::vector<double>&)>& f) {

                size_t n = simplex_.vertices_[0].points_.size();
                std::vector<Vertex> simplexVector = simplex_.vertices_;
                std::vector<double> x_c(n, 0); 
                double inv_n = 1.0 / double(n);
                for (size_t j = 0; j < n; ++j){
                    Vertex vj = simplexVector[j];
                    for (size_t i = 0; i < n; ++i){
                        x_c[i] +=  vj.points_[i]*inv_n;
                    }
                }
                return {x_c,f(x_c)};
            }

            inline Vertex reflectionVertex(const Vertex& centroid, const Vertex& worst, double reflectionFactor, const std::function<double(const std::vector<double>&)>& f) {

                size_t n = centroid.points_.size();
                std::vector<double> x(n); 
                for (size_t j = 0; j < n; ++j){
                    double xw = worst.points_[j];
                    double xc = centroid.points_[j]; 
                    x[j] = xc + reflectionFactor*(xc-xw);
                }
                return {x,f(x)};
            }

            inline Vertex expansionVertex(const Vertex& centroid, const Vertex& reflection, double expansionFactor, const std::function<double(const std::vector<double>&)>& f) {

                size_t n = centroid.points_.size();
                std::vector<double> x(n); 
                for (size_t j = 0; j < n; ++j){
                    double xr = reflection.points_[j];
                    double xc = centroid.points_[j]; 
                    x[j] = xc + expansionFactor*(xr-xc);
                }
                return {x,f(x)};
            }

            inline Vertex contractionVertex(const Vertex& centroid, const Vertex& worst, double contractionFactor, const std::function<double(const std::vector<double>&)>& f) {

                size_t n = centroid.points_.size();
                std::vector<double> x(n,0); 
                for (size_t j = 0; j < n; ++j){
                    double xw = worst.points_[j];
                    double xc = centroid.points_[j]; 
                    x[j] = xc + contractionFactor*(xw-xc);
                }
                return {x,f(x)};
            }

            inline Simplex shrinkSimplex(const Simplex& simplex_, double shrinkFactor, const std::function<double(const std::vector<double>&)>& f) {

                size_t n = simplex_.vertices_[0].points_.size();
                std::vector<Vertex> simplexOut;
                simplexOut.reserve(n+1);
                std::vector<Vertex> simplexIn = simplex_.vertices_;
                Vertex bestVertex = simplex_.best();
                simplexOut.push_back(bestVertex); 
                for (size_t j = 1; j < n; j++){
                    std::vector<double> x(n); 
                    for (size_t i = 0; i < n; i++){
                        double value = bestVertex.points_[i];
                        x[i] = value + shrinkFactor*(simplexIn[j].points_[i] - value);
                    }
                    simplexOut.push_back({x,f(x)}); 
                } 
                Simplex output{simplexOut};
                return output;
            }

            inline bool checkConvergence(const Simplex& simplex_, double toleranceThreshold) {

                size_t n = simplex_.vertices_[0].points_.size();
                double maxDistance = 0.0,maxValue = 0.0;
                Vertex bestVertex = simplex_.best();
                for (size_t i = 1; i < n; i++) {
                    double distance = 0.0;
                    Vertex iVertex = simplex_.vertices_[i];
                    for (int j = 0; j < n; j++) {
                        distance += std::pow(iVertex.points_[j] - bestVertex.points_[j], 2);
                    }
                    maxDistance = std::max(maxDistance, std::sqrt(distance));
                    maxValue = std::max(maxValue, std::abs(iVertex.value_ - bestVertex.value_));
                }
                return (maxDistance < toleranceThreshold || maxValue < toleranceThreshold);
            }


        }

        struct NelderMeadParameters {
            double perturbationFactor = 0.05, reflectionFactor=1.0, expansionFactor=2.0,contractionFactor=.5,shrinkFactor=.5;
            neldermeadtools::SimplexInitializationMethod initSimplexMethod = neldermeadtools::SimplexInitializationMethod::BASIC;
        };

        struct NelderMeadResult {std::vector<double> x; double f; int iter;};

        inline NelderMeadResult nelderMead(const std::vector<double>& x0, const std::function<double(const std::vector<double>&)>& f, int maxIter = 100, double toleranceThreshold = 1e-12, const NelderMeadParameters& nmParams = {}) {

            neldermeadtools::Simplex simplex = neldermeadtools::initialSimplex(x0,f,nmParams.perturbationFactor,nmParams.initSimplexMethod);
            int nIter = 0;
            for (int iter = 1; iter <= maxIter; iter++) { 

                nIter = iter;
                if (neldermeadtools::checkConvergence(simplex,toleranceThreshold)){break;}
                neldermeadtools::Vertex centroid = neldermeadtools::centroidVertex(simplex,f);
                neldermeadtools::Vertex worst = simplex.worst();
                neldermeadtools::Vertex reflection = neldermeadtools::reflectionVertex(centroid, worst, nmParams.reflectionFactor, f);
                neldermeadtools::Vertex best = simplex.best();
                neldermeadtools::Vertex second_worst = simplex.secondWorst();

                if (reflection.value_<best.value_){
                    neldermeadtools::Vertex expansion = neldermeadtools::expansionVertex(centroid, reflection,nmParams.expansionFactor,f); 
                    if (expansion.value_<reflection.value_){
                        simplex.setVertex(expansion, simplex.vertices_.size()-1);
                    }else{
                        simplex.setVertex(reflection, simplex.vertices_.size()-1);
                    }
                }
                else if (reflection.value_>=best.value_ and reflection.value_<second_worst.value_){
                    simplex.setVertex(reflection, simplex.vertices_.size()-1);
                }
                else {
                    if (reflection.value_<worst.value_){
                        neldermeadtools::Vertex contraction = neldermeadtools::contractionVertex(centroid, reflection,nmParams.contractionFactor,f); 
                        if (contraction.value_<reflection.value_){
                            simplex.setVertex(contraction, simplex.vertices_.size()-1);
                        }else{
                            simplex = neldermeadtools::shrinkSimplex(simplex,nmParams.shrinkFactor,f);
                        }
                    }else{
                        neldermeadtools::Vertex contraction = neldermeadtools::contractionVertex(centroid, worst,nmParams.contractionFactor,f); 
                        if (contraction.value_<worst.value_){
                            simplex.setVertex(contraction, simplex.vertices_.size()-1);
                        }else{
                            simplex = neldermeadtools::shrinkSimplex(simplex,nmParams.shrinkFactor,f);
                        }
                    }
                }
                if (iter==maxIter){break;}

            }
            
            return {simplex.vertices_[0].points_,simplex.vertices_[0].value_,nIter};
        }


    }
}