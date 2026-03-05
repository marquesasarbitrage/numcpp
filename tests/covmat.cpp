#include "numcpp/objects/covmat.hpp"
#include <iostream>
#include "utils.hpp"
#include <cassert>

void testCorrelationMatrix3x3() {

    using namespace numcpp::objects; 

    Eigen::VectorXd corrCoeffs(3) ; 
    corrCoeffs << -0.006748268411, -0.2459694651, 0.03689723423; 

    Eigen::MatrixXd corrMatrix = Eigen::MatrixXd::Zero(3,3); 
    corrMatrix << 1.0, -0.006748268411,  -0.2459694651,
                -0.006748268411, 1.0, 0.03689723423,
                -0.2459694651, 0.03689723423, 1.0;
    
    Eigen::MatrixXd cholesky = Eigen::MatrixXd::Zero(3,3);  
    cholesky << 1.0, 0.0, 0.0,
                -0.006748268411, 1.0, 0.0,
                -0.2459694651, 0.03524, 0.9686;

    CorrelationMatrix corrMat(corrCoeffs); 

    assetIsMatrixesSimilar(corrMatrix,corrMat.matrix_);

    assert(corrMat.matrix_.rows()==3); 
    assert(corrMat.getCorrelation(0,0)==1); 
    assert(corrMat.getCorrelation(0,2)==-0.2459694651);

    assetIsMatrixesSimilar(cholesky,corrMat.getCholeskyDecomposition(), 1e-4);

}

void testCorrelationMatrix2x2() {

    using namespace numcpp::objects; 

    double corr = -0.006748268411;

    Eigen::MatrixXd corrMatrix(2,2); 
    corrMatrix << 1.0, -0.006748268411,
                -0.006748268411, 1.0;
    
    Eigen::MatrixXd cholesky(2,2); 
    cholesky << 1.0, 0.0, 
                -0.007, 1.0;

    CorrelationMatrix corrMat(corr); 
    assetIsMatrixesSimilar(corrMatrix,corrMat.matrix_);

    assert(corrMat.matrix_.rows()==2); 
    assert(corrMat.getCorrelation(0,0)==1); 
    assert(corrMat.getCorrelation(0,1)==-0.006748268411);
  
    assetIsMatrixesSimilar(cholesky,corrMat.getCholeskyDecomposition(), 1e-3);

}

void testCorrelationMatrix3x3MatrixConstructor() {

    using namespace numcpp::objects; 
    Eigen::MatrixXd corrMatrix = Eigen::MatrixXd::Zero(3,3); 
    corrMatrix << 1.0, -0.006748268411, -0.2459694651, 
        -0.006748268411, 1.0, 0.03689723423, 
        -0.2459694651, 0.03689723423, 1.0;

    CorrelationMatrix corrMat(corrMatrix);

    Eigen::MatrixXd cholesky = Eigen::MatrixXd::Zero(3,3); 
    cholesky <<  1.0, 0.0, 0.0, 
        -0.006748268411, 1.0, 0.0, 
        -0.2459694651, 0.03524, 0.9686; 
    

    assetIsMatrixesSimilar(corrMatrix,corrMat.matrix_);

    assert(corrMat.getCorrelation(0,0)==1); 
    assert(corrMat.getCorrelation(0,2)==-0.2459694651);

    assetIsMatrixesSimilar(cholesky,corrMat.getCholeskyDecomposition(), 1e-4);


}

void testCorrelationMatrixNonSquareError() {

    using namespace numcpp::objects; 
    try {
        Eigen::MatrixXd corrMatrix(3,4); 
        corrMatrix << 1.0, -0.006748268411, -0.2459694651,1.0, 
            -0.006748268411, 1.0, 0.03689723423, 1.0, 
            -0.2459694651, 0.03689723423, 1.0, 1.0;

        CorrelationMatrix corrMat(corrMatrix);
        assert(false);

    } catch (const numcpp::errors::Error& e) {

        assert(true);
    } catch (const std::exception& e) {

        assert(false);
    }

    try {
        Eigen::MatrixXd corrMatrix(3,3); 
        corrMatrix << 1.0, -0.03689723423, -0.2459694651, 
            -0.006748268411, 1.0, 0.03689723423, 
            -0.2459694651, 0.03689723423, 1.0;

        CorrelationMatrix corrMat(corrMatrix);
        assert(false);

    } catch (const numcpp::errors::Error& e) {

        assert(true);
    } catch (const std::exception& e) {

        assert(false);
    }

    try {

        Eigen::MatrixXd corrMatrix(3,3); 
        corrMatrix << 1.0, -0.03689723423, -0.2459694651, 
            -0.006748268411, 1.0, 0.03689723423, 
            -0.2459694651, 0.03689723423, 1.0;

        CorrelationMatrix corrMat(corrMatrix);
        assert(false);

    } catch (const numcpp::errors::Error& e) {

        assert(true);
    } catch (const std::exception& e) {

        assert(false);
    }
}

void testCorrelationMatrixNonSymmetricError() {

    try {
        using namespace numcpp::objects; 
        Eigen::MatrixXd corrMatrix(3,3);
        corrMatrix << 
            1.0, -0.03689723423, -0.2459694651, 
            -0.006748268411, 1.0, 0.03689723423, 
            -0.2459694651, 0.03689723423, 1.0;

        CorrelationMatrix corrMat(corrMatrix);
        assert(false);

    } catch (const numcpp::errors::Error& e) {

        assert(true);
    } catch (const std::exception& e) {

        assert(false);
    }

}

void testCorrelationMatrixNonPositiveSemiDefiniteError() {

    try {
        using namespace numcpp::objects; 
        Eigen::MatrixXd corrMatrix(3,3);
        corrMatrix << 1.0, -1.0, 0.0, 
            -1.0, 1.0, 0.0, 
            0.0, 0.0, 1.0;

        CorrelationMatrix corrMat(corrMatrix);
        assert(false);

    } catch (const numcpp::errors::Error& e) {

        assert(true);
    } catch (const std::exception& e) {

        assert(false);
    }
}

void testCorrelationMatrixInvalidCoefficients() {

    using namespace numcpp::objects; 
    double corr = -2.0;
    try {
        
        CorrelationMatrix corrMat(-2.0);
        assert(false);  
    } catch (const numcpp::errors::Error& e) {

        assert(true);

    } catch (const std::exception& e ) {

        assert(false);
    }


    try {
        Eigen::VectorXd coefs(4);
        coefs << 0.6, 0.6, 0.6, 0.6;
        CorrelationMatrix corrMat(coefs);
        assert(false);  
    } catch (const numcpp::errors::Error& e) {

        assert(true);

    } catch (const std::exception& e ) {

        assert(false);
    }

    try {
        Eigen::VectorXd coefs(6);
        coefs << 0.6, 0.6, 0.6, 0.6, 0.6,0.6;
        //Eigen::VectorXd coefs = {0.6, 0.6, 0.6, 0.6, 0.6, 0.6};
        CorrelationMatrix corrMat(coefs);
        assert(true);  
    } catch (const numcpp::errors::Error& e) {

        assert(false);

    } catch (const std::exception& e ) {

        assert(false);
    }



    
}

void testCovarianceMatrixConstructor() {

    using namespace numcpp::objects; 
    Eigen::MatrixXd covMatrix(3,3);
    Eigen::MatrixXd corrMatrix(3,3);

    covMatrix << 0.0004026385246,-0.000005150756644,-0.0001141394503, 
                -0.000005150756644,0.00144691174,0.00003245725223, 
                -0.0001141394503,	0.00003245725223, 0.0005348029993;

    corrMatrix << 1.0, -0.006748268411, -0.245969, 
                -0.006748268411, 1.0, 0.03689723423, 
                -0.245969, 0.03689723423, 1.0;

    Eigen::VectorXd variances(3);  
    variances << 0.0004026385246,0.00144691174,0.0005348029993;

    CovarianceMatrix cov1(covMatrix);
    assert(cov1.correlationMatrix_.matrix_.rows()==3);
    assert(cov1.getCovariance(0,1) == -0.000005150756644);
    assert(cov1.getVariance(1) == 0.00144691174);
    assert(cov1.getStandardDeviation(1) == std::sqrt(0.00144691174));
    assetIsMatrixesSimilar(covMatrix,cov1.get(), 1e-3);
    assetIsMatrixesSimilar(corrMatrix,cov1.correlationMatrix_.matrix_, 1e-3);
    CovarianceMatrix cov2(CorrelationMatrix(corrMatrix),variances);

    assert(cov2.correlationMatrix_.matrix_.rows()==3);
    assert(isClose(cov2.getCovariance(0,1), -0.000005150756644, 1e-9));
    assert(cov2.getVariance(1) == 0.00144691174);
    assert(cov2.getStandardDeviation(1) == std::sqrt(0.00144691174));
    assetIsMatrixesSimilar(covMatrix,cov2.get(), 1e-3);
    assetIsMatrixesSimilar(corrMatrix,cov2.correlationMatrix_.matrix_, 1e-3);
}

void testCovarianceMatrixMismatchError() {

    using namespace numcpp::objects; 
    try {
        Eigen::MatrixXd corrMatrix(3,3);
        corrMatrix << 
            1.0, -0.006748268411, -0.2459694651, 
            -0.006748268411, 1.0, 0.03689723423, 
            -0.2459694651, 0.03689723423, 1.0;

        Eigen::VectorXd variances(4);  
        variances << 0.0004026385246,0.00144691174,0.0005348029993, 0.0005348029993;
        
        CovarianceMatrix cov1(corrMatrix,variances);
        assert(false);

    } catch (const numcpp::errors::Error& e) {

        assert(true);

    } catch (const std::exception& e) {

        assert(false);
    }
}

void testCovarianceMatrixNegativeVarianceError() {
    using namespace numcpp::objects; 
    try {
        Eigen::MatrixXd corrMatrix(3,3);
        corrMatrix << 
            1.0, -0.006748268411, -0.2459694651, 
            -0.006748268411, 1.0, 0.03689723423, 
            -0.2459694651, 0.03689723423, 1.0;

        Eigen::VectorXd variances(3);
        variances << 0.0004026385246,0.00144691174,-0.0005348029993;
        
        CovarianceMatrix cov1(corrMatrix,variances);
        assert(false);

    } catch (const numcpp::errors::Error& e) {

        assert(true);

    } catch (const std::exception& e) {

        assert(false);
    }
}

void testInvalidCorrelationMatrix() {

    using namespace numcpp::objects; 
    Eigen::MatrixXd corrMatrix(3,3);
    corrMatrix << 1.0, -0.006748268411, -0.245969, 
                -0.006748268411, 1.0, 0.03689723423, 
                -0.245969, 0.98, .03689723423;

    try {

        CorrelationMatrix cov1(corrMatrix);
        assert(false); 
    } catch (const numcpp::errors::Error& e) {

        assert(true);

    } catch (const std::exception& e) {

        assert(false);
    }
}

int main() {
    testCorrelationMatrix3x3(); 
    testCorrelationMatrix2x2(); 
    testCorrelationMatrix3x3MatrixConstructor(); 
    testCorrelationMatrixNonSquareError(); 
    testCorrelationMatrixNonSymmetricError(); 
    testCorrelationMatrixNonPositiveSemiDefiniteError(); 
    testCorrelationMatrixInvalidCoefficients(); 
    testCovarianceMatrixConstructor(); 
    testCovarianceMatrixMismatchError(); 
    testCovarianceMatrixNegativeVarianceError(); 
    testInvalidCorrelationMatrix();
    std::cout << "All tests for objects module have been passed successfully!" << std::endl;
    return 0; 
}