#include <cstdint>
#include <vector>
#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;

// [[Rcpp::plugins("cpp11")]]

arma::vec sumX(const arma::cube& targetArr, std::int32_t vecLen);
arma::vec sumY(const arma::cube& targetArr, std::int32_t vecLen);
arma::vec sumZ(const arma::cube& targetArr, std::int32_t vecLen);

arma::cube repMatX(const arma::vec& sourceVec, std::int32_t len1, std::int32_t len2);
arma::cube repMatY(const arma::vec& sourceVec, std::int32_t len1, std::int32_t len2);
arma::cube repMatZ(const arma::vec& sourceVec, std::int32_t len1, std::int32_t len2);

// [[Rcpp::export]]
List hoge(
    arma::cube sourceCube,
    arma::vec xk,
    arma::vec yk,
    arma::vec zk,
    std::int32_t numIter
) {
    std::int32_t dimX = sourceCube.n_rows;
    std::int32_t dimY = sourceCube.n_cols;
    std::int32_t dimZ = sourceCube.n_slices;
    arma::vec recArrX = arma::vec(dimX);
    arma::vec recArrY = arma::vec(dimY);
    arma::vec recArrZ = arma::vec(dimZ);
    arma::cube reconstArray = sourceCube;

    for (std::int32_t i = 0; i < numIter; i++)
    {
        recArrX = sumX(reconstArray, dimX);
        reconstArray = reconstArray % 
            repMatX(xk / recArrX, dimY, dimZ);
        reconstArray.elem( find_nonfinite(reconstArray) ).zeros();
        
        recArrY = sumY(reconstArray, dimY);
        reconstArray = reconstArray %
            repMatY(yk / recArrY, dimX, dimZ);
        reconstArray.elem( find_nonfinite(reconstArray) ).zeros();
        
        recArrZ = sumZ(reconstArray, dimZ);
        reconstArray = reconstArray %
            repMatZ(zk / recArrZ, dimX, dimY);
        reconstArray.elem( find_nonfinite(reconstArray) ).zeros();
    }
    return List::create(reconstArray);
}

inline arma::vec sumX(const arma::cube& targetArr, std::int32_t vecLen) {
    // std::int32_t vecLen = targetArr.n_rows;
    arma::vec retVec = arma::vec(vecLen);
    for (std::int32_t i = 0; i < vecLen; i++)
    {
        retVec(i) = arma::accu(
            targetArr.subcube(
                arma::span(i),
                arma::span::all,
                arma::span::all
            )
        );
    }
    return retVec;
}

inline arma::vec sumY(const arma::cube& targetArr, std::int32_t vecLen) {
    // std::int32_t vecLen = targetArr.n_cols;
    arma::vec retVec = arma::vec(vecLen);
    for (std::int32_t i = 0; i < vecLen; i++)
    {
        retVec(i) = arma::accu(
            targetArr.subcube(
                arma::span::all,
                arma::span(i),
                arma::span::all
            )
        );
    }
    return retVec;
}

inline arma::vec sumZ(const arma::cube& targetArr, std::int32_t vecLen) {
    // std::int32_t vecLen = targetArr.n_slices;
    arma::vec retVec = arma::vec(vecLen);
    for (std::int32_t i = 0; i < vecLen; i++)
    {
        retVec(i) = arma::accu(
            targetArr.subcube(
                arma::span::all,
                arma::span::all,
                arma::span(i)
            )
        );
    }
    return retVec;
}

inline arma::cube repMatX(const arma::vec& sourceVec, std::int32_t len1, std::int32_t len2) {
    arma::cube retCube = arma::cube(sourceVec.n_rows, len1, len2);
    arma::mat sourceMat = repmat(sourceVec, 1, len1);
    retCube.each_slice() = sourceMat;
    return retCube;
}

inline arma::cube repMatY(const arma::vec& sourceVec, std::int32_t len1, std::int32_t len2) {
    arma::cube retCube = arma::cube(len1, sourceVec.n_rows, len2);
    arma::mat sourceMat = repmat(sourceVec, 1, len1);
    retCube.each_slice() = sourceMat.t();
    return retCube;
}

inline arma::cube repMatZ(const arma::vec& sourceVec, std::int32_t len1, std::int32_t len2) {
    std::int32_t zLen = sourceVec.n_rows;
    arma::cube retCube = arma::cube(len1, len2, zLen);
    // arma::mat sourceMat = arma::mat(len1, len2);
    for (std::int32_t i = 0; i < zLen; i++)
    {
        // sourceMat.fill(sourceVec(i));
        retCube.slice(i).fill(sourceVec(i));// = sourceMat;
    }
    return retCube;
}