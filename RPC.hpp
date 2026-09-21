#pragma once

#include <iostream>
#include <vector>
#include <Eigen/Dense>
#include "core/rpc.hpp"


// Eigen-facing wrapper around core/rpc.hpp
struct RPC {

  static Eigen::MatrixXd createProjectionMatrix(size_t windowSize, size_t nDimensions) {
    Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> projectionMatrix(nDimensions, windowSize);
    cccrt::rpc::makeProjectionMatrix(projectionMatrix.data(), nDimensions, windowSize, 42);
    return projectionMatrix;
  }

  static inline size_t calcXdFlatArrayIndex(const Eigen::VectorXd &indexTuple, const size_t bound) {
    return cccrt::rpc::flatIndex(indexTuple.data(), indexTuple.size(), bound);
  }

  static double calc(const Eigen::MatrixXd &projectionMatrix, const Eigen::VectorXd &data, const size_t resolution, double hop=0.5) {
    const size_t nDim = projectionMatrix.rows();
    const size_t windowSize = projectionMatrix.cols();
    // core expects a row-major matrix; MatrixXd is column-major
    Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> rowMajor = projectionMatrix;
    const size_t nHops = cccrt::rpc::numHops(data.size(), windowSize, cccrt::rpc::hopSizeInSamples(windowSize, hop));
    std::vector<double> proj(nDim * nHops);
    std::vector<uint64_t> cells(nHops);
    return cccrt::rpc::calc(rowMajor.data(), nDim, windowSize, data.data(), data.size(), resolution, hop, proj.data(), cells.data());
  }


};
