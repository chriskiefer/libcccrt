#pragma once

#include <iostream>
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <EigenRand>


struct RPC {

  static Eigen::MatrixXd createProjectionMatrix(size_t windowSize, size_t nDimensions) {
    Eigen::Rand::Vmt19937_64 urng{ 42 };
    Eigen::MatrixXd projectionMatrix = Eigen::Rand::normal<Eigen::MatrixXd>(nDimensions, windowSize, urng, 0, sqrt(1.0/nDimensions));
    return projectionMatrix;
  }

  static Eigen::MatrixXd calcProjection(const Eigen::MatrixXd &projectionMatrix, const Eigen::VectorXd &data, double hop=0.5) {
    size_t hopSize = std::max((size_t)1, static_cast<size_t>(projectionMatrix.cols() * hop));
    size_t nHops = static_cast<size_t>((data.size() - projectionMatrix.cols()) / hopSize) + 1;
    // cout << "Hopsize: " << hopSize << ", Hops: " << nHops << endl;
    // cout << projectionMatrix << endl;
    Eigen::MatrixXd projections(projectionMatrix.rows(), nHops);
    // cout << projections << endl;
    for(size_t i=0; i < nHops; i++) {
      // cout << "dw: " << i;
      // auto dataWindow = data(Eigen::seqN(i*hopSize, projectionMatrix.cols()));
      const auto dataWindow = data.segment(i*hopSize,projectionMatrix.cols());

      // cout << "\n" << dataWindow << endl;
      const auto projection = projectionMatrix * dataWindow;
      // cout << projection << endl;
      projections.col(i) = projection;
    }
    // cout << projections << endl; 
    return projections;
  }

  static inline size_t calcXdFlatArrayIndex(const Eigen::VectorXd &indexTuple, const size_t bound) {
    size_t index = indexTuple(0);
    for (size_t i =1; i < indexTuple.size(); i++) {
        index *= bound;
        index += indexTuple(i);
    }
    return index;
  }

  static double calculateProjectionArea(Eigen::MatrixXd &projections, const size_t resolution) {
    // cout << "area calc" << endl;
    size_t dims = projections.rows();
    // cout << projections << endl;
    //translate array to histogram bin indexes
    for(size_t row=0; row < projections.rows(); row++) {
      projections.row(row) = projections.row(row).array() - projections.row(row).minCoeff();
      // cout << "p:\n" << projections << endl;
      if (projections.row(row).maxCoeff() > 0) {
        projections.row(row) /=  (projections.row(row).maxCoeff() * 1.000001);
      }
      // cout << "p:\n" << projections << endl;
    }
    projections =  projections * resolution;
    // cout << "p:\n" << projections << endl;
    projections =  projections.array().floor();
    // cout << "p:\n" << projections << endl;

    //flat storage of multidimensional histogram, fill with zeros
    size_t histoSize = pow(resolution,dims);


    // Eigen::SparseVector<size_t> histo(histoSize);

    // // size_t area=0;
    // for(size_t i=0; i < projections.cols(); i++) {
    //   Eigen::VectorXd indexTuple = projections.col(i);
    //   size_t index = RPC::calcXdFlatArrayIndex(indexTuple, resolution);
    //   // cout << index << endl;
    //   // histo(index) = true;
    //   // if (histo.coeffRef(index)==0)
    //   //   area++;
    //   histo.coeffRef(index) = 1;
    //   // histo[index] = 1;
    // }


    // double area = histo.sum();
    // // cout << "Area: " << area << ", sum: " << histo.sum() << endl;
    // // return static_cast<double>(histo.size());
    // return area;
    
    std::vector<size_t> indexes(projections.cols());
    for(size_t i=0; i < projections.cols(); i++) {
      Eigen::VectorXd indexTuple = projections.col(i);
      size_t index = RPC::calcXdFlatArrayIndex(indexTuple, resolution);
      indexes[i] = index;
    }
    std::sort(indexes.begin(), indexes.end());
    auto last = std::unique(indexes.begin(), indexes.end());
    indexes.erase(last, indexes.end());
    return static_cast<double>(indexes.size());


  }

  static double calc(const Eigen::MatrixXd &projectionMatrix, const Eigen::VectorXd &data, const size_t resolution, double hop=0.5) {
    double res = 0;
    auto projections = RPC::calcProjection(projectionMatrix, data, hop);
    res = calculateProjectionArea(projections, resolution);
    return res;
  }


};