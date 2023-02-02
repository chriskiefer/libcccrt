#pragma once

#include <iostream>
#include <Eigen/Dense>


namespace fractal {
  struct sevcik {
    /* references:
    https://arxiv.org/pdf/1003.5266.pdf
    https://neuropsychology.github.io/NeuroKit/_modules/neurokit2/complexity/fractal_sevcik.html#fractal_sevcik
    */
    static double calc(const Eigen::VectorXd &sequence) {
      Eigen::ArrayXd y(sequence.size()-1);
      
      //scale 0 - 1
      Eigen::VectorXd sequence_scaled = sequence.array() - sequence.minCoeff();
      sequence_scaled = sequence_scaled / sequence_scaled.maxCoeff();

      //compute y
      for(size_t i=1; i < sequence_scaled.size(); i++) {
        y(i-1) = sequence_scaled(i) - sequence_scaled(i-1);
      }

      double x = 1.0 / (sequence_scaled.size()-1);
      x = x * x;

      y = (y * y);
      y = y + x;
      y =y.sqrt();
      double L = y.sum();
      double sfd = 1 + (log(L) / log(2 * y.size()));
      return sfd;
    }
  };
};