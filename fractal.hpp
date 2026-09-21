#pragma once

#include <iostream>
#include <Eigen/Dense>
#include "core/sevcik.hpp"


namespace fractal {
  // Eigen-facing wrapper around core/sevcik.hpp
  struct sevcik {
    /* references:
    https://arxiv.org/pdf/1003.5266.pdf
    https://neuropsychology.github.io/NeuroKit/_modules/neurokit2/complexity/fractal_sevcik.html#fractal_sevcik
    */
    static double calc(const Eigen::VectorXd &sequence) {
      return cccrt::sevcik(sequence.data(), sequence.size());
    }
  };
};
