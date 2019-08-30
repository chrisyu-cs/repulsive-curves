#pragma once

#include "vector_multiplier.h"
#include <vector>
#include "Eigen/Dense"

namespace LWS {

    class DenseMatrixMult : public VectorMultiplier<DenseMatrixMult> {
        public:
        DenseMatrixMult(Eigen::MatrixXd mat);
        void Multiply(std::vector<double> &v, std::vector<double> &b) const;

        private:
        Eigen::MatrixXd A;

    };
}
