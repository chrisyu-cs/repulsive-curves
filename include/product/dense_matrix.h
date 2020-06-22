#pragma once

#include "libgmultigrid/matrix_free.h"
#include "libgmultigrid/vector_multiplier.h"
#include <vector>
#include "Eigen/Dense"

namespace LWS {

    class DenseMatrixMult : public VectorMultiplier<DenseMatrixMult> {
        public:
        DenseMatrixMult(Eigen::MatrixXd mat) {
            A = mat;
        }

        template<typename V, typename Dest>
        void Multiply(V &v, Dest &b) const;

        private:
        Eigen::MatrixXd A;

    };

    template<typename V, typename Dest>
    void DenseMatrixMult::Multiply(V &v, Dest &b) const {
        Eigen::VectorXd x;
        int nRows = v.size();
        x.setZero(nRows);

        for (int i = 0; i < nRows; i++) {
            x(i) = v[i];
        }

        Eigen::VectorXd Ax = A * x;

        for (int i = 0; i < nRows; i++) {
            b[i] = Ax(i);
        }
    }

    using WrappedMatrix = Product::MatrixReplacement<DenseMatrixMult>;
}
