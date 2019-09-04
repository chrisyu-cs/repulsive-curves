#pragma once

#include "vector_multiplier.h"
#include <vector>
#include "Eigen/Dense"

namespace LWS {

    class DenseMatrixMult : public VectorMultiplier<DenseMatrixMult> {
        public:
        DenseMatrixMult(Eigen::MatrixXd mat);
        template<typename V, typename Dest>
        void Multiply(V &v, Dest &b) const;

        private:
        Eigen::MatrixXd A;

    };

    template<typename V, typename Dest>
    void DenseMatrixMult::Multiply(V &v, Dest &b) const {
        Eigen::VectorXd x;
        x.setZero(v.size());

        for (size_t i = 0; i < v.size(); i++) {
            x(i) = v[i];
        }

        Eigen::VectorXd Ax = A * x;

        for (size_t i = 0; i < v.size(); i++) {
            b[i] = Ax(i);
        }
    }
}
