#include "product/dense_matrix.h"

namespace LWS {

    DenseMatrixMult::DenseMatrixMult(Eigen::MatrixXd mat) {
        A = mat;
    }

    void DenseMatrixMult::Multiply(std::vector<double> &v, std::vector<double> &b) const {
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
