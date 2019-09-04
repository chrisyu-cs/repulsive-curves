#include "product/dense_matrix.h"

namespace LWS {

    DenseMatrixMult::DenseMatrixMult(Eigen::MatrixXd mat) {
        A = mat;
    }
}
