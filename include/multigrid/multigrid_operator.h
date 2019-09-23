#pragma once

#include <vector>
#include <Eigen/Core>
#include <Eigen/Sparse>

namespace LWS {
    struct IndexedMatrix {
        Eigen::SparseMatrix<double> M;
        int fineOffset;
        int coarseOffset;
    };

    enum class MultigridMode {
        MatrixOnly,
        Barycenter,
        EdgeLengths
    };

    class MultigridOperator {
        public:
        MultigridOperator();
        int lowerSize;
        int upperSize;
        std::vector<IndexedMatrix> matrices;
        Eigen::VectorXd prolong(Eigen::VectorXd v, MultigridMode mode);
        Eigen::VectorXd restrictWithTranspose(Eigen::VectorXd v, MultigridMode mode);
        Eigen::VectorXd restrictWithPinv(Eigen::VectorXd v, MultigridMode mode);
    };

    Eigen::VectorXd ApplyPinv(Eigen::SparseMatrix<double> &J, Eigen::VectorXd &x);
    Eigen::MatrixXd ApplyPinv(Eigen::SparseMatrix<double> &J, Eigen::MatrixXd &xs);
}