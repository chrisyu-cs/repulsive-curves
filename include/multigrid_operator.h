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
        Eigen::VectorXd mapUpward(Eigen::VectorXd v, MultigridMode mode);
        Eigen::VectorXd mapDownward(Eigen::VectorXd v, MultigridMode mode);
    };
}