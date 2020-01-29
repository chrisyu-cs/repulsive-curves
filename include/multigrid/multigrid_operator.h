#pragma once

#include <vector>
#include <Eigen/Core>
#include <Eigen/Sparse>
#include "nullspace_projector.h"

namespace LWS {
    struct IndexedMatrix {
        Eigen::SparseMatrix<double> M;
        int fineOffset;
        int coarseOffset;
    };
    
    template<typename V>
    V ApplyPinv(Eigen::SparseMatrix<double> &J, V &x) {
        Eigen::SparseMatrix<double> JTJ = J.transpose() * J;
        Eigen::SparseLU<Eigen::SparseMatrix<double>, Eigen::COLAMDOrdering<int>> JTJ_solver;
        JTJ_solver.analyzePattern(JTJ);
        JTJ_solver.factorize(JTJ);

        return JTJ_solver.solve(J.transpose() * x);
    }

    template<> inline
    Eigen::MatrixXd ApplyPinv(Eigen::SparseMatrix<double> &J, Eigen::MatrixXd &xs) {
        Eigen::MatrixXd out(J.cols(), xs.cols());
        for (int i = 0; i < xs.cols(); i++) {
            Eigen::VectorXd c_i = xs.col(i);
            out.col(i) = ApplyPinv(J, c_i);
        }
        return out;
    }

    class MultigridOperator {
        public:
        MultigridOperator() {}
        virtual ~MultigridOperator() {}

        int lowerSize;
        int upperSize;
        NullSpaceProjector* lowerP;
        NullSpaceProjector* upperP;
        
        std::vector<IndexedMatrix> matrices;
        std::vector<IndexedMatrix> edgeMatrices;
        
        virtual Eigen::VectorXd prolong(Eigen::VectorXd v) = 0;
        virtual Eigen::VectorXd restrictWithTranspose(Eigen::VectorXd v) = 0;
        virtual Eigen::VectorXd restrictWithPinv(Eigen::VectorXd v) = 0;
    };
}