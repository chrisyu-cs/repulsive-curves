#include "multigrid/multigrid_operator.h"
#include <iostream>

namespace LWS {

    MultigridOperator::MultigridOperator() {}

    Eigen::VectorXd MultigridOperator::prolong(Eigen::VectorXd v, MultigridMode mode) {
        Eigen::VectorXd out;

        if (mode == MultigridMode::MatrixOnly) {
            assert(v.rows() == lowerSize);
            out.setZero(upperSize);
        }
        else if (mode == MultigridMode::Barycenter) {
            assert(v.rows() == lowerSize + 1);
            out.setZero(upperSize + 1);
        }

        for (size_t i = 0; i < matrices.size(); i++) {
            int outputStart = matrices[i].fineOffset;
            int inputStart = matrices[i].coarseOffset;
            int outputRows = matrices[i].M.rows();
            int inputRows = matrices[i].M.cols();

            out.block(outputStart, 0, outputRows, 1) = matrices[i].M * v.block(inputStart, 0, inputRows, 1);
        }

        if (mode == MultigridMode::Barycenter) {
            out(upperSize) = v(lowerSize);
        }

        return out;
    }

    Eigen::VectorXd MultigridOperator::restrictWithTranspose(Eigen::VectorXd v, MultigridMode mode) {
        Eigen::VectorXd out;

        if (mode == MultigridMode::MatrixOnly) {
            assert(v.rows() == upperSize);
            out.setZero(lowerSize);
        }
        else if (mode == MultigridMode::Barycenter) {
            assert(v.rows() == upperSize + 1);
            out.setZero(lowerSize + 1);
        }

        for (size_t i = 0; i < matrices.size(); i++) {
            int outputStart = matrices[i].coarseOffset;
            int inputStart = matrices[i].fineOffset;
            int outputRows = matrices[i].M.cols();
            int inputRows = matrices[i].M.rows();

            //out.block(outputStart, 0, outputRows, 1) = (v.block(inputStart, 0, inputRows, 1).transpose() *  matrices[i].M).transpose();
            out.block(outputStart, 0, outputRows, 1) = matrices[i].M.transpose() * v.block(inputStart, 0, inputRows, 1);
        }

        if (mode == MultigridMode::Barycenter) {
            out(lowerSize) = 0;
        }

        return out;
    }

    Eigen::VectorXd ApplyPinv(Eigen::SparseMatrix<double> &J, Eigen::VectorXd &x) {
        Eigen::SparseMatrix<double> JTJ = J.transpose() * J;
        Eigen::SparseLU<Eigen::SparseMatrix<double>, Eigen::COLAMDOrdering<int>> JTJ_solver;
        JTJ_solver.analyzePattern(JTJ);
        JTJ_solver.factorize(JTJ);

        Eigen::VectorXd JT_v = J.transpose() * x;
        Eigen::VectorXd pinv_v = JTJ_solver.solve(JT_v);
        return pinv_v;
    }

    Eigen::MatrixXd ApplyPinv(Eigen::SparseMatrix<double> &J, Eigen::MatrixXd &xs) {
        Eigen::MatrixXd out(J.cols(), xs.cols());
        for (int i = 0; i < xs.cols(); i++) {
            Eigen::VectorXd c_i = xs.col(i);
            out.col(i) = ApplyPinv(J, c_i);
        }
        return out;
    }

    Eigen::VectorXd MultigridOperator::restrictWithPinv(Eigen::VectorXd v, MultigridMode mode) {
        Eigen::VectorXd out;
        
        if (mode == MultigridMode::MatrixOnly) {
            assert(v.rows() == upperSize);
            out.setZero(lowerSize);
        }
        else if (mode == MultigridMode::Barycenter) {
            assert(v.rows() == upperSize + 1);
            out.setZero(lowerSize + 1);
        }

        for (size_t i = 0; i < matrices.size(); i++) {
            int outputStart = matrices[i].coarseOffset;
            int inputStart = matrices[i].fineOffset;
            int outputRows = matrices[i].M.cols();
            int inputRows = matrices[i].M.rows();

            Eigen::VectorXd vblock = v.block(inputStart, 0, inputRows, 1);

            //out.block(outputStart, 0, outputRows, 1) = (v.block(inputStart, 0, inputRows, 1).transpose() *  matrices[i].M).transpose();
            out.block(outputStart, 0, outputRows, 1) = ApplyPinv(matrices[i].M, vblock);
        }

        if (mode == MultigridMode::Barycenter) {
            out(lowerSize) = 0;
        }

        return out;
    }
}
