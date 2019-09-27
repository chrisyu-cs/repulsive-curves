#include "multigrid/multigrid_operator.h"
#include <iostream>

namespace LWS {

    MultigridOperator::MultigridOperator() {}

    Eigen::VectorXd MultigridOperator::prolong(Eigen::VectorXd v, MultigridMode mode) {
        Eigen::VectorXd out;

        if (mode == MultigridMode::MatrixOnly && lowerP) {
            std::cout << "[prolong] Barycenter before 1st projection = " << lowerP->EvaluateConstraints(v) << std::endl;
            v = lowerP->Multiply(v);
            std::cout << "[prolong] Barycenter before prolongation = " << lowerP->EvaluateConstraints(v) << std::endl;
        }

        if (mode == MultigridMode::MatrixOnly) {
            if (v.rows() != lowerSize) {
                std::cerr << "Size mismatch in prolong" << std::endl;
            }
            out.setZero(upperSize);
        }
        else if (mode == MultigridMode::Barycenter) {
            if (v.rows() != lowerSize + 1) {
                std::cerr << "Size mismatch in prolong" << std::endl;
            }
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

        if (mode == MultigridMode::MatrixOnly && upperP) {
            std::cout << "[prolong] Barycenter after prolongation = " << upperP->EvaluateConstraints(out) << std::endl;
            out = upperP->Multiply(out);
            std::cout << "[prolong] Barycenter after 2nd projection = " << upperP->EvaluateConstraints(out) << "\n" << std::endl;
        }

        return out;
    }

    Eigen::VectorXd MultigridOperator::restrictWithTranspose(Eigen::VectorXd v, MultigridMode mode) {
        Eigen::VectorXd out;

        if (mode == MultigridMode::MatrixOnly && upperP) {
            std::cout << "[restrictWithTranspose] Barycenter before 1st projection = " << upperP->EvaluateConstraints(v) << std::endl;
            v = upperP->Multiply(v);
            std::cout << "[restrictWithTranspose] Barycenter before restriction = " << upperP->EvaluateConstraints(v) << std::endl;
        }

        if (mode == MultigridMode::MatrixOnly) {
            if (v.rows() != upperSize) {
                std::cerr << "Size mismatch in restrictWithTranspose" << std::endl;
            }
            out.setZero(lowerSize);
        }
        else if (mode == MultigridMode::Barycenter) {
            if (v.rows() != upperSize + 1) {
                std::cerr << "Size mismatch in restrictWithTranspose" << std::endl;
            }
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
            out(lowerSize) = v(upperSize);
        }

        if (mode == MultigridMode::MatrixOnly && lowerP) {
            std::cout << "[restrictWithTranspose] Barycenter after restriction = " << lowerP->EvaluateConstraints(out) << std::endl;
            out = lowerP->Multiply(out);
            std::cout << "[restrictWithTranspose] Barycenter after 2nd projection = " << lowerP->EvaluateConstraints(out) << "\n" << std::endl;
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

        if (mode == MultigridMode::MatrixOnly && upperP) {
            v = upperP->Multiply(v);
        }
        
        if (mode == MultigridMode::MatrixOnly) {
            if (v.rows() != upperSize) {
                std::cerr << "Size mismatch in restrictWithPinv" << std::endl;
            }
            out.setZero(lowerSize);
        }
        else if (mode == MultigridMode::Barycenter) {
            if (v.rows() != upperSize + 1) {
                std::cerr << "Size mismatch in restrictWithPinv" << std::endl;
            }
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
            out(lowerSize) = v(upperSize);
        }

        if (mode == MultigridMode::MatrixOnly && lowerP) {
            out = lowerP->Multiply(out);
        }

        return out;
    }
}
