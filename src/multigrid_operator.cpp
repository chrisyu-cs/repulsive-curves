#include "multigrid_operator.h"
#include <iostream>

namespace LWS {

    MultigridOperator::MultigridOperator() {}


    Eigen::VectorXd MultigridOperator::mapUpward(Eigen::VectorXd v, MultigridMode mode) {
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

    Eigen::VectorXd MultigridOperator::mapDownward(Eigen::VectorXd v, MultigridMode mode) {
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
}
