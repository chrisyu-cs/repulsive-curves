#include "multigrid/multigrid_operator.h"
#include <iostream>

namespace LWS {

    MultigridOperator::MultigridOperator() {}

    Eigen::VectorXd MultigridOperator::prolong(Eigen::VectorXd v, ProlongationMode mode) {
        Eigen::VectorXd out;

        if (mode == ProlongationMode::MatrixAndProjector) {
            // std::cout << "[prolong] Barycenter before 1st projection = " << lowerP->EvaluateConstraints(v) << std::endl;
            v = lowerP->ProjectToNullspace(v);
            // std::cout << "[prolong] Barycenter before prolongation = " << lowerP->EvaluateConstraints(v) << std::endl;
        }

        if (mode == ProlongationMode::MatrixOnly || mode == ProlongationMode::MatrixAndProjector) {
            if (v.rows() != lowerSize) {
                std::cerr << "Size mismatch in prolong" << std::endl;
            }
            out.setZero(upperSize);
        }
        else if (mode == ProlongationMode::Barycenter) {
            if (v.rows() != lowerSize + 1) {
                std::cerr << "Size mismatch in prolong" << std::endl;
            }
            out.setZero(upperSize + 1);
        }
        else if (mode == ProlongationMode::Matrix3AndEdgeConstraints) {
            out.setZero(4 * upperSize + 3);
        }

        if (mode == ProlongationMode::Matrix3AndEdgeConstraints) {
            prolongVerts3X(v, out);
            prolongBarycenter3X(v, out);
            prolongEdgeConstraints(v, out);
        }
        else {
            prolongVerts1X(v, out); 
        }

        if (mode == ProlongationMode::MatrixAndProjector) {
            out = upperP->ProjectToNullspace(out);
        }

        if (mode == ProlongationMode::Barycenter) {
            out(upperSize) = v(lowerSize);
        }

        return out;
    }

    Eigen::VectorXd MultigridOperator::restrictWithTranspose(Eigen::VectorXd v, ProlongationMode mode) {
        Eigen::VectorXd out;

        if (mode == ProlongationMode::MatrixAndProjector) {
            v = upperP->ProjectToNullspace(v);
        }

        if (mode == ProlongationMode::MatrixOnly || mode == ProlongationMode::MatrixAndProjector) {
            if (v.rows() != upperSize) {
                std::cerr << "Size mismatch in restrictWithTranspose" << std::endl;
            }
            out.setZero(lowerSize);
        }
        else if (mode == ProlongationMode::Barycenter) {
            if (v.rows() != upperSize + 1) {
                std::cerr << "Size mismatch in restrictWithTranspose" << std::endl;
            }
            out.setZero(lowerSize + 1);
        }
        else if (mode == ProlongationMode::Matrix3AndEdgeConstraints) {
            out.setZero(4 * lowerSize + 3);
            std::cout << "lowerSize = " << lowerSize << std::endl;
        }

        if (mode == ProlongationMode::Matrix3AndEdgeConstraints) {
            restrictVertsTranspose3X(v, out);
            restrictBarycenter3X(v, out);
            restrictEdgeConstraintsTranspose(v, out);
        }
        else {
            restrictVertsTranspose1X(v, out);
        }

        if (mode == ProlongationMode::MatrixAndProjector) {
            out = lowerP->ProjectToNullspace(out);
        }

        if (mode == ProlongationMode::Barycenter) {
            out(lowerSize) = v(upperSize);
        }

        return out;
    }

    Eigen::VectorXd MultigridOperator::restrictWithPinv(Eigen::VectorXd v, ProlongationMode mode) {
        Eigen::VectorXd out;

        if (mode == ProlongationMode::MatrixAndProjector) {
            v = upperP->ProjectToNullspace(v);
        }
        
        if (mode == ProlongationMode::MatrixOnly || mode == ProlongationMode::MatrixAndProjector) {
            if (v.rows() != upperSize) {
                std::cerr << "Size mismatch in restrictWithPinv" << std::endl;
            }
            out.setZero(lowerSize);
        }
        else if (mode == ProlongationMode::Barycenter) {
            if (v.rows() != upperSize + 1) {
                std::cerr << "Size mismatch in restrictWithPinv" << std::endl;
            }
            out.setZero(lowerSize + 1);
        }
        else if (mode == ProlongationMode::Matrix3AndEdgeConstraints) {
            out.setZero(4 * lowerSize + 3);
        }
        
        if (mode == ProlongationMode::Matrix3AndEdgeConstraints) {
            restrictVertsPinv3X(v, out);
            restrictBarycenter3X(v, out);
            restrictEdgeConstraintsPinv(v, out);
        }
        else {
            restrictVertsPinv1X(v, out);
        }

        if (mode == ProlongationMode::MatrixAndProjector) {
            out = lowerP->ProjectToNullspace(out);
        }

        if (mode == ProlongationMode::Barycenter) {
            out(lowerSize) = v(upperSize);
        }

        return out;
    }
}
