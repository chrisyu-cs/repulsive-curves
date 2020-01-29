#include "multigrid/multigrid_operator.h"
#include <iostream>

namespace LWS {

    MultigridOperator::MultigridOperator() {}

    void MultigridOperator::checkAndSetOutputSize(Eigen::VectorXd &in, Eigen::VectorXd &out, int inputVerts, int outputVerts, ProlongationMode mode) {

        int expectedInput = 0;
        int expectedOutput = 0;

        if (mode == ProlongationMode::Matrix3AndEdgeConstraints) {
            expectedInput = 4 * inputVerts + 3;
            expectedOutput = 4 * outputVerts + 3;
        }
        else if (mode == ProlongationMode::Matrix3AndProjector) {
            expectedInput = 3 * inputVerts;
            expectedOutput = 3 * outputVerts;
        }
        else {
            std::cerr << "Mode not implemented" << std::endl;
            throw 1;
        }

        if (in.rows() != expectedInput) {
            std::cerr << "Input doesn't match expected size: " << in.rows() << " vs " << expectedInput << std::endl;
            throw 1;
        }

        out.setZero(expectedOutput);

    }

    Eigen::VectorXd MultigridOperator::prolong(Eigen::VectorXd v, ProlongationMode mode) {
        Eigen::VectorXd out;

        checkAndSetOutputSize(v, out, lowerSize, upperSize, mode);

        if (mode == ProlongationMode::Matrix3AndEdgeConstraints) {
            prolongVerts3X(v, out);
            prolongBarycenter3X(v, out);
            prolongEdgeConstraints(v, out);
        }
        else if (mode == ProlongationMode::Matrix3AndProjector) {
            v = lowerP->ProjectToNullspace(v);
            prolongVerts3X(v, out);
            out = upperP->ProjectToNullspace(out);
        }
        else {
            prolongVerts1X(v, out); 
        }

        return out;
    }

    Eigen::VectorXd MultigridOperator::restrictWithTranspose(Eigen::VectorXd v, ProlongationMode mode) {
        Eigen::VectorXd out;

        checkAndSetOutputSize(v, out, upperSize, lowerSize, mode);

        if (mode == ProlongationMode::Matrix3AndEdgeConstraints) {
            restrictVertsTranspose3X(v, out);
            restrictBarycenter3X(v, out);
            restrictEdgeConstraintsTranspose(v, out);
        }
        else if (mode == ProlongationMode::Matrix3AndProjector) {
            v = upperP->ProjectToNullspace(v);
            restrictVertsTranspose3X(v, out);
            out = lowerP->ProjectToNullspace(out);
        }
        else {
            restrictVertsTranspose1X(v, out);
        }

        return out;
    }

    Eigen::VectorXd MultigridOperator::restrictWithPinv(Eigen::VectorXd v, ProlongationMode mode) {
        Eigen::VectorXd out;
        
        checkAndSetOutputSize(v, out, upperSize, lowerSize, mode);
        
        if (mode == ProlongationMode::Matrix3AndEdgeConstraints) {
            restrictVertsPinv3X(v, out);
            restrictBarycenter3X(v, out);
            restrictEdgeConstraintsPinv(v, out);
        }
        else if (mode == ProlongationMode::Matrix3AndProjector) {
            v = upperP->ProjectToNullspace(v);
            restrictVertsPinv3X(v, out);
            out = lowerP->ProjectToNullspace(out);
        }
        else {
            restrictVertsPinv1X(v, out);
        }

        return out;
    }
}
