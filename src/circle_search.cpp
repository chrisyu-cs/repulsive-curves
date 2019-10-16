#include "circle_search.h"

namespace LWS {

    Eigen::MatrixXd CircleSearch::DifferentiateRHS(PolyCurveNetwork* curves, Eigen::MatrixXd &projectedGradient,
    Eigen::MatrixXd &l2Gradient, BVHNode3D* bvh, BlockClusterTree* mult, double alpha, double beta, double epsilon) {
        mult->SetBlockTreeMode(BlockTreeMode::Matrix3Only);
        // Space for the difference quotient for (derivative Gram) * gradient
        Eigen::MatrixXd gramTimesGradient, gramTimesGradientEps;
        int nVerts = curves->NumVertices();
        gramTimesGradient.setZero(nVerts, 3);
        gramTimesGradientEps.setZero(nVerts, 3);

        // Multiply Gram * gradient now
        mult->Multiply(projectedGradient, gramTimesGradient);

        // We have the initial l2 gradient already, so now we just need
        // the epsilon-perturbed vectors
        Eigen::MatrixXd origPositions = curves->positions;
        double projNorm = projectedGradient.norm();
        // Move by epsilon in the normalized projected gradient direction
        curves->positions += (epsilon  / projNorm) * (projectedGradient);
        bvh->recomputeCentersOfMass(curves);
        mult->refreshEdgeWeights();
        
        // Get the epsilon-perturbed L2 gradient
        Eigen::MatrixXd l2GradientEps;
        l2GradientEps.setZero(nVerts, 3);
        SpatialTree::TPEGradientBarnesHut(curves, bvh, l2GradientEps, alpha, beta);

        // Multiply (Gram + eps) * gradient
        mult->Multiply(projectedGradient, gramTimesGradientEps);
        // Numerical derivative of Gram * gradient
        gramTimesGradientEps = (gramTimesGradientEps - gramTimesGradient) / epsilon;

        // Numerical derivative of Barnes-Hut
        l2GradientEps = (l2GradientEps - l2Gradient) / epsilon;
        curves->positions = origPositions;
        
        mult->SetBlockTreeMode(BlockTreeMode::Matrix3AndProjector);

        return l2GradientEps + gramTimesGradientEps;
    }
}