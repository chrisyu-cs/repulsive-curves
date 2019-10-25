#include "circle_search.h"
#include "flow/gradient_constraint_types.h"

namespace LWS {

    Eigen::VectorXd CircleSearch::DifferentiateRHS(PolyCurveNetwork* curves, Eigen::VectorXd &projectedGradient,
    Eigen::MatrixXd &l2Gradient, BVHNode3D* bvh, BlockClusterTree* mult, Eigen::VectorXd &lambda,
    NullSpaceProjector* projector, double alpha, double beta, double epsilon) {
        mult->SetBlockTreeMode(BlockTreeMode::Matrix3Only);
        // Space for the difference quotient for (derivative Gram) * gradient
        Eigen::VectorXd gramTimesGradient, gramTimesGradientEps;
        int nVerts = curves->NumVertices();
        int nEdges = curves->NumEdges();
        int nConstraints = lambda.rows();;
        gramTimesGradient.setZero(nVerts * 3 + nConstraints);
        gramTimesGradientEps.setZero(nVerts * 3 + nConstraints);

        // Multiply Gram * gradient now
        mult->Multiply(projectedGradient, gramTimesGradient);

        // Add constraint blocks using lambda
        gramTimesGradient.block(0, 0, nVerts * 3, 1) += projector->ApplyBTranspose(lambda);
        gramTimesGradient.block(nVerts * 3, 0, nConstraints, 1) += projector->ApplyB(projectedGradient);

        // We have the initial l2 gradient already, so now we just need
        // the epsilon-perturbed vectors
        Eigen::MatrixXd origPositions = curves->positions;
        double projNorm = projectedGradient.norm();
        // Move by epsilon in the search direction
        curves->positions -= (epsilon  / projNorm) * (projectedGradient);
        
        // Get the epsilon-perturbed L2 gradient
        Eigen::MatrixXd l2GradientEps;
        l2GradientEps.setZero(nVerts, 3);
        bvh->recomputeCentersOfMass(curves);
        SpatialTree::TPEGradientBarnesHut(curves, bvh, l2GradientEps, alpha, beta);

        // Multiply (Gram + eps) * gradient
        mult->refreshEdgeWeights();
        mult->Multiply(projectedGradient, gramTimesGradientEps);
        // Assemble new constraint block
        Eigen::SparseMatrix<double> B_eps;
        EdgeLengthConstraint constraint(curves);
        constraint.FillConstraintMatrix(B_eps);
        // Use new block to multiply with lambda
        gramTimesGradientEps.block(0, 0, nVerts * 3, 1) += B_eps.transpose() * lambda;
        gramTimesGradientEps.block(nVerts * 3, 0, nConstraints, 1) += B_eps * projectedGradient;

        // Numerical derivative of Gram * gradient
        gramTimesGradientEps = (gramTimesGradientEps - gramTimesGradient) / epsilon;

        // Numerical derivative of Barnes-Hut
        l2GradientEps = (l2GradientEps - l2Gradient) / epsilon;
        Eigen::VectorXd sum;
        sum.setZero(nVerts * 3 + nConstraints);
        MatrixIntoVectorX3(l2GradientEps, sum);

        // Reset original values
        curves->positions = origPositions;
        bvh->recomputeCentersOfMass(curves);
        mult->refreshEdgeWeights();
        mult->SetBlockTreeMode(BlockTreeMode::Matrix3AndProjector);
        
        sum = sum + gramTimesGradientEps;
        return sum;
    }
}