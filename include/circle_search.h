#pragma once

#include <Eigen/Core>
#include "spatial/tpe_bvh.h"
#include "product/block_cluster_tree.h"

namespace LWS {

    class CircleSearch {
        // Compute the right-hand side of the second derivative system,
        // which is comprised of directional derivative of energy gradient,
        // and (directional derivative of Gram matrix) * energy gradient.
        public:
        static Eigen::VectorXd DifferentiateRHS(PolyCurveNetwork* curves, Eigen::VectorXd &projectedGradient,
        Eigen::MatrixXd &l2Gradient, BVHNode3D* bvh, BlockClusterTree* mult, Eigen::VectorXd &lambda,
        NullSpaceProjector* projector, double alpha, double beta, double epsilon);

        template<typename Solver, typename Smoother>
        static Eigen::VectorXd SolveSecondDerivative(PolyCurveNetwork* curves, Eigen::MatrixXd &projectedGradient,
        Eigen::MatrixXd &l2Gradient, BVHNode3D* bvh, Solver* solver, double alpha, double beta, double epsilon);

        template<typename Solver, typename Smoother>
        static double CircleSearchStep(PolyCurveNetwork* curves, Eigen::MatrixXd &projectedGradient,
        Eigen::MatrixXd &l2Gradient, BVHNode3D* bvh, Solver* solver, double gradDot, double alpha, double beta, double epsilon);

    };

    template<typename Solver, typename Smoother>
    Eigen::VectorXd CircleSearch::SolveSecondDerivative(PolyCurveNetwork* curves, Eigen::MatrixXd &projectedGradient,
    Eigen::MatrixXd &l2Gradient, BVHNode3D* bvh, Solver* solver, double alpha, double beta, double epsilon) {
        int nVerts = curves->NumVertices();
        BlockClusterTree* mult = solver->levels[0]->GetMultiplier();

        // Recover the value of lambda
        mult->SetBlockTreeMode(BlockTreeMode::Matrix3Only);
        Eigen::VectorXd lambda(nVerts * 3);
        lambda.setZero();
        Eigen::VectorXd projGrad3x(nVerts * 3);
        Eigen::VectorXd l2grad3x(nVerts * 3);
        projGrad3x.setZero();
        l2grad3x.setZero();
        MatrixIntoVectorX3(projectedGradient, projGrad3x);
        MatrixIntoVectorX3(l2Gradient, l2grad3x);
        // lambda now contains G * u
        mult->Multiply(projGrad3x, lambda);
        NullSpaceProjector* projector = solver->levels[0]->GetConstraintProjector();
        // l2grad3x now contains B * phi
        projector->ApplyB(l2grad3x, l2grad3x);
        // lambda now contains B * G * u
        projector->ApplyB(lambda, lambda);
        // lambda now contains B * phi - B * G * u
        lambda = l2grad3x - lambda;
        // Solve (BB^T) lambda = B phi - B G u to get lambda
        projector->SolveBBT(lambda, lambda);

        Eigen::VectorXd rhs = DifferentiateRHS(curves, projGrad3x, l2Gradient, bvh, mult, lambda,
            projector, alpha, beta, epsilon);

        int nConstraints = curves->NumEdges() + 3;
        Eigen::VectorXd w = rhs.block(0, 0, nVerts * 3, 1);
        Eigen::VectorXd mu = rhs.block(nVerts * 3, 0, nConstraints, 1);

        Eigen::VectorXd B_dagger(mu.rows());
        projector->ApplyBPinv(mu, B_dagger);
        Eigen::VectorXd GB_dagger(mu.rows());
        GB_dagger.setZero();
        mult->Multiply(B_dagger, GB_dagger);

        w = -w + GB_dagger;

        mult->SetBlockTreeMode(BlockTreeMode::Matrix3AndProjector);

        w = curves->constraintProjector->ProjectToNullspace(w);
        Eigen::VectorXd v = solver->template VCycleSolve<Smoother>(w);
        v = B_dagger - v;

        return v;
    }

    inline double alphaOfStep(double t, double R, double alpha_0, double alpha_1) {
        return (1.0 / R) * (alpha_0 * t + alpha_1 * t * t);
    }

    inline void setCircleStep(PolyCurveNetwork* curves, Eigen::MatrixXd &origPos, Eigen::MatrixXd &p_dot, Eigen::MatrixXd &K,
    double t, double R, double alpha_0, double alpha_1) {
        double alpha_of_t = alphaOfStep(t, R, alpha_0, alpha_1);
        double p_dot_coeff = sin(alpha_of_t) / alpha_0;
        double K_coeff = (1.0 - cos(alpha_of_t));

        curves->positions = origPos + R * (p_dot_coeff * -p_dot + R * K_coeff * K);
    }

    template<typename Solver, typename Smoother>
    double CircleSearch::CircleSearchStep(PolyCurveNetwork* curves, Eigen::MatrixXd &projectedGradient,
    Eigen::MatrixXd &l2Gradient, BVHNode3D* bvh, Solver* solver, double gradDot, double alpha, double beta, double epsilon) {
        int nVerts = curves->NumVertices();
        Eigen::VectorXd p_dotdot = SolveSecondDerivative<Solver, Smoother>(curves, projectedGradient, l2Gradient,
            bvh, solver, alpha, beta, epsilon);

        Eigen::VectorXd p_dot(nVerts * 3);
        p_dot.setZero();
        MatrixIntoVectorX3(projectedGradient, p_dot);
        
        BlockClusterTree* mult = solver->levels[0]->GetMultiplier();
        mult->SetBlockTreeMode(BlockTreeMode::Matrix3Only);
        double G_pdot_pdot = mult->DotProduct(p_dot, p_dot);
        double G_pdot_pdotdot = mult->DotProduct(p_dot, p_dotdot);

        Eigen::VectorXd K = (1.0 / G_pdot_pdot) * (p_dotdot - p_dot * (G_pdot_pdotdot / G_pdot_pdot));
        double R = 1.0 / sqrt(mult->DotProduct(K, K));
        double alpha_0 = sqrt(G_pdot_pdot);
        double alpha_1 = 0.5 * (G_pdot_pdotdot / alpha_0);
        double alpha_max = M_PI / 2;

        double gradNorm = projectedGradient.norm();
        double t_guess = (-alpha_0 + sqrt(alpha_0 * alpha_0 + 4 * R * alpha_1 * alpha_max)) / (2 * alpha_1);
        // t_guess = fmin(t_guess, 1.0 / gradNorm);

        Eigen::MatrixXd K_matrix(K.rows() / 3, 3);
        VectorXdIntoMatrix(K, K_matrix);
        Eigen::MatrixXd origPos = curves->positions;

        mult->SetBlockTreeMode(BlockTreeMode::Matrix3AndProjector);

        double initialEnergy = SpatialTree::TPEnergyBH(curves, bvh, alpha, beta);
        double newEnergy = initialEnergy;
        double sigma = 0.01;
        int numBacktracks = 0;
        double ls_step_threshold = 1e-10;

        // Line search
        while (fabs(t_guess) > ls_step_threshold) {
            // Search along this direction
            setCircleStep(curves, origPos, projectedGradient, K_matrix, t_guess, R, alpha_0, alpha_1);
            // Check the new energy
            bvh->recomputeCentersOfMass(curves);
            newEnergy = SpatialTree::TPEnergyBH(curves, bvh, alpha, beta);
            double decrease = initialEnergy - newEnergy;
            double targetDecrease = sigma * t_guess * gradNorm * gradDot;

            // If the energy hasn't decreased enough to meet the Armijo condition,
            // halve the step size.
            if (decrease < targetDecrease) {
                t_guess /= 2;
                numBacktracks++;
            }
            // Otherwise, accept the current step.
            else {
                // Update the centers of mass to reflect the new positions
                bvh->recomputeCentersOfMass(curves);
                break;
            }

            if (fabs(t_guess) <= ls_step_threshold) {
                std::cout << "Failed to find a non-trivial step after " << numBacktracks << " backtracks" << std::endl;
                // Restore initial positions if step size goes to 0
                curves->positions = origPos;
                return 0;
            }
        }

        // Backprojection

        std::cout << "t_guess = " << t_guess << " (" << numBacktracks << " backtracks)" << std::endl;
        std::cout << "Energy: " << initialEnergy << " -> " << newEnergy << std::endl;

        return fabs(t_guess);
    }
}
