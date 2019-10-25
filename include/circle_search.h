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
        Eigen::MatrixXd &dE, BVHNode3D* bvh, BlockClusterTree* mult, Eigen::VectorXd &lambda,
        NullSpaceProjector* projector, double alpha, double beta, double epsilon);

        template<typename Solver, typename Smoother>
        static Eigen::VectorXd SolveSecondDerivative(PolyCurveNetwork* curves, Eigen::MatrixXd &projectedGradient,
        Eigen::MatrixXd &dE, BVHNode3D* bvh, Solver* solver, double alpha, double beta, double epsilon);

        template<typename Solver, typename Smoother>
        static double CircleSearchStep(PolyCurveNetwork* curves, Eigen::MatrixXd &projectedGradient,
        Eigen::MatrixXd &dE, BVHNode3D* bvh, Solver* solver, std::vector<double> &targetLengths,
        double gradDot, double alpha, double beta, double epsilon);

    };

    template<typename Solver, typename Smoother>
    Eigen::VectorXd CircleSearch::SolveSecondDerivative(PolyCurveNetwork* curves, Eigen::MatrixXd &projectedGradient,
    Eigen::MatrixXd &dE, BVHNode3D* bvh, Solver* solver, double alpha, double beta, double epsilon) {
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
        MatrixIntoVectorX3(dE, l2grad3x);
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

        Eigen::VectorXd rhs = DifferentiateRHS(curves, projGrad3x, dE, bvh, mult, lambda,
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
    double t, double R, double alpha_0, double alpha_1, Eigen::MatrixXd &p_dotdot) {
        double alpha_of_t = alphaOfStep(t, R, alpha_0, alpha_1);
        // alpha_of_t = M_PI;
        double p_dot_coeff = sin(alpha_of_t) / alpha_0;
        double K_coeff = (1.0 - cos(alpha_of_t));

        curves->positions = origPos + R * (p_dot_coeff * p_dot + R * K_coeff * K);
        // curves->positions = origPos + -p_dot * t + 0.5 * t * t * -p_dotdot;
        // curves->positions = origPos + -p_dot * t;
    }

    inline void plotSteps(PolyCurveNetwork* curves, Eigen::MatrixXd &origPos, Eigen::MatrixXd &p_dot, Eigen::MatrixXd &K,
    double R, double alpha_0, double alpha_1, BVHNode3D* bvh, std::vector<double> &targetLengths, double gradDot, Eigen::MatrixXd &dE) {
        for (double t = 0; t < 2e-2; t += 1e-4) {
            // Eigen::MatrixXd realDE(curves->NumVertices(), 3);
            // realDE.setZero();
            // TPESC::FillGradientVectorDirect(curves, realDE, 2, 4);

            setCircleStep(curves, origPos, p_dot, K, t, R, alpha_0, alpha_1, p_dot);
            bvh->recomputeCentersOfMass(curves);

            Eigen::VectorXd constraints(curves->NumEdges() + 3);
            constraints.setZero();
            curves->FillConstraintViolations(constraints, targetLengths);
            double constraintNorm = constraints.norm();

            double dotLengths = 0;

            for (int i = 0; i < curves->NumVertices(); i++) {
                for (int j = 0; j < 3; j++) {
                    dotLengths += p_dot(i, j) * dE(i, j);
                }
            }

            double targetDecrease = dotLengths * t;

            double energy = SpatialTree::TPEnergyBH(curves, bvh, 2, 4);
            double realEnergy = 0; //TPESC::tpe_total(curves, 2, 4);

            std::cout << t << ", " << energy << ", " << -targetDecrease << ", " << realEnergy << std::endl;
        }
    }

    template<typename Solver, typename Smoother>
    double CircleSearch::CircleSearchStep(PolyCurveNetwork* curves, Eigen::MatrixXd &projectedGradient,
    Eigen::MatrixXd &dE, BVHNode3D* bvh, Solver* solver, std::vector<double> &targetLengths,
    double gradDot, double alpha, double beta, double epsilon) {
        int nVerts = curves->NumVertices();
        Eigen::VectorXd p_dotdot = SolveSecondDerivative<Solver, Smoother>(curves, projectedGradient, dE,
            bvh, solver, alpha, beta, epsilon);

        Eigen::MatrixXd p_dotdot_mat;
        p_dotdot_mat.setZero(nVerts, 3);
        VectorXdIntoMatrix(p_dotdot, p_dotdot_mat);

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
        // double t_guess = (gradNorm > 1) ? 1.0 / gradNorm : 1.0 / sqrt(gradNorm);
        double t_guess = -2 * alpha_0 / alpha_1;

        if (alphaOfStep(t_guess, R, alpha_0, alpha_1) > alpha_max) {
            double alpha_t = alphaOfStep(t_guess, R, alpha_0, alpha_1);
            std::cout << "Initial guess " << t_guess << " had alpha " << alpha_t << " greater than max " << alpha_max << std::endl;
            t_guess = (-alpha_0 + sqrt(alpha_0 * alpha_0 + 4 * R * alpha_1 * alpha_max)) / (2 * alpha_1);
        }
        if (t_guess < 0) {
            std::cout << "Guess " << t_guess << " was negative" << std::endl;
            t_guess = (gradNorm > 1) ? 1.0 / gradNorm : 1.0 / sqrt(gradNorm);
        }
        t_guess = fmin(t_guess, 1.0 / gradNorm);

        std::cout << "p_dot norm = " << gradNorm << std::endl;
        std::cout << "p_dot G-norm = " << G_pdot_pdot << std::endl;
        std::cout << "p_ddot norm = " << p_dotdot.norm() << std::endl;
        std::cout << "K norm = " << K.norm() << std::endl;
        std::cout << "t_guess = " << t_guess << std::endl;
        std::cout << "alpha_0 = " << alpha_0 << std::endl;
        std::cout << "alpha_1 = " << alpha_1 << std::endl;

        Eigen::MatrixXd K_matrix(K.rows() / 3, 3);
        VectorXdIntoMatrix(K, K_matrix);
        Eigen::MatrixXd origPos = curves->positions;

        std::cout << "Dot product of pddot with K = " << mult->DotProduct(p_dot, K) << std::endl;
        setCircleStep(curves, origPos, projectedGradient, K_matrix, t_guess, R, alpha_0, alpha_1, p_dotdot_mat);
        Eigen::MatrixXd newPos = curves->positions;
        Eigen::MatrixXd diff = newPos - origPos;
        Eigen::VectorXd diffvec(curves->NumVertices(), 3);
        diffvec.setZero();
        MatrixIntoVectorX3(diff, diffvec);

        std::cout << "Distance between old and new pos = " << diff.norm() << std::endl;
        std::cout << "Radius of 'circle' = " << 1.0 / K.norm() << std::endl;
        std::cout << "Dot product between axis and diff = " << mult->DotProduct(diffvec, K) << std::endl;

        mult->SetBlockTreeMode(BlockTreeMode::Matrix3AndProjector);

        double initialEnergy = SpatialTree::TPEnergyBH(curves, bvh, alpha, beta);
        double newEnergy = initialEnergy;
        double sigma = 0.01;
        int numBacktracks = 0;
        double ls_step_threshold = 1e-10;

        // plotSteps(curves, origPos, projectedGradient, K_matrix, R, alpha_0, alpha_1, bvh, targetLengths, gradDot, dE);

        // Line search
        while (fabs(t_guess) > ls_step_threshold) {
            // Search along this direction
            setCircleStep(curves, origPos, projectedGradient, K_matrix, t_guess, R, alpha_0, alpha_1, p_dotdot_mat);
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

        double backproj_threshold = 1e-3;
        Eigen::MatrixXd correction(nVerts, 3);
        correction.setZero();
        bool first = true;

        // Backprojection
        while (fabs(t_guess) > ls_step_threshold) {
            break;
            if (first) {
                first = false;
            }
            else {
                setCircleStep(curves, origPos, projectedGradient, K_matrix, t_guess, R, alpha_0, alpha_1, p_dotdot_mat);
            }
            // Update the centers of mass to reflect the new positions
            bvh->recomputeCentersOfMass(curves);

            Eigen::VectorXd phi(curves->NumEdges() + 3);
            double maxBefore = curves->FillConstraintViolations(phi, targetLengths);        

            // Compute and apply the correction
            solver->template BackprojectMultigrid<Smoother>(curves, phi, correction);
            for (int i = 0; i < nVerts; i++) {
                CurveVertex* p = curves->GetVertex(i);
                Vector3 cur = p->Position();
                p->SetPosition(cur + SelectRow(correction, i));
            }

            // Add length violations to RHS
            double maxViolation = curves->FillConstraintViolations(phi, targetLengths);
            std::cout << "  Constraint: " << maxBefore << " -> " << maxViolation << std::endl;

            if (maxViolation < backproj_threshold) {
                break;
            }
            else {
                std::cout << "Couldn't backproject; decreased step to " << t_guess << std::endl;
                t_guess /= 2;
            }
        }

        std::cout << "t_guess = " << t_guess << " (" << numBacktracks << " backtracks)" << std::endl;
        std::cout << "Total length = " << curves->TotalLength() << std::endl;
        std::cout << "Energy: " << initialEnergy << " -> " << newEnergy << std::endl;

        return fabs(t_guess);
    }
}
