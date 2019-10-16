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
        static Eigen::MatrixXd DifferentiateRHS(PolyCurveNetwork* curves, Eigen::MatrixXd &projectedGradient,
        Eigen::MatrixXd &l2Gradient, BVHNode3D* bvh, BlockClusterTree* mult, double alpha, double beta, double epsilon);

        template<typename Solver, typename Smoother>
        static Eigen::MatrixXd SolveSecondDerivative(PolyCurveNetwork* curves, Eigen::MatrixXd &projectedGradient,
        Eigen::MatrixXd &l2Gradient, BVHNode3D* bvh, Solver* solver, double alpha, double beta, double epsilon);

        template<typename Solver, typename Smoother>
        static void CircleSearchStep(PolyCurveNetwork* curves, Eigen::MatrixXd &projectedGradient,
        Eigen::MatrixXd &l2Gradient, BVHNode3D* bvh, Solver* solver, double alpha, double beta, double epsilon);

    };

    template<typename Solver, typename Smoother>
    Eigen::MatrixXd CircleSearch::SolveSecondDerivative(PolyCurveNetwork* curves, Eigen::MatrixXd &projectedGradient,
    Eigen::MatrixXd &l2Gradient, BVHNode3D* bvh, Solver* solver, double alpha, double beta, double epsilon) {
        BlockClusterTree* mult = solver->levels[0]->GetMultiplier();
        Eigen::MatrixXd rhs = DifferentiateRHS(curves, projectedGradient, l2Gradient, bvh, mult, alpha, beta, epsilon);

        Eigen::VectorXd rhs3x;
        rhs3x.setZero(rhs.rows() * 3);
        MatrixIntoVectorX3(rhs, rhs3x);

        rhs3x = curves->constraintProjector->ProjectToNullspace(rhs3x);
        rhs3x = solver->template VCycleSolve<Smoother>(rhs3x);

        rhs.setZero();
        VectorXdIntoMatrix(rhs3x, rhs);
        return rhs;
    }

    inline double alphaOfStep(double t, double R, double alpha_0, double alpha_1) {
        return (1.0 / R) * (alpha_0 * t + alpha_1 * t * t);
    }

    inline void setStep(PolyCurveNetwork* curves, Eigen::MatrixXd &origPos, Eigen::MatrixXd &p_dot, Eigen::MatrixXd &K,
    double t, double R, double alpha_0, double alpha_1) {
        double alpha_of_t = alphaOfStep(t, R, alpha_0, alpha_1);
        double p_dot_coeff = sin(alpha_of_t) / alpha_0;
        double K_coeff = (1.0 - cos(alpha_of_t));

        curves->positions = origPos + R * (p_dot_coeff * p_dot + R * K_coeff * K);
    }

    template<typename Solver, typename Smoother>
    void CircleSearch::CircleSearchStep(PolyCurveNetwork* curves, Eigen::MatrixXd &projectedGradient,
    Eigen::MatrixXd &l2Gradient, BVHNode3D* bvh, Solver* solver, double alpha, double beta, double epsilon) {
        Eigen::MatrixXd pdd_matrix = SolveSecondDerivative(curves, projectedGradient, l2Gradient,
            bvh, solver, alpha, beta, epsilon);
        
        Eigen::VectorXd p_dot(projectedGradient.rows() * 3);
        p_dot.setZero();
        MatrixIntoVectorX3(projectedGradient, p_dot);

        Eigen::VectorXd p_dotdot(secondDeriv.rows() * 3);
        p_dotdot.setZero();
        MatrixIntoVectorX3(pdd_matrix, p_dotdot);
        
        BlockClusterTree* mult = solver->levels[0]->GetMultiplier();
        mult->SetBlockTreeMode(BlockTreeMode::Matrix3Only);
        double G_pdot_pdot = mult->DotProduct(p_dot, p_dot);
        double G_pdot_pdotdot = mult->DotProduct(p_dot, p_dotdot);

        Eigen::VectorXd K = (1.0 / G_pdot_pdot) * (p_dotdot - p_dot * (G_pdot_pdotdot / G_pdot_pdot));
        double R = mult->DotProduct(K, K);
        double alpha_0 = sqrt(G_pdot_pdot);
        double alpha_1 = 0.5 * (G_pdot_pdotdot / alpha_0);

        double t = -2 * (alpha_0 / alpha_1);

        Eigen::MatrixXd K_matrix(K.rows() / 3, 3);
        VectorXdIntoMatrix(K, K_matrix);
        Eigen::MatrixXd origPos = curves->positions;

        // TODO: actually search along this direction
        setStep(curves, origPos, projectedGradient, K_matrix, t, R, alpha_0, alpha_1);

        mult->SetBlockTreeMode(BlockTreeMode::Matrix3AndProjector);
    }
}
