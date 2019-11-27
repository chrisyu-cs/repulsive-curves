#pragma once

#include <Eigen/Core>
#include "tpe_energy_sc.h"
#include "sobo_slobo.h"
#include "poly_curve_network.h"
#include "flow/gradient_constraints.h"

namespace LWS {

    class CircleSearch {
        // Compute the right-hand side of the second derivative system,
        // which is comprised of directional derivative of energy gradient,
        // and (directional derivative of Gram matrix) * energy gradient.
        public:
        CircleSearch(PolyCurveNetwork* c, double a, double b, double eps);

        template<typename T>
        Eigen::VectorXd DifferentiateRHS(Eigen::VectorXd &projectedGradient,
            GradientConstraints<T> &constraints, Eigen::MatrixXd &A);

        template<typename T>
        Eigen::VectorXd SolveSecondDerivative(Eigen::VectorXd &projectedGradient,
            GradientConstraints<T> &constraints, Eigen::MatrixXd &l2_dE, Eigen::MatrixXd &A);

        template<typename T>
        double CircleSearchStep(Eigen::MatrixXd &projectedGradient, GradientConstraints<T> &constraints,
            Eigen::MatrixXd &l2_dE, Eigen::MatrixXd &A, Eigen::VectorXd &constraintTargets);

        PolyCurveNetwork* curves;
        double alpha;
        double beta;
        double epsilon;

    };

    inline void printVectorWithPrecision(Eigen::VectorXd &v, int precision) {
        std::cout.precision(precision);
        for (int i = 0; i < v.rows(); i++) {
            std::cout << v(i) << std::endl;
        }
    }

    inline void printVectorWithPrecision(Eigen::MatrixXd &v, int precision) {
        Eigen::VectorXd vec(v.rows() * 3);
        MatrixIntoVectorX3(v, vec);
        std::cout.precision(precision);
        for (int i = 0; i < vec.rows(); i++) {
            std::cout << vec(i) << std::endl;
        }
    }

    template<typename T>
    Eigen::VectorXd CircleSearch::SolveSecondDerivative(Eigen::VectorXd &projectedGradient,
    GradientConstraints<T> &constraints, Eigen::MatrixXd &l2_dE, Eigen::MatrixXd &A) {
        int nVerts = curves->NumVertices();
        auto lu = A.partialPivLu();
        Eigen::VectorXd rhs = DifferentiateRHS(projectedGradient, constraints, A);
        Eigen::VectorXd secondDeriv = lu.solve(rhs);
        // std::cout << secondDeriv << std::endl;
        std::cout << secondDeriv.rows() << " rows" << std::endl;
        std::cout << "Curve has " << nVerts << " vertices" << std::endl;

        // printVectorWithPrecision(secondDeriv, 15);

        return secondDeriv;
    }

    template<typename T>
    double CircleSearch::CircleSearchStep(Eigen::MatrixXd &projectedGradient, GradientConstraints<T> &constraints,
    Eigen::MatrixXd &l2_dE, Eigen::MatrixXd &A, Eigen::VectorXd &constraintTargets) {
        return 0;
    }

    template<typename T>
    Eigen::VectorXd CircleSearch::DifferentiateRHS(Eigen::VectorXd &projectedGradient,
    GradientConstraints<T> &constraints, Eigen::MatrixXd &A) {
        int nVerts = curves->NumVertices();
        int nEdges = curves->NumEdges();

        Eigen::MatrixXd origPositions = curves->positions;

        Eigen::MatrixXd pdot_mat(nVerts, 3);
        pdot_mat.setZero();
        for (int i = 0; i < nVerts; i++) {
            Vector3 p_i{projectedGradient(3 * i), projectedGradient(3 * i + 1), projectedGradient(3 * i + 2)};
            SetRow(pdot_mat, i, p_i);
        }

        // Space for the difference quotient for (derivative Gram) * gradient
        Eigen::VectorXd gramTimesGradient, gramTimesGradientEps;
        gramTimesGradient.setZero(A.rows());
        gramTimesGradientEps.setZero(A.rows());

        // Multiply Gram * gradient now
        gramTimesGradient = A * projectedGradient;

        Eigen::MatrixXd l2Gradient;
        l2Gradient.setZero(nVerts, 3);
        TPESC::FillGradientVectorDirect(curves, l2Gradient, alpha, beta);

        // We have the initial l2 gradient already, so now we just need
        // the epsilon-perturbed vectors
        double projNorm = projectedGradient.norm();

        // Move by epsilon in the search direction
        curves->positions -= epsilon * pdot_mat;

        Eigen::VectorXd posVec;
        posVec.setZero(3 * nVerts);
        MatrixIntoVectorX3(curves->positions, posVec);
        // printVectorWithPrecision(posVec, 15);
        
        // Get the epsilon-perturbed L2 gradient
        Eigen::MatrixXd l2GradientEps;
        l2GradientEps.setZero(nVerts, 3);
        TPESC::FillGradientVectorDirect(curves, l2GradientEps, alpha, beta);

        // Multiply (Gram + eps) * gradient
        Eigen::MatrixXd A_eps;
        A_eps.setZero(A.rows(), A.cols());
        SobolevCurves::Sobolev3XWithConstraints(curves, constraints, alpha, beta, A_eps);
        gramTimesGradientEps = A_eps * projectedGradient;

        std::cout << "Norm of diff As = " << (A - A_eps).norm() << std::endl;

        // Numerical derivative of Gram * gradient
        gramTimesGradientEps = (gramTimesGradientEps - gramTimesGradient) / epsilon;

        printVectorWithPrecision(l2GradientEps, 15);

        // Numerical derivative of L2 gradient
        l2GradientEps = (l2GradientEps - l2Gradient) / epsilon;
        Eigen::VectorXd sum;
        sum.setZero(A.rows());
        MatrixIntoVectorX3(l2GradientEps, sum);


        // Reset original values
        curves->positions = origPositions;
        
        sum = sum - gramTimesGradientEps;
        return sum;
    }
}
