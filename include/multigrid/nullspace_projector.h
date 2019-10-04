#pragma once

#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <iostream>
#include "flow/gradient_constraints.h"

namespace LWS {

    class NullSpaceProjector {
        private:
        Eigen::SparseMatrix<double> B;
        Eigen::SparseMatrix<double> BBT;
        bool prefactored;
        Eigen::SparseLU<Eigen::SparseMatrix<double>, Eigen::COLAMDOrdering<int>> BBT_solver;

        public:
        template<typename T>
        NullSpaceProjector(GradientConstraints<T> &constraints) {
            constraints.FillConstraintMatrix(B);
            // Pre-factorize B*B^T
            Eigen::SparseMatrix<double> BBT = B * B.transpose();
            BBT_solver.analyzePattern(BBT);
            BBT_solver.factorize(BBT);
        }

        ~NullSpaceProjector() {}

        template<typename V, typename Dest>
        void ProjectToNullspace(V &v, Dest &out) {
            // We want to get (B^T) (B B^T)^{-1} B v, so start from the right
            out = B * v;
            // Multiply with (B B^T) inverse 
            out = BBT_solver.solve(out);
            // Lastly multiply with B^T
            out = B.transpose() * out;
            // The full operator is P = I - (B^T) (B B^T)^{-1} B, so
            // subtract the above from the input
            out = v - out;
        }

        template<typename V, typename Dest>
        void ApplyBPinv(V &v, Dest &out) {
            out = BBT_solver.solve(v);
            out = B.transpose() * out;
        }

        Eigen::VectorXd ProjectToNullspace(Eigen::VectorXd &v) {
            Eigen::VectorXd out(v.rows());
            ProjectToNullspace(v, out);
            return out;
        }

        Eigen::VectorXd EvaluateConstraints(Eigen::VectorXd &v) {
            return B * v;
        }
    };

    /*
        template<typename V, typename Dest>
        void Multiply(V &v, Dest &out) {
            // We want to get (B^T) (B B^T)^{-1} B v, so start from the right
            out = B * v;
            Eigen::SparseMatrix<double> BBT = B * B.transpose();
            // Multiply with (B B^T) inverse
            out = BBT.partialPivLu().solve(out);
            // Lastly multiply with B^T
            out = B.transpose() * out;
            // The full operator is P = I - (B^T) (B B^T)^{-1} B, so
            // subtract the above from the input
            out = v - out;
        }
     */

    /*
    void NullSpaceProjector::Multiply(Eigen::VectorXd &v, Eigen::VectorXd &out) {
        // We want to get (B^T) (B B^T)^{-1} B v, so start from the right
        out = B * v;
        Eigen::SparseMatrix<double> BBT = B * B.transpose();
        // Multiply with (B B^T) inverse            
        Eigen::SparseLU<Eigen::SparseMatrix<double>, Eigen::COLAMDOrdering<int>> BBT_solver;
        BBT_solver.analyzePattern(BBT);
        BBT_solver.factorize(BBT);
        out = BBT_solver.solve(out);
        // Lastly multiply with B^T
        out = B.transpose() * out;
        // The full operator is P = I - (B^T) (B B^T)^{-1} B, so
        // subtract the above from the input
        out = v - out;
    }
    */

}