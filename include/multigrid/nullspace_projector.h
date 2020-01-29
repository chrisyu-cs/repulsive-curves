#pragma once

#include <Eigen/Dense>
#include <Eigen/Sparse>
#include "domain_constraints.h"

namespace LWS {

    class NullSpaceProjector {
        private:
        Eigen::SparseMatrix<double> B;
        Eigen::SparseMatrix<double> BBT;
        bool prefactored;
        Eigen::SparseLU<Eigen::SparseMatrix<double>, Eigen::COLAMDOrdering<int>> BBT_solver;

        public:
        template<typename T>
        NullSpaceProjector(DomainConstraints<T> &constraints) {
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

        inline Eigen::MatrixXd ProjectorMatrix() {
            // We want to get (B^T) (B B^T)^{-1} B v, so start from the right
            Eigen::MatrixXd BBT_inv = (B * B.transpose()).toDense().inverse();
            BBT_inv = B.transpose() * BBT_inv * B;
            return BBT_inv;
        }

        // Multiply by (BB^T)^(-1).
        template<typename V, typename Dest>
        void SolveBBT(V &v, Dest &out) {
            out = BBT_solver.solve(v);
        }

        // Multiply by B^(dagger).
        template<typename V, typename Dest>
        void ApplyBPinv(V &v, Dest &out) {
            out = BBT_solver.solve(v);
            out = B.transpose() * out;
        }

        // Multiply by B^(dagger).
        template<typename V>
        Eigen::VectorXd ApplyBPinv(V &v) {
            Eigen::VectorXd out = BBT_solver.solve(v);
            return B.transpose() * out;
        }

        // Multiply by B.
        template<typename V, typename Dest>
        void ApplyB(V &v, Dest &out) {
            out = B * v;
        }

        // Multiply by B.
        template<typename V>
        Eigen::VectorXd ApplyB(V &v) {
            return B * v;
        }

        // Multiply by B transposed.
        template<typename V, typename Dest>
        void ApplyBTranspose(V &v, Dest &out) {
            out = B.transpose() * v;
        }

        // Multiply by B transposed.
        template<typename V>
        Eigen::VectorXd ApplyBTranspose(V &v) {
            return B.transpose() * v;
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
}