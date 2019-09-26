#pragma once

#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <iostream>

namespace LWS {

    class NullSpaceProjector {
        private:
        Eigen::MatrixXd B;

        public:
        NullSpaceProjector(Eigen::MatrixXd constraints) {
            B = constraints;
        }

        ~NullSpaceProjector() {}

        template<typename V, typename Dest>
        void Multiply(V &v, Dest &out) {
            // We want to get (B^T) (B B^T)^{-1} B v, so start from the right
            out = B * v;
            Eigen::MatrixXd BBT = B * B.transpose();
            // Multiply with (B B^T) inverse
            out = BBT.partialPivLu().solve(out);
            // Lastly multiply with B^T
            out = B.transpose() * out;
            // The full operator is P = I - (B^T) (B B^T)^{-1} B, so
            // subtract the above from the input
            out = v - out;
        }

        Eigen::VectorXd Multiply(Eigen::VectorXd &v) {
            Eigen::VectorXd out(v.rows());
            Multiply(v, out);
            return out;
        }
    };

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