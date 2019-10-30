#pragma once

#include <Eigen/Core>
#include <Eigen/Sparse>
#include <iostream>

namespace LWS {
    
    template<typename T>
    class GradientConstraints {
        public:
        int NumConstraintRows() const {
            return static_cast<const T&>(*this).NumConstraintRows();
        }

        int NumExpectedCols() const {
            return static_cast<const T&>(*this).NumExpectedCols();
        }

        int SaddleNumRows() const {
            return NumConstraintRows() + NumExpectedCols();
        }

        void FillConstraintMatrix(Eigen::SparseMatrix<double> &B) const {
            std::vector<Eigen::Triplet<double>> triplets;
            int nRows = NumConstraintRows();
            int nCols = NumExpectedCols();
            static_cast<const T&>(*this).AddTriplets(triplets);

            B.resize(nRows, nCols);
            B.setFromTriplets(triplets.begin(), triplets.end());
        }

        void FillDenseBlock(Eigen::MatrixXd &A) const {
            std::vector<Eigen::Triplet<double>> triplets;
            static_cast<const T&>(*this).AddTriplets(triplets);
            int offset = NumExpectedCols();

            for (auto &t : triplets) {
                // Copy into lower-left block
                A(offset + t.row(), t.col()) = t.value();
                // Copy transpose into upper-right block
                A(t.col(), offset + t.row()) = t.value();
            }
        }
        
        void FillConstraintViolations(Eigen::VectorXd &b, Eigen::VectorXd &targets) {
            // TODO: fill values of constraint function in a generic way
        }
    };
}