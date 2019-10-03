#pragma once

#include <Eigen/Sparse>

namespace LWS {
    
    template<typename T>
    class GradientConstraints {
        public:
        void FillConstraintMatrix(Eigen::SparseMatrix<double> &B) {
            return static_cast<T&>(*this).FillConstraintMatrix(B);
        }
    };

}