#pragma once

#include <Eigen/Sparse>
#include "poly_curve.h"

namespace LWS {
    
    template<typename T>
    class GradientConstraints {
        public:
        void FillConstraintMatrix(Eigen::SparseMatrix<double> &B) {
            return static_cast<T&>(*this).FillConstraintMatrix(B);
        }
    };

    class BarycenterConstraint : public GradientConstraints<BarycenterConstraint> {
        private:
        PolyCurveGroup* curves;
        int nVerts;

        public:
        BarycenterConstraint(PolyCurveGroup* c) {
            curves = c;
            nVerts = curves->NumVertices();
        }

        void FillConstraintMatrix(Eigen::SparseMatrix<double> &B) {
            std::vector<Eigen::Triplet<double>> triplets;
            B.resize(1, nVerts);
            double totalLength = curves->TotalLength();
            // Fill a single row with normalized vertex weights
            for (int i = 0; i < nVerts; i++) {
                double wt = curves->GetCurvePoint(i).DualLength() / totalLength;
                triplets.push_back(Eigen::Triplet<double>(0, i, wt));
            }
            B.setFromTriplets(triplets.begin(), triplets.end());
        }
    };

}