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

    class EdgeLengthConstraint : public GradientConstraints<EdgeLengthConstraint> {
        private:
        PolyCurveGroup* curves;
        int nVerts;

        public:
        EdgeLengthConstraint(PolyCurveGroup* c) {
            curves = c;
            nVerts = curves->NumVertices();
        }

        void FillConstraintMatrix(Eigen::SparseMatrix<double> &B) {
            std::vector<Eigen::Triplet<double>> triplets;

            // Add the barycenter row
            double totalLength = curves->TotalLength();
            // Fill a single row with normalized vertex weights
            for (int i = 0; i < nVerts; i++) {
                double wt = curves->GetCurvePoint(i).DualLength() / totalLength;
                triplets.push_back(Eigen::Triplet<double>(0, 3 * i, wt));
                triplets.push_back(Eigen::Triplet<double>(1, 3 * i + 1, wt));
                triplets.push_back(Eigen::Triplet<double>(2, 3 * i + 2, wt));
            }

            // Add the edge length rows
            for (int i = 0; i < nVerts; i++) {
                PointOnCurve pt1 = curves->GetCurvePoint(i);
                PointOnCurve pt2 = pt1.Next();

                // This is the gradient of edge length wrt pt1; the gradient wrt pt2 is just negative of this.
                Vector3 grad1 = pt1.Position() - pt2.Position();
                grad1 = grad1.normalize();

                int j1 = curves->GlobalIndex(pt1);
                int j2 = curves->GlobalIndex(pt2);
                int start = i + 3;

                // Write the three gradient entries for pt1 into the row
                triplets.push_back(Eigen::Triplet<double>(start, 3 * j1,     grad1.x));
                triplets.push_back(Eigen::Triplet<double>(start, 3 * j1 + 1, grad1.y));
                triplets.push_back(Eigen::Triplet<double>(start, 3 * j1 + 2, grad1.z));

                // Similarly write the three gradient entries for pt2 into the same row
                triplets.push_back(Eigen::Triplet<double>(start, 3 * j2,     -grad1.x));
                triplets.push_back(Eigen::Triplet<double>(start, 3 * j2 + 1, -grad1.y));
                triplets.push_back(Eigen::Triplet<double>(start, 3 * j2 + 2, -grad1.z));
            }

            B.resize(nVerts + 3, 3 * nVerts);
            B.setFromTriplets(triplets.begin(), triplets.end());
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