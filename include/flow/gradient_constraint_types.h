#pragma once

#include <Eigen/Sparse>
#include "gradient_constraints.h"
#include "poly_curve_network.h"

namespace LWS {

    class EdgeLengthConstraint : public GradientConstraints<EdgeLengthConstraint> {
        private:
        PolyCurveNetwork* curves;

        public:
        EdgeLengthConstraint(PolyCurveNetwork* c) {
            curves = c;
        }

        void FillConstraintMatrix(Eigen::SparseMatrix<double> &B) {
            std::vector<Eigen::Triplet<double>> triplets;
            int nVerts = curves->NumVertices();
            int nEdges = curves->NumEdges();

            // Add the barycenter row
            double totalLength = curves->TotalLength();
            // Fill a single row with normalized vertex weights
            for (int i = 0; i < nVerts; i++) {
                double wt = curves->GetVertex(i)->DualLength() / totalLength;
                triplets.push_back(Eigen::Triplet<double>(0, 3 * i, wt));
                triplets.push_back(Eigen::Triplet<double>(1, 3 * i + 1, wt));
                triplets.push_back(Eigen::Triplet<double>(2, 3 * i + 2, wt));
            }

            // Add the edge length rows
            for (int i = 0; i < nEdges; i++) {
                CurveEdge* edge = curves->GetEdge(i);
                CurveVertex* pt1 = edge->prevVert;
                CurveVertex* pt2 = edge->nextVert;

                // This is the gradient of edge length wrt pt1; the gradient wrt pt2 is just negative of this.
                Vector3 grad1 = pt1->Position() - pt2->Position();
                grad1 = grad1.normalize();

                int j1 = pt1->GlobalIndex();
                int j2 = pt2->GlobalIndex();
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

    class BarycenterConstraint3X : public GradientConstraints<BarycenterConstraint3X> {
        private:
        PolyCurveNetwork* curves;
        int nVerts;

        public:
        BarycenterConstraint3X(PolyCurveNetwork* c) {
            curves = c;
            nVerts = curves->NumVertices();
        }

        void FillConstraintMatrix(Eigen::SparseMatrix<double> &B) {
            std::vector<Eigen::Triplet<double>> triplets;
            B.resize(3, 3 * nVerts);
            double totalLength = curves->TotalLength();
            // Fill a single row with normalized vertex weights
            for (int i = 0; i < nVerts; i++) {
                double wt = curves->GetVertex(i)->DualLength() / totalLength;
                triplets.push_back(Eigen::Triplet<double>(0, 3 * i, wt));
                triplets.push_back(Eigen::Triplet<double>(1, 3 * i + 1, wt));
                triplets.push_back(Eigen::Triplet<double>(2, 3 * i + 2, wt));
            }
            B.setFromTriplets(triplets.begin(), triplets.end());
        }
    };

    class BarycenterConstraint : public GradientConstraints<BarycenterConstraint> {
        private:
        PolyCurveNetwork* curves;
        int nVerts;

        public:
        BarycenterConstraint(PolyCurveNetwork* c) {
            curves = c;
            nVerts = curves->NumVertices();
        }

        void FillConstraintMatrix(Eigen::SparseMatrix<double> &B) {
            std::vector<Eigen::Triplet<double>> triplets;
            B.resize(1, nVerts);
            double totalLength = curves->TotalLength();
            // Fill a single row with normalized vertex weights
            for (int i = 0; i < nVerts; i++) {
                double wt = curves->GetVertex(i)->DualLength() / totalLength;
                triplets.push_back(Eigen::Triplet<double>(0, i, wt));
            }
            B.setFromTriplets(triplets.begin(), triplets.end());
        }
    };
}