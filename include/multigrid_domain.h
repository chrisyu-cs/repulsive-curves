#pragma once

#include "poly_curve.h"
#include "multigrid_operator.h"
#include "product/matrix_free.h"
#include "product/dense_matrix.h"
#include "product/test_matrices.h"

#include <Eigen/Core>

namespace LWS {

    template<typename T, typename Mult>
    class MultigridDomain {
        public:

        MultigridDomain<T, Mult>* Coarsen(MultigridOperator &prolongOp, MultigridOperator &sparsifyOp) const {
            return static_cast<T const&>(*this).Coarsen(prolongOp, sparsifyOp);
        }

        Product::MatrixReplacement<Mult> GetMultiplier() const {
            return static_cast<T const&>(*this).GetMultiplier();
        }

        Eigen::MatrixXd GetFullMatrix() const {
            return static_cast<T const&>(*this).GetFullMatrix() ;
        }
        
        int NumVertices() const {
            return static_cast<T const&>(*this).NumVertices();
        }
    };

    class PolyCurveDomain : public MultigridDomain<PolyCurveDomain, DenseMatrixMult> {
        public:
        PolyCurveGroup* curves;
        DenseMatrixMult* multiplier;
        WrappedMatrix hMatrix;
        
        PolyCurveDomain(PolyCurveGroup* c) {
            curves = c;
            int nVerts = c->NumVertices();
            multiplier = new DenseMatrixMult(TestMatrices::LaplacianDirichlet1D(nVerts));
            hMatrix = WrappedMatrix(multiplier, nVerts);
        }

        ~PolyCurveDomain() {
            delete multiplier;
        }

        MultigridDomain<PolyCurveDomain, DenseMatrixMult>* Coarsen(MultigridOperator &prolongOp, MultigridOperator &sparsifyOp) const {
            PolyCurveGroup* coarsened = curves->Coarsen(prolongOp, sparsifyOp);
            return new PolyCurveDomain(coarsened);
        }

        Product::MatrixReplacement<DenseMatrixMult> GetMultiplier() const {
            return hMatrix;
        }

        Eigen::MatrixXd GetFullMatrix() const {
            int nVerts = curves->NumVertices();
            return TestMatrices::LaplacianDirichlet1D(nVerts);
        }

        int NumVertices() const {
            return curves->NumVertices();
        }
    };

    class Interval1DDomain : public MultigridDomain<Interval1DDomain, DenseMatrixMult> {
        public:
        int nVerts;
        DenseMatrixMult* multiplier;
        WrappedMatrix hMatrix;

        Interval1DDomain(int n) {
            nVerts = n;
            multiplier = new DenseMatrixMult(TestMatrices::LaplacianDirichlet1D(nVerts));
            hMatrix = WrappedMatrix(multiplier, nVerts);
        }

        ~Interval1DDomain() {
            delete multiplier;
        }

        MultigridDomain<Interval1DDomain, DenseMatrixMult>* Coarsen(MultigridOperator &prolongOp, MultigridOperator &sparsifyOp) const {
            
            bool isOdd = (nVerts % 2) == 1;
            int sparseVerts = nVerts / 2 + 1;
            
            std::vector<Eigen::Triplet<double>> triplets;

            if (isOdd) {
                for (int i = 0; i < sparseVerts; i++) {
                    if (i == 0) {
                        triplets.push_back(Eigen::Triplet<double>(0, i, 1));
                        triplets.push_back(Eigen::Triplet<double>(1, i, 0.5));
                    }
                    else if (i == sparseVerts - 1) {
                        triplets.push_back(Eigen::Triplet<double>(nVerts - 1, i, 1));
                        triplets.push_back(Eigen::Triplet<double>(nVerts - 2, i, 0.5));
                    }
                    else {
                        int oldI = 2 * i;
                        triplets.push_back(Eigen::Triplet<double>(oldI - 1, i, 0.5));
                        triplets.push_back(Eigen::Triplet<double>(oldI, i, 1));
                        triplets.push_back(Eigen::Triplet<double>(oldI + 1, i, 0.5));
                    }
                }
            }
            else {
                for (int i = 0; i < sparseVerts; i++) {
                    if (i == 0) {
                        triplets.push_back(Eigen::Triplet<double>(0, i, 1));
                        triplets.push_back(Eigen::Triplet<double>(1, i, 0.5));
                    }
                    else if (i == sparseVerts - 2) {
                        triplets.push_back(Eigen::Triplet<double>(nVerts - 3, i, 0.5));
                        triplets.push_back(Eigen::Triplet<double>(nVerts - 2, i, 1));
                    }
                    else if (i == sparseVerts - 1) {
                        triplets.push_back(Eigen::Triplet<double>(nVerts - 1, i, 1));
                    }
                    else {
                        int oldI = 2 * i;
                        triplets.push_back(Eigen::Triplet<double>(oldI - 1, i, 0.5));
                        triplets.push_back(Eigen::Triplet<double>(oldI, i, 1));
                        triplets.push_back(Eigen::Triplet<double>(oldI + 1, i, 0.5));
                    }
                }
            }

            Eigen::SparseMatrix<double> prolongMatrix;
            prolongMatrix.resize(nVerts, sparseVerts);
            prolongMatrix.setFromTriplets(triplets.begin(), triplets.end());
            prolongOp.matrices.push_back(IndexedMatrix{prolongMatrix, 0, 0});

            return new Interval1DDomain(sparseVerts);
        }

        Product::MatrixReplacement<DenseMatrixMult> GetMultiplier() const {
            return hMatrix;
        }

        Eigen::MatrixXd GetFullMatrix() const {
            return TestMatrices::LaplacianDirichlet1D(nVerts);
        }
        
        int NumVertices() const {
            return nVerts;
        }
    };

}