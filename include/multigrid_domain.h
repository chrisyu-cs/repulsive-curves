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

        Mult* GetMultiplier() const {
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
        double alpha, beta;
        
        PolyCurveDomain(PolyCurveGroup* c, double a, double b) {
            curves = c;
            alpha = a;
            beta = b;
            int nVerts = c->NumVertices();

            Eigen::MatrixXd A = GetFullMatrix();
            multiplier = new DenseMatrixMult(A);
        }

        ~PolyCurveDomain() {
            delete multiplier;
        }

        MultigridDomain<PolyCurveDomain, DenseMatrixMult>* Coarsen(MultigridOperator &prolongOp, MultigridOperator &sparsifyOp) const {
            PolyCurveGroup* coarsened = curves->Coarsen(prolongOp, sparsifyOp);
            return new PolyCurveDomain(coarsened, alpha, beta);
        }

        DenseMatrixMult* GetMultiplier() const {
            return multiplier;
        }

        Eigen::MatrixXd GetFullMatrix() const {
            int nVerts = curves->NumVertices();
            // Eigen::MatrixXd A = TestMatrices::CurveLaplacian(nVerts);
            Eigen::MatrixXd A(nVerts, nVerts);
            A.setZero();
            SobolevCurves::SobolevLengthScaled(curves, alpha, beta, A, 0.01);
            return A;
        }

        int NumVertices() const {
            return curves->NumVertices();
        }
    };

    class Interval1DDomain : public MultigridDomain<Interval1DDomain, DenseMatrixMult> {
        public:
        int nVerts;
        DenseMatrixMult* multiplier;

        Interval1DDomain(int n) {
            nVerts = n;
            multiplier = new DenseMatrixMult(GetFullMatrix());
        }

        ~Interval1DDomain() {
            delete multiplier;
        }

        DenseMatrixMult* GetMultiplier() const {
            return multiplier;
        }

        Eigen::MatrixXd GetFullMatrix() const {
            return TestMatrices::LaplacianNeumann1D(nVerts);
        }

        MultigridDomain<Interval1DDomain, DenseMatrixMult>* Coarsen(MultigridOperator &prolongOp,
                                                                    MultigridOperator &restrictOp) const {
            
            bool isOdd = (nVerts % 2) == 1;
            int sparseVerts = nVerts / 2 + 1;
            
            std::vector<Eigen::Triplet<double>> triplets;
            std::vector<Eigen::Triplet<double>> coarsenTriplets;

            for (int i = 0; i < nVerts; i++) {
                if (i == 0) {
                    triplets.push_back(Eigen::Triplet<double>(0, 0, 1));
                    coarsenTriplets.push_back(Eigen::Triplet<double>(0, 0, 1));
                }
                else if (i == nVerts - 1) {
                    triplets.push_back(Eigen::Triplet<double>(i, sparseVerts - 1, 1));
                    coarsenTriplets.push_back(Eigen::Triplet<double>(i, sparseVerts - 1, 1));
                }
                else {
                    int oldI = i / 2;
                    if (i % 2 == 0) {
                        triplets.push_back(Eigen::Triplet<double>(i, oldI - 1, 0.25));
                        triplets.push_back(Eigen::Triplet<double>(i, oldI, 0.5));
                        triplets.push_back(Eigen::Triplet<double>(i, oldI + 1, 0.25));

                        coarsenTriplets.push_back(Eigen::Triplet<double>(i, oldI, 0.5));
                    }
                    else {
                        triplets.push_back(Eigen::Triplet<double>(i, oldI, 0.5));
                        triplets.push_back(Eigen::Triplet<double>(i, oldI + 1, 0.5));

                        coarsenTriplets.push_back(Eigen::Triplet<double>(i, oldI, 0.25));
                        coarsenTriplets.push_back(Eigen::Triplet<double>(i, oldI + 1, 0.25));
                    }
                }
            }

            Eigen::SparseMatrix<double> prolongMatrix;
            prolongMatrix.resize(nVerts, sparseVerts);
            prolongMatrix.setFromTriplets(triplets.begin(), triplets.end());
            prolongOp.matrices.push_back(IndexedMatrix{prolongMatrix, 0, 0});

            Eigen::SparseMatrix<double> restrictMatrix;
            restrictMatrix.resize(nVerts, sparseVerts);
            restrictMatrix.setFromTriplets(coarsenTriplets.begin(), coarsenTriplets.end());
            restrictOp.matrices.push_back(IndexedMatrix{restrictMatrix, 0, 0});

            return new Interval1DDomain(sparseVerts);
        }
        
        int NumVertices() const {
            return nVerts;
        }
    };

}