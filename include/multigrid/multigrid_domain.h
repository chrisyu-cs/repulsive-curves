#pragma once

#include "spatial/tpe_bvh.h"
#include "poly_curve.h"
#include "multigrid_operator.h"
#include "product/matrix_free.h"
#include "product/dense_matrix.h"
#include "product/test_matrices.h"
#include "nullspace_projector.h"

#include <Eigen/Core>

namespace LWS {

    template<typename T, typename Mult>
    class MultigridDomain {
        public:
        typedef Mult MultType;

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
        
        NullSpaceProjector* GetConstraintProjector() const {
            return static_cast<T const&>(*this).GetConstraintProjector();
        }
    };

    class PolyCurveHMatrixDomain : public MultigridDomain<PolyCurveHMatrixDomain, BlockClusterTree> {
        public:
        PolyCurveGroup* curves;
        BVHNode3D* bvh;
        BlockClusterTree* tree;
        double alpha, beta;
        double epsilon;
        double sepCoeff;
        BlockTreeMode mode;

        PolyCurveHMatrixDomain(PolyCurveGroup* c, double a, double b, double sep, double e, BlockTreeMode m) {
            curves = c;
            alpha = a;
            beta = b;
            epsilon = e;
            sepCoeff = sep;
            int nVerts = c->NumVertices();

            bvh = CreateEdgeBVHFromCurve(curves);
            tree = new BlockClusterTree(curves, bvh, sepCoeff, alpha, beta, epsilon);
            mode = m;
            tree->SetBlockTreeMode(mode);
        }

        ~PolyCurveHMatrixDomain() {
            delete tree;
            delete bvh;
        }

        MultigridDomain<PolyCurveHMatrixDomain, BlockClusterTree>* Coarsen(MultigridOperator &prolongOp, MultigridOperator &sparsifyOp) const {
            PolyCurveGroup* coarsened = curves->Coarsen(prolongOp, sparsifyOp);
            return new PolyCurveHMatrixDomain(coarsened, alpha, beta, sepCoeff, epsilon, mode);
        }

        BlockClusterTree* GetMultiplier() const {
            return tree;
        }

        Eigen::MatrixXd GetFullMatrix() const {
            int rows = (mode == BlockTreeMode::Barycenter) ? curves->NumVertices() + 1 : curves->NumVertices();
            Eigen::MatrixXd A;
            A.setZero(rows, rows);
            if (mode == BlockTreeMode::Barycenter) {
                SobolevCurves::SobolevPlusBarycenter(curves, alpha, beta, A);
            }
            else {
                SobolevCurves::SobolevGramMatrix(curves, alpha, beta, A, epsilon);
            }
            return A;
        }

        int NumVertices() const {
            return curves->NumVertices();
        }

        NullSpaceProjector* GetConstraintProjector() const {
            return curves->constrP;
        }
    };

    class PolyCurveDenseDomain : public MultigridDomain<PolyCurveDenseDomain, DenseMatrixMult> {
        public:
        PolyCurveGroup* curves;
        DenseMatrixMult* multiplier;
        double alpha, beta;
        double epsilon;
        
        PolyCurveDenseDomain(PolyCurveGroup* c, double a, double b, double e) {
            curves = c;
            alpha = a;
            beta = b;
            epsilon = e;
            int nVerts = c->NumVertices();

            Eigen::MatrixXd A = GetFullMatrix();
            multiplier = new DenseMatrixMult(A);
        }

        ~PolyCurveDenseDomain() {
            delete multiplier;
        }

        MultigridDomain<PolyCurveDenseDomain, DenseMatrixMult>* Coarsen(MultigridOperator &prolongOp, MultigridOperator &sparsifyOp) const {
            PolyCurveGroup* coarsened = curves->Coarsen(prolongOp, sparsifyOp);
            return new PolyCurveDenseDomain(coarsened, alpha, beta, epsilon);
        }

        DenseMatrixMult* GetMultiplier() const {
            return multiplier;
        }

        Eigen::MatrixXd GetFullMatrix() const {
            int nVerts = curves->NumVertices();
            Eigen::MatrixXd A;
            A.setZero(nVerts, nVerts);
            SobolevCurves::SobolevGramMatrix(curves, alpha, beta, A, epsilon);
            // Eigen::MatrixXd A = TestMatrices::CurveMetricLaplacian(curves, 0.01);
            return A;
        }

        int NumVertices() const {
            return curves->NumVertices();
        }

        NullSpaceProjector* GetConstraintProjector() const {
            return 0;
        }
    };
}