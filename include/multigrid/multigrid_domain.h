#pragma once

#include "spatial/tpe_bvh.h"
#include "poly_curve.h"
#include "multigrid_operator.h"
#include "product/matrix_free.h"
#include "product/dense_matrix.h"
#include "product/test_matrices.h"
#include "nullspace_projector.h"
#include "flow/gradient_constraints.h"

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
            return static_cast<T const&>(*this).GetFullMatrix();
        }

        Eigen::VectorXd DirectSolve(Eigen::VectorXd &b) const {
            return static_cast<T const&>(*this).DirectSolve(b);
        }

        int NumVertices() const {
            return static_cast<T const&>(*this).NumVertices();
        }
        
        int NumRows() const {
            return static_cast<T const&>(*this).NumRows();
        }

        MultigridMode GetMode() const {
            return static_cast<T const&>(*this).GetMode();
        }
        
        NullSpaceProjector* GetConstraintProjector() const {
            return static_cast<T const&>(*this).GetConstraintProjector();
        }
    };

    class PolyCurveNullProjectorDomain : public MultigridDomain<PolyCurveNullProjectorDomain, BlockClusterTree> {
        public:
        PolyCurveGroup* curves;
        BVHNode3D* bvh;
        BlockClusterTree* tree;
        double alpha, beta;
        double sepCoeff;
        int nVerts;

        PolyCurveNullProjectorDomain(PolyCurveGroup* c, double a, double b, double sep) {
            curves = c;
            alpha = a;
            beta = b;
            sepCoeff = sep;
            nVerts = curves->NumVertices();

            bvh = CreateEdgeBVHFromCurve(curves);
            tree = new BlockClusterTree(curves, bvh, sepCoeff, alpha, beta);
            tree->SetBlockTreeMode(BlockTreeMode::MatrixAndProjector);

            curves->AddConstraints();
        }

        ~PolyCurveNullProjectorDomain() {
            delete tree;
            delete bvh;
        }

        MultigridDomain<PolyCurveNullProjectorDomain, BlockClusterTree>* Coarsen(MultigridOperator &prolongOp, MultigridOperator &sparsifyOp) const {
            PolyCurveGroup* coarsened = curves->Coarsen(prolongOp, sparsifyOp);
            return new PolyCurveNullProjectorDomain(coarsened, alpha, beta, sepCoeff);
        }

        BlockClusterTree* GetMultiplier() const {
            return tree;
        }

        Eigen::MatrixXd GetFullMatrix() const {
            int rows = nVerts + 1;
            Eigen::MatrixXd A;
            A.setZero(rows, rows);
            SobolevCurves::SobolevPlusBarycenter(curves, alpha, beta, A);
            return A;
        }

        Eigen::VectorXd DirectSolve(Eigen::VectorXd &b) const {
            Eigen::VectorXd b_aug(nVerts + 1);
            b_aug.block(0, 0, nVerts, 1) = b;
            b_aug(nVerts) = 0;

            Eigen::MatrixXd A = GetFullMatrix();
            b_aug = A.partialPivLu().solve(b_aug);
            Eigen::VectorXd x = b_aug.block(0, 0, nVerts, 1);
            return x;
        }

        int NumVertices() const {
            return nVerts;
        }

        int NumRows() const {
            return nVerts;
        }

        MultigridMode GetMode() const {
            return MultigridMode::MatrixAndProjector;
        }

        NullSpaceProjector* GetConstraintProjector() const {
            return curves->constraints;
        }
    };

    class PolyCurveSaddleDomain : public MultigridDomain<PolyCurveSaddleDomain, BlockClusterTree> {
        public:
        PolyCurveGroup* curves;
        BVHNode3D* bvh;
        BlockClusterTree* tree;
        double alpha, beta;
        double sepCoeff;
        int nVerts;

        PolyCurveSaddleDomain(PolyCurveGroup* c, double a, double b, double sep) {
            curves = c;
            alpha = a;
            beta = b;
            sepCoeff = sep;
            nVerts = curves->NumVertices();

            bvh = CreateEdgeBVHFromCurve(curves);
            tree = new BlockClusterTree(curves, bvh, sepCoeff, alpha, beta);
            tree->SetBlockTreeMode(BlockTreeMode::MatrixAndConstraints);
            BarycenterConstraint bc(curves);
            tree->SetConstraints(bc);
        }

        ~PolyCurveSaddleDomain() {
            delete tree;
            delete bvh;
        }

        MultigridDomain<PolyCurveSaddleDomain, BlockClusterTree>* Coarsen(MultigridOperator &prolongOp, MultigridOperator &sparsifyOp) const {
            PolyCurveGroup* coarsened = curves->Coarsen(prolongOp, sparsifyOp);
            return new PolyCurveSaddleDomain(coarsened, alpha, beta, sepCoeff);
        }

        BlockClusterTree* GetMultiplier() const {
            return tree;
        }

        Eigen::MatrixXd GetFullMatrix() const {
            int rows = nVerts + 1;
            Eigen::MatrixXd A;
            A.setZero(rows, rows);
            SobolevCurves::SobolevPlusBarycenter(curves, alpha, beta, A);
            return A;
        }

        Eigen::VectorXd DirectSolve(Eigen::VectorXd &b) const {
            Eigen::MatrixXd A = GetFullMatrix();
            return A.partialPivLu().solve(b);
        }

        int NumVertices() const {
            return nVerts;
        }

        int NumRows() const {
            return nVerts + 1;
        }

        MultigridMode GetMode() const {
            return MultigridMode::Barycenter;
        }

        NullSpaceProjector* GetConstraintProjector() const {
            return 0;
        }
    };

    class PolyCurveGramDomain : public MultigridDomain<PolyCurveGramDomain, BlockClusterTree> {
        public:
        PolyCurveGroup* curves;
        BVHNode3D* bvh;
        BlockClusterTree* tree;
        double alpha, beta;
        double epsilon;
        double sepCoeff;
        int nVerts;

        PolyCurveGramDomain(PolyCurveGroup* c, double a, double b, double sep, double e) {
            curves = c;
            alpha = a;
            beta = b;
            epsilon = e;
            sepCoeff = sep;
            nVerts = curves->NumVertices();

            bvh = CreateEdgeBVHFromCurve(curves);
            tree = new BlockClusterTree(curves, bvh, sepCoeff, alpha, beta, epsilon);
            tree->SetBlockTreeMode(BlockTreeMode::MatrixOnly);
        }

        ~PolyCurveGramDomain() {
            delete tree;
            delete bvh;
        }

        MultigridDomain<PolyCurveGramDomain, BlockClusterTree>* Coarsen(MultigridOperator &prolongOp, MultigridOperator &sparsifyOp) const {
            PolyCurveGroup* coarsened = curves->Coarsen(prolongOp, sparsifyOp);
            return new PolyCurveGramDomain(coarsened, alpha, beta, sepCoeff, epsilon);
        }

        BlockClusterTree* GetMultiplier() const {
            return tree;
        }

        Eigen::MatrixXd GetFullMatrix() const {
            int rows = curves->NumVertices();
            Eigen::MatrixXd A;
            A.setZero(rows, rows);
            SobolevCurves::SobolevGramMatrix(curves, alpha, beta, A, epsilon);
            return A;
        }

        Eigen::VectorXd DirectSolve(Eigen::VectorXd &b) const {
            Eigen::MatrixXd A = GetFullMatrix();
            return A.partialPivLu().solve(b);
        }

        int NumVertices() const {
            return curves->NumVertices();
        }

        int NumRows() const {
            return curves->NumVertices();
        }

        MultigridMode GetMode() const {
            return MultigridMode::MatrixOnly;
        }

        NullSpaceProjector* GetConstraintProjector() const {
            return 0;
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

        Eigen::VectorXd DirectSolve(Eigen::VectorXd &b) const {
            Eigen::MatrixXd A = GetFullMatrix();
            return A.partialPivLu().solve(b);
        }

        int NumVertices() const {
            return curves->NumVertices();
        }

        int NumRows() const {
            return curves->NumVertices();
        }

        MultigridMode GetMode() const {
            return MultigridMode::MatrixOnly;
        }

        NullSpaceProjector* GetConstraintProjector() const {
            return 0;
        }
    };
}