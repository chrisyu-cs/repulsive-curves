#pragma once

#include "spatial/tpe_bvh.h"
#include "poly_curve.h"
#include "multigrid_operator.h"
#include "product/matrix_free.h"
#include "product/dense_matrix.h"
#include "product/test_matrices.h"
#include "nullspace_projector.h"
#include "flow/gradient_constraint_types.h"

#include <Eigen/Core>

namespace LWS {

    template<typename T, typename Mult>
    class MultigridDomain {
        public:
        typedef Mult MultType;

        MultigridDomain<T, Mult>* Coarsen(MultigridOperator &prolongOp) const {
            return static_cast<T const&>(*this).Coarsen(prolongOp);
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

        ProlongationMode GetMode() const {
            return static_cast<T const&>(*this).GetMode();
        }
        
        NullSpaceProjector* GetConstraintProjector() const {
            return static_cast<T const&>(*this).GetConstraintProjector();
        }
    };

    class EdgeLengthSaddleDomain : public MultigridDomain<EdgeLengthSaddleDomain, BlockClusterTree> {
        public:
        PolyCurveGroup* curves;
        BVHNode3D* bvh;
        BlockClusterTree* tree;
        double alpha, beta;
        double sepCoeff;
        int nVerts;
        bool isTopLevel;

        EdgeLengthSaddleDomain(PolyCurveGroup* c, double a, double b, double sep) {
            curves = c;
            alpha = a;
            beta = b;
            sepCoeff = sep;
            nVerts = curves->NumVertices();

            bvh = CreateEdgeBVHFromCurve(curves);
            tree = new BlockClusterTree(curves, bvh, sepCoeff, alpha, beta);
            tree->SetBlockTreeMode(BlockTreeMode::Matrix3AndConstraints);

            EdgeLengthConstraint constraint(curves);
            tree->SetConstraints(constraint);
            isTopLevel = true;
        }

        ~EdgeLengthSaddleDomain() {
            if (!isTopLevel) {
                delete curves;
            }
            delete tree;
            delete bvh;
        }

        MultigridDomain<EdgeLengthSaddleDomain, BlockClusterTree>* Coarsen(MultigridOperator &prolongOp) const {
            PolyCurveGroup* coarsened = curves->Coarsen(prolongOp);
            curves->EdgeProlongation(prolongOp);
            EdgeLengthSaddleDomain* sub = new EdgeLengthSaddleDomain(coarsened, alpha, beta, sepCoeff);
            sub->isTopLevel = false;
            return sub;
        }

        BlockClusterTree* GetMultiplier() const {
            return tree;
        }

        Eigen::MatrixXd GetFullMatrix() const {
            int rows = nVerts + 1;
            Eigen::MatrixXd A;
            A.setZero(4 * nVerts + 3, 4 * nVerts + 3);
            SobolevCurves::SobolevPlusBarycenter3X(curves, alpha, beta, A);
            SobolevCurves::AddEdgeLengthConstraints(curves, A, 3 * nVerts + 3);
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
            return 4 * nVerts + 3;
        }

        ProlongationMode GetMode() const {
            return ProlongationMode::Matrix3AndEdgeConstraints;
        }

        NullSpaceProjector* GetConstraintProjector() const {
            return 0;
        }
    };

    class EdgeLengthNullProjectorDomain : public MultigridDomain<EdgeLengthNullProjectorDomain, BlockClusterTree> {
        public:
        PolyCurveGroup* curves;
        BVHNode3D* bvh;
        BlockClusterTree* tree;
        double alpha, beta;
        double sepCoeff;
        int nVerts;

        EdgeLengthNullProjectorDomain(PolyCurveGroup* c, double a, double b, double sep) {
            curves = c;
            alpha = a;
            beta = b;
            sepCoeff = sep;
            nVerts = curves->NumVertices();

            bvh = CreateEdgeBVHFromCurve(curves);
            tree = new BlockClusterTree(curves, bvh, sepCoeff, alpha, beta);
            tree->SetBlockTreeMode(BlockTreeMode::Matrix3AndProjector);

            EdgeLengthConstraint constraint(curves);
            curves->AddConstraintProjector(constraint);
        }

        ~EdgeLengthNullProjectorDomain() {
            delete tree;
            delete bvh;
        }

        MultigridDomain<EdgeLengthNullProjectorDomain, BlockClusterTree>* Coarsen(MultigridOperator &prolongOp) const {
            PolyCurveGroup* coarsened = curves->Coarsen(prolongOp);
            EdgeLengthNullProjectorDomain* coarseDomain = new EdgeLengthNullProjectorDomain(coarsened, alpha, beta, sepCoeff);
            prolongOp.lowerP = coarseDomain->GetConstraintProjector();
            prolongOp.upperP = GetConstraintProjector();
            return coarseDomain;
        }

        BlockClusterTree* GetMultiplier() const {
            return tree;
        }

        Eigen::MatrixXd GetFullMatrix() const {
            int rows = nVerts + 1;
            Eigen::MatrixXd A;
            A.setZero(4 * nVerts + 3, 4 * nVerts + 3);
            SobolevCurves::SobolevPlusBarycenter3X(curves, alpha, beta, A);
            SobolevCurves::AddEdgeLengthConstraints(curves, A, 3 * nVerts + 3);
            return A;
        }

        Eigen::VectorXd DirectSolve(Eigen::VectorXd &b) const {
            Eigen::VectorXd b_aug(4 * nVerts + 3);
            int constraintRows = nVerts + 3;
            b_aug.block(0, 0, b.rows(), 1) = b;
            b_aug.block(b.rows(), 0, constraintRows, 1).setZero();

            Eigen::MatrixXd A = GetFullMatrix();
            b_aug = A.partialPivLu().solve(b_aug);
            return b_aug.block(0, 0, b.rows(), 1);
        }

        int NumVertices() const {
            return nVerts;
        }

        int NumRows() const {
            return 3 * nVerts;
        }

        ProlongationMode GetMode() const {
            return ProlongationMode::Matrix3AndProjector;
        }

        NullSpaceProjector* GetConstraintProjector() const {
            return curves->constraintProjector;
        }

    };

    class Barycenter3Domain : public MultigridDomain<Barycenter3Domain, BlockClusterTree> {
        public:
        PolyCurveGroup* curves;
        BVHNode3D* bvh;
        BlockClusterTree* tree;
        double alpha, beta;
        double sepCoeff;
        int nVerts;

        Barycenter3Domain(PolyCurveGroup* c, double a, double b, double sep) {
            curves = c;
            alpha = a;
            beta = b;
            sepCoeff = sep;
            nVerts = curves->NumVertices();

            bvh = CreateEdgeBVHFromCurve(curves);
            tree = new BlockClusterTree(curves, bvh, sepCoeff, alpha, beta);
            tree->SetBlockTreeMode(BlockTreeMode::Matrix3AndConstraints);

            BarycenterConstraint3X constraint(curves);
            tree->SetConstraints(constraint);
        }

        ~Barycenter3Domain() {
            delete tree;
            delete bvh;
        }

        MultigridDomain<Barycenter3Domain, BlockClusterTree>* Coarsen(MultigridOperator &prolongOp) const {
            PolyCurveGroup* coarsened = curves->Coarsen(prolongOp);
            return new Barycenter3Domain(coarsened, alpha, beta, sepCoeff);
        }

        BlockClusterTree* GetMultiplier() const {
            return tree;
        }

        Eigen::MatrixXd GetFullMatrix() const {
            Eigen::MatrixXd A;
            A.setZero(3 * nVerts + 3, 3 * nVerts + 3);
            SobolevCurves::SobolevPlusBarycenter3X(curves, alpha, beta, A);
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
            return 3 * nVerts + 3;
        }

        ProlongationMode GetMode() const {
            return ProlongationMode::Matrix3AndBarycenter;
        }

        NullSpaceProjector* GetConstraintProjector() const {
            return 0;
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

            BarycenterConstraint constraint(curves);
            curves->AddConstraintProjector(constraint);
        }

        ~PolyCurveNullProjectorDomain() {
            delete tree;
            delete bvh;
        }

        MultigridDomain<PolyCurveNullProjectorDomain, BlockClusterTree>* Coarsen(MultigridOperator &prolongOp) const {
            PolyCurveGroup* coarsened = curves->Coarsen(prolongOp);
            PolyCurveNullProjectorDomain* coarseDomain = new PolyCurveNullProjectorDomain(coarsened, alpha, beta, sepCoeff);
            prolongOp.lowerP = coarseDomain->GetConstraintProjector();
            prolongOp.upperP = GetConstraintProjector();
            return coarseDomain;
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

        ProlongationMode GetMode() const {
            return ProlongationMode::MatrixAndProjector;
        }

        NullSpaceProjector* GetConstraintProjector() const {
            return curves->constraintProjector;
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

        MultigridDomain<PolyCurveSaddleDomain, BlockClusterTree>* Coarsen(MultigridOperator &prolongOp) const {
            PolyCurveGroup* coarsened = curves->Coarsen(prolongOp);
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

        ProlongationMode GetMode() const {
            return ProlongationMode::Barycenter;
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

        MultigridDomain<PolyCurveGramDomain, BlockClusterTree>* Coarsen(MultigridOperator &prolongOp) const {
            PolyCurveGroup* coarsened = curves->Coarsen(prolongOp);
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

        ProlongationMode GetMode() const {
            return ProlongationMode::MatrixOnly;
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

        MultigridDomain<PolyCurveDenseDomain, DenseMatrixMult>* Coarsen(MultigridOperator &prolongOp) const {
            PolyCurveGroup* coarsened = curves->Coarsen(prolongOp);
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

        ProlongationMode GetMode() const {
            return ProlongationMode::MatrixOnly;
        }

        NullSpaceProjector* GetConstraintProjector() const {
            return 0;
        }
    };
}