#pragma once

#include "spatial/tpe_bvh.h"
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
        PolyCurveNetwork* curves;
        BVHNode3D* bvh;
        BlockClusterTree* tree;
        double alpha, beta;
        double sepCoeff;
        int nVerts;
        bool isTopLevel;

        EdgeLengthSaddleDomain(PolyCurveNetwork* c, double a, double b, double sep) {
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
            PolyCurveNetwork* coarsened = curves->Coarsen(prolongOp, true);
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
            SobolevCurves::Sobolev3XWithConstraints<EdgeLengthConstraint>(curves, alpha, beta, A);
            // SobolevCurves::SobolevPlusBarycenter3X(curves, alpha, beta, A);
            // SobolevCurves::AddEdgeLengthConstraints(curves, A, 3 * nVerts + 3);
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

    class EdgeLengthDenseProjectorDomain : public MultigridDomain<EdgeLengthDenseProjectorDomain, DenseMatrixMult> {
        public:
        PolyCurveNetwork* curves;
        double alpha, beta;
        double sepCoeff;
        int nVerts;
        double epsilon;
        DenseMatrixMult* M;

        EdgeLengthDenseProjectorDomain(PolyCurveNetwork* c, double a, double b, double sep, double diagEps = 0) {
            curves = c;
            alpha = a;
            beta = b;
            sepCoeff = sep;
            nVerts = curves->NumVertices();
            epsilon = diagEps;

            EdgeLengthConstraint constraint(curves);
            curves->AddConstraintProjector(constraint);

            Eigen::MatrixXd A = GetFullMatrix();
            M = new DenseMatrixMult(A);
        }

        ~EdgeLengthDenseProjectorDomain() {
        }

        MultigridDomain<EdgeLengthDenseProjectorDomain, DenseMatrixMult>* Coarsen(MultigridOperator &prolongOp) const {
            PolyCurveNetwork* coarsened = curves->Coarsen(prolongOp);
            EdgeLengthDenseProjectorDomain* coarseDomain = new EdgeLengthDenseProjectorDomain(coarsened, alpha, beta, sepCoeff, epsilon);
            prolongOp.lowerP = coarseDomain->GetConstraintProjector();
            prolongOp.upperP = GetConstraintProjector();
            return coarseDomain;
        }

        DenseMatrixMult* GetMultiplier() const {
            return M;
        }

        Eigen::MatrixXd GetFullMatrix() const {
            int rows = nVerts + 1;
            Eigen::MatrixXd A;
            SobolevCurves::Sobolev3XWithConstraints<EdgeLengthConstraint>(curves, alpha, beta, A, epsilon);
            // SobolevCurves::SobolevPlusBarycenter3X(curves, alpha, beta, A, epsilon);
            // SobolevCurves::AddEdgeLengthConstraints(curves, A, 3 * nVerts + 3);
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

    template<typename Constraint>
    class ConstraintProjectorDomain : public MultigridDomain<ConstraintProjectorDomain<Constraint>, BlockClusterTree> {
        public:
        PolyCurveNetwork* curves;
        BVHNode3D* bvh;
        BlockClusterTree* tree;
        double alpha, beta;
        double sepCoeff;
        int nVerts;
        double epsilon;
        Constraint constraint;

        ConstraintProjectorDomain<Constraint>(PolyCurveNetwork* c, double a, double b, double sep, double diagEps = 0)
        : constraint(c) {
            curves = c;
            alpha = a;
            beta = b;
            sepCoeff = sep;
            nVerts = curves->NumVertices();
            epsilon = diagEps;

            bvh = CreateEdgeBVHFromCurve(curves);
            tree = new BlockClusterTree(curves, bvh, sepCoeff, alpha, beta, epsilon);
            tree->SetBlockTreeMode(BlockTreeMode::Matrix3AndProjector);

            curves->AddConstraintProjector(constraint);
            // std::cout << "Made level with " << nVerts << std::endl;
        }

        ~ConstraintProjectorDomain<Constraint>() {
            delete tree;
            delete bvh;
        }

        MultigridDomain<ConstraintProjectorDomain<Constraint>, BlockClusterTree>* Coarsen(MultigridOperator &prolongOp) const {
            PolyCurveNetwork* coarsened = curves->Coarsen(prolongOp);
            ConstraintProjectorDomain<Constraint>* coarseDomain = new ConstraintProjectorDomain<Constraint>(coarsened, alpha, beta, sepCoeff, epsilon);
            prolongOp.lowerP = coarseDomain->GetConstraintProjector();
            prolongOp.upperP = GetConstraintProjector();
            return coarseDomain;
        }

        BlockClusterTree* GetMultiplier() const {
            return tree;
        }

        Eigen::MatrixXd GetFullMatrix() const {
            Eigen::MatrixXd A;
            SobolevCurves::Sobolev3XWithConstraints<Constraint>(curves, alpha, beta, A);
            return A;
        }

        Eigen::VectorXd DirectSolve(Eigen::VectorXd &b) const {
            int fullRows = constraint.SaddleNumRows();
            int constraintRows = constraint.NumConstraintRows();

            Eigen::VectorXd b_aug(fullRows);
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

    class EdgeLengthNullProjectorDomain : public MultigridDomain<EdgeLengthNullProjectorDomain, BlockClusterTree> {
        public:
        PolyCurveNetwork* curves;
        BVHNode3D* bvh;
        BlockClusterTree* tree;
        double alpha, beta;
        double sepCoeff;
        int nVerts;
        double epsilon;

        EdgeLengthNullProjectorDomain(PolyCurveNetwork* c, double a, double b, double sep, double diagEps = 0) {
            curves = c;
            alpha = a;
            beta = b;
            sepCoeff = sep;
            nVerts = curves->NumVertices();
            epsilon = diagEps;

            bvh = CreateEdgeBVHFromCurve(curves);
            tree = new BlockClusterTree(curves, bvh, sepCoeff, alpha, beta, epsilon);
            tree->SetBlockTreeMode(BlockTreeMode::Matrix3AndProjector);

            EdgeLengthConstraint constraint(curves);
            curves->AddConstraintProjector(constraint);

            // std::cout << "Made level with " << nVerts << std::endl;
        }

        ~EdgeLengthNullProjectorDomain() {
            delete tree;
            delete bvh;
        }

        MultigridDomain<EdgeLengthNullProjectorDomain, BlockClusterTree>* Coarsen(MultigridOperator &prolongOp) const {
            PolyCurveNetwork* coarsened = curves->Coarsen(prolongOp);
            EdgeLengthNullProjectorDomain* coarseDomain = new EdgeLengthNullProjectorDomain(coarsened, alpha, beta, sepCoeff, epsilon);
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
            SobolevCurves::Sobolev3XWithConstraints<EdgeLengthConstraint>(curves, alpha, beta, A);
            // SobolevCurves::SobolevPlusBarycenter3X(curves, alpha, beta, A, epsilon);
            // SobolevCurves::AddEdgeLengthConstraints(curves, A, 3 * nVerts + 3);
            return A;
        }

        Eigen::VectorXd DirectSolve(Eigen::VectorXd &b) const {
            EdgeLengthConstraint constraint(curves);
            int fullRows = constraint.SaddleNumRows();
            int constraintRows = constraint.NumConstraintRows();

            Eigen::VectorXd b_aug(fullRows);
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

    class DenseEdgeNullProjectorDomain : public MultigridDomain<DenseEdgeNullProjectorDomain, DenseMatrixMult> {
        public:
        PolyCurveNetwork* curves;
        double alpha, beta;
        double sepCoeff;
        int nVerts;
        double epsilon;
        Eigen::MatrixXd PAP;
        DenseMatrixMult* mult;

        DenseEdgeNullProjectorDomain(PolyCurveNetwork* c, double a, double b, double sep, double diagEps = 0) {
            curves = c;
            alpha = a;
            beta = b;
            sepCoeff = sep;
            nVerts = curves->NumVertices();
            epsilon = diagEps;

            EdgeLengthConstraint constraint(curves);
            curves->AddConstraintProjector(constraint);

            Eigen::MatrixXd P = curves->constraintProjector->ProjectorMatrix();
            PAP.setZero(nVerts * 3, nVerts * 3);
            SobolevCurves::SobolevGramMatrix3X(curves, alpha, beta, PAP, epsilon);
            PAP = P * PAP * P;

            // std::cout << "Created domain with " << nVerts << std::endl;

            mult = new DenseMatrixMult(PAP);
        }

        ~DenseEdgeNullProjectorDomain() {
            delete mult;
        }

        MultigridDomain<DenseEdgeNullProjectorDomain, DenseMatrixMult>* Coarsen(MultigridOperator &prolongOp) const {
            PolyCurveNetwork* coarsened = curves->Coarsen(prolongOp);
            DenseEdgeNullProjectorDomain* coarseDomain = new DenseEdgeNullProjectorDomain(coarsened, alpha, beta, sepCoeff, epsilon);
            prolongOp.lowerP = coarseDomain->GetConstraintProjector();
            prolongOp.upperP = GetConstraintProjector();
            return coarseDomain;
        }

        DenseMatrixMult* GetMultiplier() const {
            return mult;
        }

        Eigen::MatrixXd GetFullMatrix() const {
            int rows = nVerts + 1;
            Eigen::MatrixXd A;
            SobolevCurves::Sobolev3XWithConstraints<EdgeLengthConstraint>(curves, alpha, beta, A);
            // SobolevCurves::SobolevPlusBarycenter3X(curves, alpha, beta, A, epsilon);
            // SobolevCurves::AddEdgeLengthConstraints(curves, A, 3 * nVerts + 3);
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
        PolyCurveNetwork* curves;
        BVHNode3D* bvh;
        BlockClusterTree* tree;
        double alpha, beta;
        double sepCoeff;
        int nVerts;

        Barycenter3Domain(PolyCurveNetwork* c, double a, double b, double sep) {
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
            PolyCurveNetwork* coarsened = curves->Coarsen(prolongOp);
            return new Barycenter3Domain(coarsened, alpha, beta, sepCoeff);
        }

        BlockClusterTree* GetMultiplier() const {
            return tree;
        }

        Eigen::MatrixXd GetFullMatrix() const {
            Eigen::MatrixXd A;
            
            SobolevCurves::Sobolev3XWithConstraints<BarycenterConstraint3X>(curves, alpha, beta, A);
            // SobolevCurves::SobolevPlusBarycenter3X(curves, alpha, beta, A);
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
        PolyCurveNetwork* curves;
        BVHNode3D* bvh;
        BlockClusterTree* tree;
        double alpha, beta;
        double sepCoeff;
        int nVerts;

        PolyCurveNullProjectorDomain(PolyCurveNetwork* c, double a, double b, double sep) {
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
            PolyCurveNetwork* coarsened = curves->Coarsen(prolongOp);
            PolyCurveNullProjectorDomain* coarseDomain = new PolyCurveNullProjectorDomain(coarsened, alpha, beta, sepCoeff);
            prolongOp.lowerP = coarseDomain->GetConstraintProjector();
            prolongOp.upperP = GetConstraintProjector();
            return coarseDomain;
        }

        BlockClusterTree* GetMultiplier() const {
            return tree;
        }

        Eigen::MatrixXd GetFullMatrix() const {
            int nC = curves->NumComponents();
            int rows = nVerts + nC;
            Eigen::MatrixXd A;
            A.setZero(rows, rows);
            SobolevCurves::SobolevGramMatrix(curves, alpha, beta, A);

            BarycenterConstraint constraint(curves);
            Eigen::SparseMatrix<double> B_s;
            constraint.FillConstraintMatrix(B_s);
            Eigen::MatrixXd B = B_s.toDense();

            A.block(nVerts, 0, nC, nVerts) = B;
            A.block(0, nVerts, nVerts, nC) = B.transpose();

            return A;
        }

        Eigen::VectorXd DirectSolve(Eigen::VectorXd &b) const {
            int nC = curves->NumComponents();
            Eigen::VectorXd b_aug(nVerts + nC);
            b_aug.block(0, 0, nVerts, 1) = b;
            for (int c = 0; c < nC; c++ ){
                b_aug(nVerts + c) = 0;
            }

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
        PolyCurveNetwork* curves;
        BVHNode3D* bvh;
        BlockClusterTree* tree;
        double alpha, beta;
        double sepCoeff;
        int nVerts;

        PolyCurveSaddleDomain(PolyCurveNetwork* c, double a, double b, double sep) {
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
            PolyCurveNetwork* coarsened = curves->Coarsen(prolongOp);
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
        PolyCurveNetwork* curves;
        BVHNode3D* bvh;
        BlockClusterTree* tree;
        double alpha, beta;
        double epsilon;
        double sepCoeff;
        int nVerts;

        PolyCurveGramDomain(PolyCurveNetwork* c, double a, double b, double sep, double e) {
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
            PolyCurveNetwork* coarsened = curves->Coarsen(prolongOp);
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
        PolyCurveNetwork* curves;
        DenseMatrixMult* multiplier;
        double alpha, beta;
        double epsilon;
        
        PolyCurveDenseDomain(PolyCurveNetwork* c, double a, double b, double e) {
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
            PolyCurveNetwork* coarsened = curves->Coarsen(prolongOp);
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