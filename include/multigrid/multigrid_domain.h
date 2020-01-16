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
        virtual ~MultigridDomain() {}

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
        bool isTopLevel;

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
            isTopLevel = true;
        }

        virtual ~ConstraintProjectorDomain<Constraint>() {
            if (!isTopLevel) {
                delete curves;
            }
            delete tree;
            delete bvh;
        }

        MultigridDomain<ConstraintProjectorDomain<Constraint>, BlockClusterTree>* Coarsen(MultigridOperator &prolongOp) const {
            PolyCurveNetwork* coarsened = curves->Coarsen(prolongOp);
            ConstraintProjectorDomain<Constraint>* coarseDomain = new ConstraintProjectorDomain<Constraint>(coarsened, alpha, beta, sepCoeff, epsilon);
            prolongOp.lowerP = coarseDomain->GetConstraintProjector();
            prolongOp.upperP = GetConstraintProjector();
            coarseDomain->isTopLevel = false;
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

        virtual ~EdgeLengthSaddleDomain() {
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

        virtual ~EdgeLengthNullProjectorDomain() {
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
}