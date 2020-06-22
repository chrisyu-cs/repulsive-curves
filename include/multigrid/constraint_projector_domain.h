#pragma once

#include "libgmultigrid/multigrid_domain.h"
#include "spatial/tpe_bvh.h"
#include "libgmultigrid/matrix_free.h"
#include "product/dense_matrix.h"
#include "constraint_projector_operator.h"

namespace LWS {

    template<typename Constraint>
    class ConstraintProjectorDomain : public MultigridDomain<BlockClusterTree, MatrixProjectorOperator> {
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

        virtual MultigridDomain<BlockClusterTree, MatrixProjectorOperator>* Coarsen(MatrixProjectorOperator* prolongOp) const {
            PolyCurveNetwork* coarsened = curves->Coarsen(prolongOp);
            ConstraintProjectorDomain<Constraint>* coarseDomain = new ConstraintProjectorDomain<Constraint>(coarsened, alpha, beta, sepCoeff, epsilon);
            prolongOp->lowerP = coarseDomain->GetConstraintProjector();
            prolongOp->upperP = GetConstraintProjector();
            coarseDomain->isTopLevel = false;
            return coarseDomain;
        }

        virtual BlockClusterTree* GetMultiplier() const {
            return tree;
        }

        virtual Eigen::MatrixXd GetFullMatrix() const {
            Eigen::MatrixXd A;
            SobolevCurves::Sobolev3XWithConstraints<Constraint>(curves, alpha, beta, A);
            return A;
        }

        virtual Eigen::VectorXd DirectSolve(Eigen::VectorXd &b) const {
            int fullRows = constraint.SaddleNumRows();
            int constraintRows = constraint.NumConstraintRows();

            Eigen::VectorXd b_aug(fullRows);
            b_aug.block(0, 0, b.rows(), 1) = b;
            b_aug.block(b.rows(), 0, constraintRows, 1).setZero();

            Eigen::MatrixXd A = GetFullMatrix();
            b_aug = A.partialPivLu().solve(b_aug);
            return b_aug.block(0, 0, b.rows(), 1);
        }

        virtual int NumVertices() const {
            return nVerts;
        }

        virtual int NumRows() const {
            return 3 * nVerts;
        }

        virtual MatrixProjectorOperator* MakeNewOperator() const {
            return new MatrixProjectorOperator();
        }

        virtual NullSpaceProjector* GetConstraintProjector() const {
            return curves->constraintProjector;
        }
    };
}