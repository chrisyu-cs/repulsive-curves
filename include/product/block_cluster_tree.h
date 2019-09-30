#pragma once

#include "../spatial/vertex_body.h"
#include "../spatial/tpe_bvh.h"
#include "vector_multiplier.h"
#include "sobo_slobo.h"
#include "flow/gradient_constraints.h"

#include "Eigen/Dense"

namespace LWS {

    struct ClusterPair {
        BVHNode3D* cluster1;
        BVHNode3D* cluster2;
    };

    enum class BlockTreeMode {
        MatrixOnly,
        MatrixAndProjector,
        MatrixAndConstraints
    };

    class BlockClusterTree : public VectorMultiplier<BlockClusterTree> {
        public:
        BlockClusterTree(PolyCurveGroup* cg, BVHNode3D* tree, double sepCoeff, double a, double b, double e = 0.0);
        // Loop over all currently inadmissible cluster pairs
        // and subdivide them to their children.
        void splitInadmissibleNodes();
        static bool isPairAdmissible(ClusterPair pair, double coeff);
        static bool isPairSmallEnough(ClusterPair pair);

        void PrintData();

        // Multiplies v and stores in b. Dispatches to the specific multiplication case below.
        template<typename V, typename Dest>
        void Multiply(V &v, Dest &b) const;

        // Multiplies A * v and stores it in b.
        template<typename V, typename Dest>
        void MultiplyVector(V &v, Dest &b) const;

        // Multiplies A * v, where v holds a vector3 at each vertex in a flattened column,
        //  and stores it in b.
        template<typename V3, typename Dest>
        void MultiplyVector3(V3 &v, Dest &b) const;

        // Multiplies A' * v, where A' is a matrix A augmented
        // by a constraint block.
        template<typename V, typename Dest>
        void MultiplyWithConstraints(V &v, Dest &b) const;

        template<typename T>
        void SetConstraints(GradientConstraints<T> &constraints) {
            constraints.FillConstraintMatrix(B);
            constraintsSet = true;
        }

        void AfFullProduct_hat(ClusterPair pair, Eigen::MatrixXd &v_hat, Eigen::MatrixXd &result) const;
        void AfApproxProduct_hat(ClusterPair pair, Eigen::MatrixXd &v_hat, Eigen::MatrixXd &result) const;

        void CompareBlocks();

        Eigen::MatrixXd AfFullBlock(ClusterPair pair);
        Eigen::MatrixXd AfApproxBlock(ClusterPair pair);

        void SetBlockTreeMode(BlockTreeMode m);

        private:
        BlockTreeMode mode;
        int nVerts;
        double alpha, beta, separationCoeff;
        double epsilon;
        PolyCurveGroup* curves;
        BVHNode3D* tree_root;
        std::vector<ClusterPair> admissiblePairs;
        std::vector<ClusterPair> unresolvedPairs;
        std::vector<ClusterPair> inadmissiblePairs;
        bool constraintsSet;
        Eigen::SparseMatrix<double> B;
    };

    template<typename V, typename Dest>
    void BlockClusterTree::Multiply(V &v, Dest &b) const {
        if (mode == BlockTreeMode::MatrixOnly) {
            MultiplyVector(v, b);
        }
        else if (mode == BlockTreeMode::MatrixAndProjector) {
            Eigen::VectorXd tmp(v.rows());
            tmp.setZero();
            curves->constraints->ProjectToNullspace(v, tmp);

            Eigen::VectorXd tmp2(v.rows());
            tmp2.setZero();
            MultiplyVector(tmp, tmp2);

            curves->constraints->ProjectToNullspace(tmp2, b);
        }
        else if (mode == BlockTreeMode::MatrixAndConstraints) {
            MultiplyWithConstraints(v, b);
        }

        for (int i = 0; i < nVerts; i++) {
            b(i) += epsilon * curves->GetCurvePoint(i).DualLength() * v(i);
        }
    }

    template<typename V, typename Dest>
    void BlockClusterTree::MultiplyVector(V &v, Dest &b) const {
        int nVerts = curves->NumVertices();
        Eigen::MatrixXd v_hat(nVerts, 3);
        v_hat.setZero();

        SobolevCurves::ApplyDf(curves, v, v_hat);

        Eigen::MatrixXd b_hat(nVerts, 3);
        b_hat.setZero();

        for (ClusterPair pair : inadmissiblePairs) {
            AfFullProduct_hat(pair, v_hat, b_hat);
        }
        for (ClusterPair pair : admissiblePairs) {
            AfApproxProduct_hat(pair, v_hat, b_hat);
        }

        SobolevCurves::ApplyDfTranspose(curves, b_hat, b);
    }

    template<typename V3, typename Dest>
    void BlockClusterTree::MultiplyVector3(V3 &v, Dest &b) const {
        // Slice the input vector to get every x-coordinate
        Eigen::Map<const Eigen::VectorXd, 0, Eigen::InnerStride<3>> v_x(v.data(), nVerts);
        // Slice the output vector to get x-coordinates
        Eigen::Map<Eigen::VectorXd, 0, Eigen::InnerStride<3>> dest_x(b.data(), nVerts);
        // Multiply the input x-coords into the output x-coords
        MultiplyVector(v_x, dest_x);

        // Same thing for y-coordinates
        Eigen::Map<const Eigen::VectorXd, 0, Eigen::InnerStride<3>> v_y(v.data() + 1, nVerts);
        Eigen::Map<Eigen::VectorXd, 0, Eigen::InnerStride<3>> dest_y(b.data() + 1, nVerts);
        MultiplyVector(v_y, dest_y);

        // Same thing for z-coordinates
        Eigen::Map<const Eigen::VectorXd, 0, Eigen::InnerStride<3>> v_z(v.data() + 2, nVerts);
        Eigen::Map<Eigen::VectorXd, 0, Eigen::InnerStride<3>> dest_z(b.data() + 2, nVerts);
        MultiplyVector(v_z, dest_z);
    }

    template<typename V, typename Dest>
    void BlockClusterTree::MultiplyWithConstraints(V &v, Dest &b) const {
        if (!constraintsSet) {
            std::cerr << "Called MultiplyWithConstraints without calling SetConstraints beforehand" << std::endl;
            throw 1;
        }
        int nConstraints = B.rows();
        // Do the multiplication with the top-left block
        MultiplyVector(v, b);
        // Now add the results from the constraint blocks
        // Multiply B * v
        b.block(nVerts, 0, nConstraints, 1) += B * v.block(0, 0, nVerts, 1);
        // Multiply B^T * phi (lagrange multipliers block)
        b.block(0, 0, nVerts, 1) += B.transpose() * v.block(nVerts, 0, nConstraints, 1);
    }
}