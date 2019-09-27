#pragma once

#include "../spatial/vertex_body.h"
#include "../spatial/tpe_bvh.h"
#include "vector_multiplier.h"
#include "sobo_slobo.h"

#include "Eigen/Dense"

namespace LWS {

    struct ClusterPair {
        BVHNode3D* cluster1;
        BVHNode3D* cluster2;
    };

    enum class BlockTreeMode {
        MatrixOnly,
        MatrixAndProjector,
        Barycenter,
        EdgeConstraint
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

        // Multiplies A' * v, where the last row and column of A store a barycenter constraint row.
        // Stores the result in b.
        template<typename V, typename Dest>
        void MultiplyWithBarycenter(V &v, Dest &b) const;

        void AfFullProduct_hat(ClusterPair pair, Eigen::MatrixXd &v_hat, Eigen::MatrixXd &result) const;
        void AfApproxProduct_hat(ClusterPair pair, Eigen::MatrixXd &v_hat, Eigen::MatrixXd &result) const;

        void CompareBlocks();

        Eigen::MatrixXd AfFullBlock(ClusterPair pair);
        Eigen::MatrixXd AfApproxBlock(ClusterPair pair);

        void SetBlockTreeMode(BlockTreeMode m);

        private:
        BlockTreeMode mode;
        double alpha, beta, separationCoeff;
        double epsilon;
        PolyCurveGroup* curves;
        BVHNode3D* tree_root;
        std::vector<ClusterPair> admissiblePairs;
        std::vector<ClusterPair> unresolvedPairs;
        std::vector<ClusterPair> inadmissiblePairs;
    };

    template<typename V, typename Dest>
    void BlockClusterTree::Multiply(V &v, Dest &b) const {
        if (mode == BlockTreeMode::MatrixOnly) {
            MultiplyVector(v, b);
        }
        else if (mode == BlockTreeMode::MatrixAndProjector) {
            Eigen::VectorXd tmp(v.rows());
            tmp.setZero();
            curves->constraints->Multiply(v, tmp);

            Eigen::VectorXd tmp2(v.rows());
            tmp2.setZero();
            MultiplyVector(tmp, tmp2);

            curves->constraints->Multiply(tmp2, b);
        }
        else if (mode == BlockTreeMode::Barycenter) {
            MultiplyWithBarycenter(v, b);
        }
        else if (mode == BlockTreeMode::EdgeConstraint) {
            // TODO
        }

        int nVerts = curves->NumVertices();
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

    template<typename V, typename Dest>
    void BlockClusterTree::MultiplyWithBarycenter(V &v, Dest &b) const {
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

        double totalLength = curves->TotalLength();

        // Multiply the last row and column of the matrix
        for (int i = 0; i < nVerts; i++) {
            double weight = curves->GetCurvePoint(i).DualLength() / totalLength;
            // Add last column sum
            b(i) += v(nVerts) * weight;
            // Add last row sum
            b(nVerts) += v(i) * weight;
        }
    }
}