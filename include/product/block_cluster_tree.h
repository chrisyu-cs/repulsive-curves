#pragma once

#include "../spatial/vertex_body.h"
#include "../spatial/tpe_bvh.h"
#include "vector_multiplier.h"

#include "Eigen/Dense"

namespace LWS {

    struct ClusterPair {
        BVHNode3D* cluster1;
        BVHNode3D* cluster2;
    };

    enum class BlockTreeMode {
        MatrixOnly,
        Barycenter,
        EdgeConstraint
    };

    class BlockClusterTree : public VectorMultiplier<BlockClusterTree> {
        public:
        BlockClusterTree(PolyCurveGroup* cg, BVHNode3D* tree, double sepCoeff, double a, double b);
        // Loop over all currently inadmissible cluster pairs
        // and subdivide them to their children.
        void splitInadmissibleNodes();
        static bool isPairAdmissible(ClusterPair pair, double coeff);
        static bool isPairSmallEnough(ClusterPair pair);

        void PrintData();

        void Multiply(std::vector<double> &v, std::vector<double> &b) const;

        // Multiplies A * v and stores it in b.
        void MultiplyVector(std::vector<double> &v, std::vector<double> &b) const;
        // Multiplies A' * v, where the last row and column of A store a barycenter constraint row.
        // Stores the result in b.
        void MultiplyWithBarycenter(std::vector<double> &v, std::vector<double> &b) const;

        void AfFullProduct_hat(ClusterPair pair, std::vector<Vector3> &v_hat, std::vector<Vector3> &result) const;
        void AfApproxProduct_hat(ClusterPair pair, std::vector<Vector3> &v_hat, std::vector<Vector3> &result) const;

        void CompareBlocks();

        Eigen::MatrixXd AfFullBlock(ClusterPair pair);
        Eigen::MatrixXd AfApproxBlock(ClusterPair pair);

        void SetBlockTreeMode(BlockTreeMode m);

        private:
        BlockTreeMode mode;
        double alpha, beta, separationCoeff;
        PolyCurveGroup* curves;
        BVHNode3D* tree_root;
        std::vector<ClusterPair> admissiblePairs;
        std::vector<ClusterPair> unresolvedPairs;
        std::vector<ClusterPair> inadmissiblePairs;
    };
}