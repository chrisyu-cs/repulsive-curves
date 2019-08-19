#pragma once

#include "vertex_body.h"
#include "tpe_bvh.h"

namespace LWS {

    struct ClusterPair {
        BVHNode3D* cluster1;
        BVHNode3D* cluster2;
    };

    class BlockClusterTree {
        public:
        BlockClusterTree(PolyCurveGroup* cg, BVHNode3D* tree, double sepCoeff, double a, double b);
        // Loop over all currently inadmissible cluster pairs
        // and subdivide them to their children.
        void splitInadmissibleNodes();
        static bool isPairAdmissible(ClusterPair pair, double coeff);
        static bool isPairSmallEnough(ClusterPair pair);

        void PrintData();

        void MultiplyVector(std::vector<double> &v, std::vector<double> &b);

        void AfFullProduct_verts(ClusterPair pair, std::vector<double> &v, std::vector<double> &result);
        void AfApproxProduct_verts(ClusterPair pair, std::vector<double> &v, std::vector<double> &result);

        void AfFullProduct_hat(ClusterPair pair, std::vector<Vector3> &v_hat, std::vector<Vector3> &result);
        void AfApproxProduct_hat(ClusterPair pair, std::vector<Vector3> &v_hat, std::vector<Vector3> &result);

        private:
        double alpha, beta, separationCoeff;
        PolyCurveGroup* curves;
        BVHNode3D* tree_root;
        std::vector<ClusterPair> admissiblePairs;
        std::vector<ClusterPair> unresolvedPairs;
        std::vector<ClusterPair> imadmissiblePairs;
    };
}