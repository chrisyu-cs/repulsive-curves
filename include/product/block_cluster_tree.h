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
        MatrixAndConstraints,
        Matrix3Only,
        Matrix3AndProjector,
        Matrix3AndConstraints
    };

    class BlockClusterTree : public VectorMultiplier<BlockClusterTree> {
        public:
        static long illSepTime;
        static long wellSepTime;
        static long traversalTime;

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

        // Multiplies the inadmissible clusters for A * v, storing it in b.
        void MultiplyInadmissible(Eigen::MatrixXd &v_hat, Eigen::MatrixXd &b_hat) const;
        // Multiplies the admissible clusters for A * v, storing it in b.
        void MultiplyAdmissible(Eigen::MatrixXd &v, Eigen::MatrixXd &b) const;

        // Multiplies A * v, where v holds a vector3 at each vertex in a flattened column,
        //  and stores it in b.
        template<typename V3, typename Dest>
        void MultiplyVector3(V3 &v, Dest &b) const;

        // Multiplies A' * v, where A' is a matrix A augmented
        // by a constraint block.
        template<typename V, typename Dest>
        void MultiplyWithConstraints(V &v, Dest &b) const;

        // Combination of the above two.
        template<typename V3, typename Dest>
        void MultiplyWithConstraints3(V3 &v, Dest &b) const;

        template<typename T>
        void SetConstraints(GradientConstraints<T> &constraints) {
            constraints.FillConstraintMatrix(B);
            constraintsSet = true;
        }

        void AfFullProduct(ClusterPair pair, Eigen::MatrixXd &v_hat, Eigen::MatrixXd &result) const;
        void AfApproxProduct(ClusterPair pair, Eigen::MatrixXd &v_hat, Eigen::MatrixXd &result) const;

        void SetVIs(BVHNode3D* root, Eigen::VectorXd &v_hat);

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
            curves->constraintProjector->ProjectToNullspace(v, tmp);

            Eigen::VectorXd tmp2(v.rows());
            tmp2.setZero();
            MultiplyVector(tmp, tmp2);

            curves->constraintProjector->ProjectToNullspace(tmp2, b);
        }
        else if (mode == BlockTreeMode::MatrixAndConstraints) {
            MultiplyWithConstraints(v, b);
        }

        else if (mode == BlockTreeMode::Matrix3Only) {
            MultiplyVector3(v, b);
        }

        else if (mode == BlockTreeMode::Matrix3AndProjector) {
            long start = Utils::currentTimeMilliseconds();
            Eigen::VectorXd tmp(v.rows());
            tmp.setZero();
            curves->constraintProjector->ProjectToNullspace(v, tmp);
            long proj1 = Utils::currentTimeMilliseconds();

            Eigen::VectorXd tmp2(v.rows());
            tmp2.setZero();
            MultiplyVector3(tmp, tmp2);
            long proj2 = Utils::currentTimeMilliseconds();

            curves->constraintProjector->ProjectToNullspace(tmp2, b);
            long end = Utils::currentTimeMilliseconds();
        }

        else if (mode == BlockTreeMode::Matrix3AndConstraints) {
            MultiplyWithConstraints3(v, b);
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

        Eigen::MatrixXd b_hat_adm(nVerts, 3);
        b_hat_adm.setZero();

        Eigen::MatrixXd b_hat_inadm(nVerts, 3);
        b_hat_inadm.setZero();

        long illSepStart = Utils::currentTimeMilliseconds();
        MultiplyInadmissible(v_hat, b_hat_inadm);
        long middle = Utils::currentTimeMilliseconds();
        MultiplyAdmissible(v_hat, b_hat_adm);
        long wellSepEnd = Utils::currentTimeMilliseconds();

        b_hat_adm += b_hat_inadm;

        long illTime = (middle - illSepStart);
        long wellTime = (wellSepEnd - middle);

        illSepTime += illTime;
        wellSepTime += wellTime;

        SobolevCurves::ApplyDfTranspose(curves, b_hat_adm, b);
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

    template<typename V3, typename Dest>
    void BlockClusterTree::MultiplyWithConstraints3(V3 &v, Dest &b) const {
        if (!constraintsSet) {
            std::cerr << "Called MultiplyWithConstraints3 without calling SetConstraints beforehand" << std::endl;
            throw 1;
        }
        int B_start = nVerts * 3;
        int nConstraints = B.rows();
        // Do the multiplication with the top-left block
        MultiplyVector3(v, b);
        // Now add the results from the constraint blocks
        // Multiply B * v
        b.block(B_start, 0, nConstraints, 1) += B * v.block(0, 0, B_start, 1);
        // Multiply B^T * phi (lagrange multipliers block)
        b.block(0, 0, B_start, 1) += B.transpose() * v.block(B_start, 0, nConstraints, 1);
    }
}