#pragma once

#include "../spatial/vertex_body.h"
#include "../spatial/tpe_bvh.h"
#include "libgmultigrid/vector_multiplier.h"
#include "sobo_slobo.h"
#include "libgmultigrid/domain_constraints.h"
#include "poly_curve_network.h"

#include "Eigen/Dense"
#include <fstream>

namespace LWS
{

    class ClusterPair
    {
    public:
        ClusterPair(BVHNode3D *c1 = NULL, BVHNode3D *c2 = NULL, int d = 0)
            : cluster1(c1), cluster2(c2), depth(d) {}

        BVHNode3D *cluster1;
        BVHNode3D *cluster2;
        int depth;
    };

    enum class BlockTreeMode
    {
        MatrixOnly,
        MatrixAndProjector,
        MatrixAndConstraints,
        Matrix3Only,
        Matrix3AndProjector,
        Matrix3AndConstraints
    };

    class BlockClusterTree : public VectorMultiplier<BlockClusterTree>
    {
    public:
        static long illSepTime;
        static long wellSepTime;
        static long traversalTime;

        BlockClusterTree(PolyCurveNetwork *cg, BVHNode3D *tree, double sepCoeff, double a, double b, double e = 0.0);
        ~BlockClusterTree();
        // Loop over all currently inadmissible cluster pairs
        // and subdivide them to their children.
        void splitInadmissibleNodes(int depth);
        static bool isPairAdmissible(ClusterPair pair, double coeff);
        static bool isPairSmallEnough(ClusterPair pair);

        inline void refreshEdgeWeights()
        {
            tree_root->refreshWeightsVector(curves, BodyType::Edge);
        }

        inline void PremultiplyAfFrac(double s)
        {
            Af_1_frac.setOnes(tree_root->numElements);
            Af_1_frac = MultiplyAfFrac(Af_1_frac, s);
        }

        void PrintData();
        void PrintAdmissibleClusters(std::ofstream &stream);
        void PrintInadmissibleClusters(std::ofstream &stream);

        // Multiplies v and stores in b. Dispatches to the specific multiplication case below.
        template <typename V, typename Dest>
        void Multiply(V &v, Dest &b) const;

        // Multiplies A * v and stores it in b.
        template <typename V, typename Dest>
        void MultiplyVector(V &v, Dest &b) const;

        // Multiplies A * v, where v holds a vector3 at each vertex in a flattened column,
        //  and stores it in b.
        template <typename V3, typename Dest>
        void MultiplyVector3(V3 &v, Dest &b) const;

        // Multiplies L^(s/2) * v, the fractional Laplacian.
        template <typename V, typename Dest>
        void MultiplyByFracLaplacian(V &v, Dest &b, double s) const;

        template <typename V3, typename Dest>
        void MultiplyByFracLaplacian3(V3 &v, Dest &b, double s) const;

        // Multiplies A' * v, where A' is a matrix A augmented
        // by a constraint block.
        template <typename V, typename Dest>
        void MultiplyWithConstraints(V &v, Dest &b) const;

        // Combination of the above two.
        template <typename V3, typename Dest>
        void MultiplyWithConstraints3(V3 &v, Dest &b) const;

        template <typename T>
        void SetConstraints(DomainConstraints<T> &constraints)
        {
            constraints.FillConstraintMatrix(B);
            constraintsSet = true;
        }

        void sum_AIJ_VJ() const;
        void sum_AIJ_VJ_Low() const;
        void sum_AIJ_VJ_Frac(double s) const;

        void AfFullProduct(ClusterPair pair, const Eigen::MatrixXd &v_hat, Eigen::MatrixXd &result) const;
        void AfFullProductLow(ClusterPair pair, const Eigen::VectorXd &v_mid, Eigen::VectorXd &result) const;
        void AfFullProductFrac(ClusterPair pair, const Eigen::VectorXd &v_mid, Eigen::VectorXd &result, double s) const;

        void AfApproxProduct(ClusterPair pair, const Eigen::MatrixXd &v_hat, Eigen::MatrixXd &result) const;
        void AfApproxProductLow(ClusterPair pair, const Eigen::VectorXd &v_mid, Eigen::VectorXd &result) const;

        template <typename V>
        Eigen::VectorXd MultiplyAf(const V &v) const;
        template <typename V>
        Eigen::VectorXd MultiplyAfLow(const V &v) const;
        template <typename V>
        Eigen::VectorXd MultiplyAfFrac(const V &v, double s) const;

        template <typename V>
        void SetVIs(BVHNode3D *node, const V &v_hat) const;
        template <typename V>
        void SetVIsLow(BVHNode3D *node, const V &v_mid) const;

        template <typename Dest>
        void SetBIs(BVHNode3D *node, Dest &b_tilde) const;
        template <typename Dest>
        void SetBIsLow(BVHNode3D *node, Dest &b_tilde) const;
        template <typename Dest>
        void SetBIsFrac(BVHNode3D *node, Dest &b_tilde, double s) const;

        template <typename Dest>
        void PropagateBIs(BVHNode3D *node, double parent_BI, Dest &b_tilde) const;

        void SetBlockTreeMode(BlockTreeMode m) const;

    private:
        // Multiplies the inadmissible clusters for A * v, storing it in b.
        void MultiplyInadmissibleParallel(const Eigen::MatrixXd &v_hat, Eigen::MatrixXd &b_hat) const;
        void MultiplyAdmissibleFast(const Eigen::MatrixXd &v_hat, Eigen::MatrixXd &b_hat) const;

        // Same, but for low-order term.
        void MultiplyAdmissibleLowFast(const Eigen::VectorXd &v_mid, Eigen::VectorXd &b_mid) const;
        void MultiplyInadmissibleLowParallel(const Eigen::VectorXd &v_mid, Eigen::VectorXd &b_mid) const;

        // Fractional laplacian
        void MultiplyAdmissibleFracFast(const Eigen::VectorXd &v_mid, Eigen::VectorXd &b_mid, double s) const;
        void MultiplyInadmissibleFracParallel(const Eigen::VectorXd &v_mid, Eigen::VectorXd &b_mid, double s) const;

        Eigen::VectorXd Af_1, Af_1_low, Af_1_frac;
        mutable BlockTreeMode mode;
        int nVerts;
        double alpha, beta, separationCoeff;
        double epsilon;
        PolyCurveNetwork *curves;
        BVHNode3D *tree_root;
        std::vector<ClusterPair> admissiblePairs;
        std::vector<std::vector<ClusterPair>> admissibleByCluster;

        std::vector<ClusterPair> unresolvedPairs;
        std::vector<ClusterPair> inadmissiblePairs;
        bool constraintsSet;
        Eigen::SparseMatrix<double> B;
    };

    template <typename V>
    Eigen::VectorXd BlockClusterTree::MultiplyAf(const V &v) const
    {
        tree_root->recursivelyZeroMVFields();
        Eigen::VectorXd result(v.rows());
        result.setZero();
        SetVIs(tree_root, v);
        SetBIs(tree_root, result);
        result = tree_root->fullMasses.asDiagonal() * result;
        return result;
    }

    template <typename V>
    Eigen::VectorXd BlockClusterTree::MultiplyAfLow(const V &v) const
    {
        tree_root->recursivelyZeroMVFields();
        Eigen::VectorXd result(v.rows());
        result.setZero();
        SetVIsLow(tree_root, v);
        SetBIsLow(tree_root, result);
        result = tree_root->fullMasses.asDiagonal() * result;
        return result;
    }

    template <typename V>
    Eigen::VectorXd BlockClusterTree::MultiplyAfFrac(const V &v, double s) const
    {
        tree_root->recursivelyZeroMVFields();
        Eigen::VectorXd result(v.rows());
        result.setZero();
        SetVIsLow(tree_root, v); // same procedure as low-order
        SetBIsFrac(tree_root, result, s);
        result = tree_root->fullMasses.asDiagonal() * result;
        return result;
    }

    template <typename V>
    void BlockClusterTree::SetVIs(BVHNode3D *node, const V &v_hat) const
    {
        if (node->IsLeaf())
        {
            int index = node->VertexIndex();
            double w_j = node->totalMass;
            double v_hat_j = v_hat(index);
            node->V_I = w_j * v_hat_j;
        }
        else
        {
            node->V_I = 0;
            // Start at the roots and propagate upward
            for (BVHNode3D *child : node->children)
            {
                SetVIs(child, v_hat);
                node->V_I += child->V_I;
            }
        }
    }

    template <typename V>
    void BlockClusterTree::SetVIsLow(BVHNode3D *node, const V &v_mid) const
    {
        if (node->IsLeaf())
        {
            int index = node->VertexIndex();
            double w_j = node->totalMass;
            double v_mid_j = v_mid(index);
            node->V_I = w_j * v_mid_j;
        }
        else
        {
            node->V_I = 0;
            // Start at the roots and propagate upward
            for (BVHNode3D *child : node->children)
            {
                SetVIsLow(child, v_mid);
                node->V_I += child->V_I;
            }
        }
    }

    template <typename Dest>
    void BlockClusterTree::SetBIs(BVHNode3D *node, Dest &b_tilde) const
    {
        // if (curves->NumVertices() > 1000) sum_AIJ_VJ_Parallel();
        sum_AIJ_VJ();
        node->aIJ_VJ = 0;
        // Now recursively propagate downward
        PropagateBIs(node, 0, b_tilde);
    }

    template <typename Dest>
    void BlockClusterTree::SetBIsLow(BVHNode3D *node, Dest &b_tilde) const
    {
        // if (curves->NumVertices() > 1000) sum_AIJ_VJ_Low_Parallel();
        sum_AIJ_VJ_Low();
        node->aIJ_VJ = 0;
        // Now recursively propagate downward
        PropagateBIs(node, 0, b_tilde);
    }

    template <typename Dest>
    void BlockClusterTree::SetBIsFrac(BVHNode3D *node, Dest &b_tilde, double s) const
    {
        sum_AIJ_VJ_Frac(s);
        node->aIJ_VJ = 0;
        // Now recursively propagate downward
        PropagateBIs(node, 0, b_tilde);
    }

    template <typename Dest>
    void BlockClusterTree::PropagateBIs(BVHNode3D *node, double parent_BI, Dest &b_tilde) const
    {
        node->B_I = parent_BI + node->aIJ_VJ;
        if (node->IsLeaf())
        {
            b_tilde(node->VertexIndex()) = node->B_I;
        }
        else
        {
            for (BVHNode3D *child : node->children)
            {
                PropagateBIs(child, node->B_I, b_tilde);
            }
        }
    }

    template <typename V, typename Dest>
    void BlockClusterTree::Multiply(V &v, Dest &b) const
    {
        if (mode == BlockTreeMode::MatrixOnly)
        {
            MultiplyVector(v, b);
        }
        else if (mode == BlockTreeMode::MatrixAndProjector)
        {
            Eigen::VectorXd tmp(v.rows());
            tmp.setZero();
            curves->constraintProjector->ProjectToNullspace(v, tmp);

            Eigen::VectorXd tmp2(v.rows());
            tmp2.setZero();
            MultiplyVector(tmp, tmp2);

            curves->constraintProjector->ProjectToNullspace(tmp2, b);
        }
        else if (mode == BlockTreeMode::MatrixAndConstraints)
        {
            MultiplyWithConstraints(v, b);
        }

        else if (mode == BlockTreeMode::Matrix3Only)
        {
            MultiplyVector3(v, b);
        }

        else if (mode == BlockTreeMode::Matrix3AndProjector)
        {
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

        else if (mode == BlockTreeMode::Matrix3AndConstraints)
        {
            MultiplyWithConstraints3(v, b);
        }

        for (int i = 0; i < nVerts; i++)
        {
            b(i) += epsilon * curves->GetVertex(i)->DualLength() * v(i);
        }
    }

    template <typename V, typename Dest>
    void BlockClusterTree::MultiplyVector(V &v, Dest &b) const
    {
        int nEdges = curves->NumEdges();
        Eigen::MatrixXd v_hat(nEdges, 3);
        v_hat.setZero();
        Eigen::VectorXd v_mid(nEdges);
        v_mid.setZero();

        // Set up input and outputs for metric
        SobolevCurves::ApplyDf(curves, v, v_hat);
        SobolevCurves::ApplyMid(curves, v, v_mid);

        Eigen::MatrixXd b_hat_adm(nEdges, 3);
        b_hat_adm.setZero();
        Eigen::MatrixXd b_hat_inadm(nEdges, 3);
        b_hat_inadm.setZero();

        // Set up inputs and outputs for low-order part
        Eigen::VectorXd b_mid_adm(nEdges);
        b_mid_adm.setZero();
        Eigen::VectorXd b_mid_inadm(nEdges);
        b_mid_inadm.setZero();

        // Multiply inadmissible blocks
        MultiplyInadmissibleParallel(v_hat, b_hat_inadm);
        MultiplyInadmissibleLowParallel(v_mid, b_mid_inadm);
        // Multiply admissible blocks
        MultiplyAdmissibleFast(v_hat, b_hat_adm);
        MultiplyAdmissibleLowFast(v_mid, b_mid_adm);

        b_hat_adm += b_hat_inadm;
        b_mid_adm += b_mid_inadm;

        SobolevCurves::ApplyDfTranspose(curves, b_hat_adm, b);
        SobolevCurves::ApplyMidTranspose(curves, b_mid_adm, b);
    }

    // Multiplies L^s * v, the fractional Laplacian.
    template <typename V, typename Dest>
    void BlockClusterTree::MultiplyByFracLaplacian(V &v, Dest &b, double s) const
    {
        int nEdges = curves->NumEdges();
        Eigen::VectorXd v_mid(nEdges);
        v_mid.setZero();

        // Set up input and outputs for metric
        SobolevCurves::ApplyMid(curves, v, v_mid);

        // Set up inputs and outputs for low-order part
        Eigen::VectorXd b_mid_adm(nEdges);
        b_mid_adm.setZero();
        Eigen::VectorXd b_mid_inadm(nEdges);
        b_mid_inadm.setZero();

        // Multiply inadmissible blocks
        MultiplyInadmissibleFracParallel(v_mid, b_mid_inadm, s);
        // Multiply admissible blocks
        MultiplyAdmissibleFracFast(v_mid, b_mid_adm, s);

        b_mid_adm += b_mid_inadm;

        SobolevCurves::ApplyMidTranspose(curves, b_mid_adm, b);
    }

    template <typename V3, typename Dest>
    void BlockClusterTree::MultiplyVector3(V3 &v, Dest &b) const
    {
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

    template <typename V3, typename Dest>
    void BlockClusterTree::MultiplyByFracLaplacian3(V3 &v, Dest &b, double s) const
    {
        Eigen::VectorXd constrVals;
        constrVals.setZero(b.rows());

        // Slice the input vector to get every x-coordinate
        Eigen::Map<const Eigen::VectorXd, 0, Eigen::InnerStride<3>> v_x(v.data(), nVerts);
        // Slice the output vector to get x-coordinates
        Eigen::Map<Eigen::VectorXd, 0, Eigen::InnerStride<3>> dest_x(b.data(), nVerts);
        // Multiply the input x-coords into the output x-coords
        MultiplyByFracLaplacian(v_x, dest_x, s);

        // Same thing for y-coordinates
        Eigen::Map<const Eigen::VectorXd, 0, Eigen::InnerStride<3>> v_y(v.data() + 1, nVerts);
        Eigen::Map<Eigen::VectorXd, 0, Eigen::InnerStride<3>> dest_y(b.data() + 1, nVerts);
        MultiplyByFracLaplacian(v_y, dest_y, s);

        // Same thing for z-coordinates
        Eigen::Map<const Eigen::VectorXd, 0, Eigen::InnerStride<3>> v_z(v.data() + 2, nVerts);
        Eigen::Map<Eigen::VectorXd, 0, Eigen::InnerStride<3>> dest_z(b.data() + 2, nVerts);
        MultiplyByFracLaplacian(v_z, dest_z, s);
    }

    template <typename V, typename Dest>
    void BlockClusterTree::MultiplyWithConstraints(V &v, Dest &b) const
    {
        if (!constraintsSet)
        {
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

    template <typename V3, typename Dest>
    void BlockClusterTree::MultiplyWithConstraints3(V3 &v, Dest &b) const
    {
        if (!constraintsSet)
        {
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
} // namespace LWS
