#pragma once 

#include "tpe_energy_sc.h"
#include <Eigen/Core>
#include <Eigen/Dense>
#include "product/block_cluster_tree.h"
#include "utils.h"

namespace LWS {
    struct EdgePositionPair {
        Vector3 f1;
        Vector3 f2;
    };

    class SobolevCurves {
        public:
        // Computes the local matrix contribution from a pair of edges sampled at
        // parametric points t1 and t2, where each edge has vertex positions (f1, f2)
        // and scalar function values (u1, u2) at its two endpoints.
        // The order of the indices is (e1.v1, e1.v2, e2.v1, e2.v2).
        // Writes the 4x4 matrix to the argument "out".
        static void LocalMatrix(EdgePositionPair e1, EdgePositionPair e2, double t1, double t2,
            double alpha, double beta, double out[4][4]);

        static void KfMatrix(EdgePositionPair e1, EdgePositionPair e2, double t1, double t2,
            double alpha, double beta, double out[4][4]);

        static void AddToFirst(double acc[4][4], double add[4][4], double weight);

        static void IntegrateLocalMatrix(EdgePositionPair e1, EdgePositionPair e2,
            double alpha, double beta, double out[4][4]);

        static void AddEdgePairContribution(PolyCurveGroup* loop, double alpha, double beta,
            PointOnCurve p1, PointOnCurve p2, Eigen::MatrixXd &A);

        static void AddEdgePairContributionLow(PolyCurveGroup* loop, double alpha, double beta,
            PointOnCurve p1, PointOnCurve p2, Eigen::MatrixXd &A);
        
        // Fills the global Sobolev-Slobodeckij Gram matrix.
        static void FillGlobalMatrix(PolyCurveGroup* loop, double alpha,
            double beta, Eigen::MatrixXd &A);
        
        // Computes the inner product of the two given vectors, under the metric
        // whose Gram matrix is given.
        static double SobolevDot(std::vector<Vector3> &as, std::vector<Vector3> &bs, Eigen::MatrixXd &J);

        static double SobolevDotMV(std::vector<Vector3> &as, std::vector<Vector3> &bs, BlockClusterTree* tree);

        // Normalizes the given vector so that it has norm 1 under the metric
        // whose Gram matrix is given.
        static void SobolevNormalize(std::vector<Vector3> &as, Eigen::MatrixXd &J);

        // Compute the projection of A onto the orthogonal complement of B under the given metric.
        static void SobolevOrthoProjection(std::vector<Vector3> &as, std::vector<Vector3> &bs, Eigen::MatrixXd &J);

        // Compute the projection of A onto the orthogonal complement of B,
        // using the block cluster tree instead of a matrix.
        static void SobolevOrthoProjectionMV(std::vector<Vector3> &as, std::vector<Vector3> &bs, BlockClusterTree *tree);

        // Assemble the Df differential operator as a matrix
        static void DfMatrix(PolyCurveGroup* loop, Eigen::MatrixXd &out);

        // Assemble the dense Sobolev Gram matrix on edges.
        static void FillAfMatrix(PolyCurveGroup* loop, double alpha, double beta, Eigen::MatrixXd &A);

        static void FillGlobalMatrixSlow(PolyCurveGroup* loop, double alpha,
            double beta, Eigen::MatrixXd &A);

        // Map a scalar-valued function on vertices to a gradient-valued function on edges
        static void ApplyDf(PolyCurveGroup* loop, std::vector<double> &as, std::vector<Vector3> &out);
        // Multiply by the transpose of the above Df (taking edge values back to vertices)
        static void ApplyDfTranspose(PolyCurveGroup* loop, std::vector<Vector3> &es, std::vector<double> &out);

        static void MultiplyComponents(Eigen::MatrixXd &A, std::vector<Vector3> &x, std::vector<Vector3> &out);

        static void ApplyGramMatrix(PolyCurveGroup* loop, std::vector<double> &as, double alpha, double beta, std::vector<double> &out);

        static double SobolevDotIterative(PolyCurveGroup* loop, std::vector<Vector3> &as, std::vector<Vector3> &bs);
    };
}