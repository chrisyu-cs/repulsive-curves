#pragma once 

#include "tpe_energy_sc.h"
#include <Eigen/Core>
#include <Eigen/Dense>
#include "utils.h"
#include "poly_curve_network.h"

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

        static void AddEdgePairContribution(PolyCurveNetwork* loop, double alpha, double beta,
            CurveEdge* p1, CurveEdge* p2, Eigen::MatrixXd &A);

        static void AddEdgePairContributionLow(PolyCurveNetwork* loop, double alpha, double beta,
            CurveEdge* p1, CurveEdge* p2, Eigen::MatrixXd &A);
        
        // Fills the global Sobolev-Slobodeckij Gram matrix.
        static void SobolevGramMatrix(PolyCurveNetwork* loop, double alpha,
            double beta, Eigen::MatrixXd &A, double diagEps = 0);

        // Fills the global Sobolev-Slobodeckij Gram matrix, plus an
        // extra row for a barycenter constraint.
        static void SobolevPlusBarycenter(PolyCurveNetwork* loop, double alpha,
            double beta, Eigen::MatrixXd &A);

        // Fills the global Sobolev-Slobodeckij Gram matrix, plus an
        // extra row for a barycenter constraint.
        static void SobolevPlusBarycenter3X(PolyCurveNetwork* loop, double alpha,
            double beta, Eigen::MatrixXd &A);

        // Fills the matrix A with edge length constraints, starting from the row baseIndex.
        static void AddEdgeLengthConstraints(PolyCurveNetwork* curves, Eigen::MatrixXd &A, int baseIndex);
        
        // Computes the inner product of the two given vectors, under the metric
        // whose Gram matrix is given.
        static double SobolevDot(Eigen::MatrixXd &as, Eigen::MatrixXd &bs, Eigen::MatrixXd &J);

        // Normalizes the given vector so that it has norm 1 under the metric
        // whose Gram matrix is given.
        static void SobolevNormalize(Eigen::MatrixXd &as, Eigen::MatrixXd &J);

        // Compute the projection of A onto the orthogonal complement of B under the given metric.
        static void SobolevOrthoProjection(Eigen::MatrixXd &as, Eigen::MatrixXd &bs, Eigen::MatrixXd &J);

        // Assemble the Df differential operator as a matrix
        static void DfMatrix(PolyCurveNetwork* loop, Eigen::MatrixXd &out);

        // Map a scalar-valued function on vertices to a gradient-valued function on edges
        template<typename V, typename M>
        static void ApplyDf(PolyCurveNetwork* loop, V &as, M &out);

        // Multiply by the transpose of the above Df (taking edge values back to vertices)
        template<typename V, typename M>
        static void ApplyDfTranspose(PolyCurveNetwork* loop, M &es, V &out);

        static void MultiplyComponents(Eigen::MatrixXd &A, std::vector<Vector3> &x, std::vector<Vector3> &out);
    };

    template<typename V, typename M>
    void SobolevCurves::ApplyDf(PolyCurveNetwork* loop, V &as, M &out) {
        int nEdges = loop->NumEdges();
        for (int i = 0; i < nEdges; i++) {
            CurveEdge* e = loop->GetEdge(i);
            int i_edge = e->GlobalIndex();
            CurveVertex* p1 = e->prevVert;
            CurveVertex* p2 = e->nextVert;
            int i1 = p1->GlobalIndex();
            int i2 = p2->GlobalIndex();

            // Positive weight on the second vertex, negative on the first
            double d12 = as(i2) - as(i1);
            Vector3 e12 = (p2->Position() - p1->Position());
            double length = norm(e12);

            Vector3 grad12 = (d12 * e12) / (length * length);

            SetRow(out, i_edge, grad12);
        }
    } 

    template<typename V, typename M>
    void SobolevCurves::ApplyDfTranspose(PolyCurveNetwork* loop, M &es, V &out) {
        int nVerts = loop->NumVertices();
        for (int i = 0; i < nVerts; i++) {
            double result = 0;
            CurveVertex* p1 = loop->GetVertex(i);
            int i_vert = p1->GlobalIndex();

            for (int e = 0; e < p1->numEdges(); e++) {
                CurveEdge* edge = p1->edge(e);
                int i_edge = edge->GlobalIndex();
                Vector3 v_edge = edge->Vector();
                double len_edge = norm(v_edge);
                double w_edge = 1.0 / len_edge;
                // If this edge is a "forward" edge,
                // then the weight is negative
                if (edge->prevVert == p1) {
                    w_edge *= -1;
                }
                v_edge /= len_edge;

                result += w_edge * dot(v_edge, SelectRow(es, i_edge));
            }
            out(i_vert) += result;
        }
    }
}