#pragma once 

#include "tpe_energy_sc.h"
#include <Eigen/Core>
#include <Eigen/Dense>
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
        static void SobolevGramMatrix(PolyCurveGroup* loop, double alpha,
            double beta, Eigen::MatrixXd &A);

        // Fills the global Sobolev-Slobodeckij Gram matrix, plus an
        // extra row for a barycenter constraint.
        static void SobolevPlusBarycenter(PolyCurveGroup* loop, double alpha,
            double beta, Eigen::MatrixXd &A);
        
        // Computes the inner product of the two given vectors, under the metric
        // whose Gram matrix is given.
        static double SobolevDot(Eigen::MatrixXd &as, Eigen::MatrixXd &bs, Eigen::MatrixXd &J);

        // Normalizes the given vector so that it has norm 1 under the metric
        // whose Gram matrix is given.
        static void SobolevNormalize(Eigen::MatrixXd &as, Eigen::MatrixXd &J);

        // Compute the projection of A onto the orthogonal complement of B under the given metric.
        static void SobolevOrthoProjection(Eigen::MatrixXd &as, Eigen::MatrixXd &bs, Eigen::MatrixXd &J);

        // Assemble the Df differential operator as a matrix
        static void DfMatrix(PolyCurveGroup* loop, Eigen::MatrixXd &out);

        // Assemble the dense Sobolev Gram matrix on edges.
        static void FillAfMatrix(PolyCurveGroup* loop, double alpha, double beta, Eigen::MatrixXd &A);

        // Map a scalar-valued function on vertices to a gradient-valued function on edges
        template<typename V, typename M>
        static void ApplyDf(PolyCurveGroup* loop, V &as, M &out);

        // Multiply by the transpose of the above Df (taking edge values back to vertices)
        template<typename V, typename M>
        static void ApplyDfTranspose(PolyCurveGroup* loop, M &es, V &out);

        static void MultiplyComponents(Eigen::MatrixXd &A, std::vector<Vector3> &x, std::vector<Vector3> &out);

        static void ApplyGramMatrix(PolyCurveGroup* loop, Eigen::VectorXd &as, double alpha, double beta, Eigen::VectorXd &out);
    };

    template<typename V, typename M>
    void SobolevCurves::ApplyDf(PolyCurveGroup* loop, V &as, M &out) {
        for (PolyCurve *c : loop->curves) {
            int nVerts = c->NumVertices();
            for (int i = 0; i < nVerts; i++) {
                PointOnCurve p1 = loop->GetCurvePoint(i);
                int i1 = loop->GlobalIndex(p1);
                PointOnCurve p2 = p1.Next();
                int i2 = loop->GlobalIndex(p2);

                // Positive weight on the second vertex, negative on the first
                double d12 = as(i2) - as(i1);
                Vector3 e12 = (p2.Position() - p1.Position());
                double length = norm(e12);

                Vector3 grad12 = (d12 * e12) / (length * length);

                SetRow(out, i1, grad12);
            }
        }
    } 

    template<typename V, typename M>
    void SobolevCurves::ApplyDfTranspose(PolyCurveGroup* loop, M &es, V &out) {
        for (PolyCurve *c : loop->curves) {
            int nVerts = c->NumVertices();
            for (int i = 0; i < nVerts; i++) {
                // This is the first vertex of the next edge, so negative weight
                PointOnCurve p1 = loop->GetCurvePoint(i);
                // This is the second vertex of the previous edge, so positive weight
                PointOnCurve prev = p1.Prev();
                int i_prev = loop->GlobalIndex(prev);
                int i_next = loop->GlobalIndex(p1);

                PointOnCurve next = p1.Next();

                // Compute the two weights
                Vector3 v_prev = (p1.Position() - prev.Position());
                double len_prev = norm(v_prev);
                double w_prev = 1.0 / len_prev;
                v_prev /= len_prev;
                
                // However, the "forward" edge had a negative weight, so it is also negative here
                Vector3 v_next = (next.Position() - p1.Position());
                double len_next = norm(v_next);
                double w_next = -1.0 / len_next;
                v_next /= len_next;

                double result = w_prev * dot(v_prev, SelectRow(es, i_prev)) + w_next * dot(v_next, SelectRow(es, i_next));
                out(i_next) += result;
            }
        }
    }
}