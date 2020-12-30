#include "preconditioned/sparse_preconditioner.h"
#include "utils.h"
#include "flow/constraint_functions.h"

namespace LWS
{
    namespace preconditioned
    {
        SparseHs::SparseHs(PolyCurveNetwork *curves_, double alpha_, double beta_, double theta_)
        {
            curves = curves_;
            alpha = alpha_;
            beta = beta_;
            bct_theta = theta_;

            initBCT();
        }

        SparseHs::~SparseHs()
        {
            delete bct;
            delete bvh;
        }

        void SparseHs::initFactorizedLaplacian() const
        {
            size_t nVerts = curves->NumVertices();
            std::vector<ETriplet> triplets, triplets3x;
            for (size_t i = 0; i < nVerts; i++)
            {
                CurveVertex *v_i = curves->GetVertex(i);
                size_t nEdges = v_i->numEdges();
                double rowSum = 0;

                for (size_t j = 0; j < nEdges; j++)
                {
                    CurveVertex *v_neighbor = v_i->edge(j)->Opposite(v_i);
                    triplets.push_back(ETriplet(i, j, -1));
                    rowSum += 1;
                }

                triplets.push_back(ETriplet(i, i, rowSum));
            }

            // Expand the matrix by 3x
            TripleTriplets(triplets, triplets3x);
            // Add barycenter constraint rows
            ConstraintFunctions::AddBarycenterTriplets3X(curves, triplets3x, 3 * nVerts);

            size_t nRows = 3 * nVerts + 3;
            // Pre-factorize the cotan Laplacian
            Eigen::SparseMatrix<double> L(nRows, nRows);
            L.setFromTriplets(triplets3x.begin(), triplets3x.end());
            factorizedLaplacian.Compute(L);
        }

        void SparseHs::initBCT() const
        {
            bvh = CreateEdgeBVHFromCurve(curves);
            bct = new BlockClusterTree(curves, bvh, bct_theta, alpha, beta);
            bct->PremultiplyAfFrac(get_s());

            VariableConstraintSet vcs(curves);
            bct->SetConstraints(vcs);
            bct->SetBlockTreeMode(BlockTreeMode::Matrix3AndConstraints);
        }

        void SparseHs::MultiplyOperator(const Eigen::VectorXd &v, Eigen::VectorXd &dst) const
        {
            bct->MultiplyVector3(v, dst);
        }

        void SparseHs::ApplyLML(const Eigen::VectorXd &v, Eigen::VectorXd &dst) const
        {
            if (!factorizedLaplacian.initialized)
            {
                initFactorizedLaplacian();
            }

            dst = factorizedLaplacian.Solve(v);
            bct->MultiplyByFracLaplacian(dst, dst, 4 - 2 * get_s());

            // Re-zero out Lagrange multipliers, since the first solve
            // will have left some junk in them
            for (int i = 3 * curves->NumVertices(); i < dst.rows(); i++)
            {
                dst(i) = 0;
            }

            dst = factorizedLaplacian.Solve(dst);
        }

    } // namespace preconditioned
} // namespace LWS
