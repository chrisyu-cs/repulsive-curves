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

            int nVerts = curves->NumVertices();
            dualLengths.setOnes(3 * nVerts + 3);

            for (int i = 0; i < nVerts; i++)
            {
                CurveVertex* v_i = curves->GetVertex(i);
                int ind = v_i->GlobalIndex();
                double len = v_i->DualLength();
                dualLengths(3 * ind) = len;
                dualLengths(3 * ind + 1) = len;
                dualLengths(3 * ind + 2) = len;
            }
        }

        SparseHs::~SparseHs()
        {
            delete bct;
            delete bvh;
        }

        void SparseHs::initFactorizedLaplacian() const
        {
            double totalLength = curves->TotalLength();            
            size_t nVerts = curves->NumVertices();
            std::vector<ETriplet> triplets, triplets3x;
            for (size_t i = 0; i < nVerts; i++)
            {
                CurveVertex *v_i = curves->GetVertex(i);
                size_t nEdges = v_i->numEdges();
                double rowSum = 0;

                double len_i = v_i->DualLength();

                for (size_t j = 0; j < nEdges; j++)
                {
                    CurveVertex *v_neighbor = v_i->edge(j)->Opposite(v_i);
                    double len = (v_i->Position() - v_neighbor->Position()).norm();
                    double wt = 1.0 / len;
                    triplets.push_back(ETriplet(v_i->GlobalIndex(), v_neighbor->GlobalIndex(), -wt));
                    rowSum += wt;
                }

                triplets.push_back(ETriplet(v_i->GlobalIndex(), v_i->GlobalIndex(), rowSum + 1e-8));

                // Add a barycenter triplet
                triplets.push_back(ETriplet(nVerts, v_i->GlobalIndex(), len_i));
                // triplets.push_back(ETriplet(nVerts, v_i->GlobalIndex(), 1));
            }

            // Expand the matrix by 3x
            TripleTriplets(triplets, triplets3x);

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
            bct->MultiplyWithConstraints3(v, dst);
        }

        void SparseHs::ApplyLML(const Eigen::VectorXd &v, Eigen::VectorXd &dst) const
        {
            if (!factorizedLaplacian.initialized)
            {
                initFactorizedLaplacian();
            }

            int nVerts = curves->NumVertices();
            // Eigen::VectorXd weighted = dualLengths.diagonal() * v;

            dst = factorizedLaplacian.Solve(v);
            double s = get_s();
            bct->MultiplyByFracLaplacian3(dst, dst, 2 - s);

            // Re-zero out Lagrange multipliers, since the first solve
            // will have left some junk in them
            for (int i = 3 * curves->NumVertices(); i < dst.rows(); i++)
            {
                dst(i) = 0;
            }

            // dst = dualLengths.diagonal() * v;
            dst = factorizedLaplacian.Solve(dst);
        }

    } // namespace preconditioned
} // namespace LWS
