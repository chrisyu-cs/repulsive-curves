#pragma once

#include "poly_curve_network.h"
#include "spatial/tpe_bvh.h"
#include "product/block_cluster_tree.h"

namespace LWS
{
    namespace preconditioned
    {
        struct SparseFactorization
        {
            Eigen::SimplicialLDLT<Eigen::SparseMatrix<double>> factor;
            size_t nRows = 0;
            bool initialized = false;

            inline void Compute(Eigen::SparseMatrix<double> M)
            {
                nRows = M.rows();
                initialized = true;
                factor.compute(M);
            }

            inline Eigen::VectorXd Solve(const Eigen::VectorXd &v)
            {
                if (!initialized)
                {
                    std::cerr << "Sparse factorization was not initialized before attempting to solve." << std::endl;
                    throw 1;
                }
                return factor.solve(v);
            }

            inline Eigen::VectorXd SolveWithMasses(const Eigen::VectorXd &v, Eigen::VectorXd &mass)
            {
                if (!initialized)
                {
                    std::cerr << "Sparse factorization was not initialized before attempting to solve." << std::endl;
                    throw 1;
                }
                // Eigen::VectorXd
                return factor.solve(v);
            }
        };

        class SparseHs
        {
        public:
            SparseHs(PolyCurveNetwork *curves_, double alpha_, double beta_, double theta_);
            ~SparseHs();

            inline double get_s() const
            {
                return (beta - 1) / alpha;
            }

            void initFactorizedLaplacian() const;
            void initBCT() const;

            void MultiplyOperator(const Eigen::VectorXd &v, Eigen::VectorXd &dst) const;
            void ApplyLML(const Eigen::VectorXd &v, Eigen::VectorXd &dst) const;

            inline size_t nRows() const
            {
                return 3 * curves->NumVertices() + 3;
            }

        private:
            Eigen::VectorXd dualLengths;
            PolyCurveNetwork *curves;
            double alpha, beta, bct_theta;
            mutable SparseFactorization factorizedLaplacian;
            mutable BVHNode3D *bvh;
            mutable BlockClusterTree *bct;
        };

        void SolveIterative(PolyCurveNetwork *curves, const Eigen::VectorXd &v, Eigen::VectorXd &dst);
    } // namespace preconditioned
} // namespace LWS