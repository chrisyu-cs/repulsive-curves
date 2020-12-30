#pragma once

#include "poly_curve_network.h"
#include "product/block_cluster_tree.h"
#include "preconditioned/sparse_preconditioner.h"

class BCTMatrixReplacement;
using Eigen::SparseMatrix;

// Adapted from https://eigen.tuxfamily.org/dox/group__MatrixfreeSolverExample.html

namespace Eigen
{
    namespace internal
    {
        // MatrixReplacement looks-like a SparseMatrix, so let's inherits its traits:
        template <>
        struct traits<BCTMatrixReplacement> : public Eigen::internal::traits<Eigen::SparseMatrix<double>>
        {
        };
    } // namespace internal
} // namespace Eigen

class BCTMatrixReplacement : public Eigen::EigenBase<BCTMatrixReplacement>
{
public:
    // Required typedefs, constants, and method:
    typedef double Scalar;
    typedef double RealScalar;
    typedef int StorageIndex;
    enum
    {
        ColsAtCompileTime = Eigen::Dynamic,
        MaxColsAtCompileTime = Eigen::Dynamic,
        IsRowMajor = false
    };

    Eigen::Index rows() const
    {
        return hs->nRows();
    }

    Eigen::Index outerSize() const
    {
        return hs->nRows();
    }

    Eigen::Index innerSize() const
    {
        return hs->nRows();
    }

    Eigen::Index cols() const
    {
        return hs->nRows();
    }

    template <typename Rhs>
    Eigen::Product<BCTMatrixReplacement, Rhs, Eigen::AliasFreeProduct> operator*(const Eigen::MatrixBase<Rhs> &x) const
    {
        return Eigen::Product<BCTMatrixReplacement, Rhs, Eigen::AliasFreeProduct>(*this, x.derived());
    }

    // Custom API:
    BCTMatrixReplacement() {}

    void addHs(LWS::preconditioned::SparseHs *hs_)
    {
        hs = hs_;
    }

    const LWS::preconditioned::SparseHs *getHs() const
    {
        return hs;
    }

private:
    const LWS::preconditioned::SparseHs *hs;
};

namespace LWS
{

    namespace preconditioned
    {
        class CustomPreconditioner
        {
            typedef Eigen::Matrix<double, Eigen::Dynamic, 1> Vector;

        public:
            typedef typename Vector::StorageIndex StorageIndex;
            enum
            {
                ColsAtCompileTime = Eigen::Dynamic,
                MaxColsAtCompileTime = Eigen::Dynamic
            };
            CustomPreconditioner() {}

            template <typename MatrixType>
            explicit CustomPreconditioner(const MatrixType &fracL)
            {
                compute(fracL);
            }

            Eigen::Index rows() const { return hs->nRows(); }
            Eigen::Index cols() const { return hs->nRows(); }

            template <typename MatrixType>
            CustomPreconditioner &analyzePattern(const MatrixType &) { return *this; }

            template <typename MatrixType>
            CustomPreconditioner &factorize(const MatrixType &) { return *this; }

            template <typename MatrixType>
            CustomPreconditioner &compute(const MatrixType &fracL)
            {
                hs = fracL.getHs();
                return *this;
            }

            /** \internal */
            template <typename Rhs, typename Dest>
            void _solve_impl(const Rhs &b, Dest &x) const
            {
                std::cout << "  * GMRES iteration " << (count++) << "...\r" << std::flush;
                hs->ApplyLML(b, x);
            }

            template <typename Rhs>
            inline const Eigen::Solve<CustomPreconditioner, Rhs>
            solve(const Eigen::MatrixBase<Rhs> &b) const
            {
                return Eigen::Solve<CustomPreconditioner, Rhs>(*this, b.derived());
            }

            mutable size_t count = 0;
            const LWS::preconditioned::SparseHs *hs;

            Eigen::ComputationInfo info() { return Eigen::Success; }
        };
    } // namespace Hs
} // namespace rsurfaces

namespace Eigen
{
    namespace internal
    {

        template <typename Rhs>
        struct generic_product_impl<BCTMatrixReplacement, Rhs, SparseShape, DenseShape, GemvProduct> // GEMV stands for matrix-vector
            : generic_product_impl_base<BCTMatrixReplacement, Rhs, generic_product_impl<BCTMatrixReplacement, Rhs>>
        {
            typedef typename Product<BCTMatrixReplacement, Rhs>::Scalar Scalar;

            template <typename Dest>
            static void scaleAndAddTo(Dest &dst, const BCTMatrixReplacement &lhs, const Rhs &rhs, const Scalar &alpha)
            {
                // This method should implement "dst += alpha * lhs * rhs" inplace,
                // however, for iterative solvers, alpha is always equal to 1, so let's not bother about it.
                assert(alpha == Scalar(1) && "scaling is not implemented");
                EIGEN_ONLY_USED_FOR_DEBUG(alpha);

                const LWS::preconditioned::SparseHs *hs = lhs.getHs();

                Eigen::VectorXd product;
                product.setZero(rhs.rows());
                hs->MultiplyOperator(rhs, product);
                dst += product;
            }
        };

    } // namespace internal
} // namespace Eigen
