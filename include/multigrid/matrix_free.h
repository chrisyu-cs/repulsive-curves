#pragma once

#include <iostream>
#include <Eigen/Core>
#include <Eigen/Dense>
#include <Eigen/IterativeLinearSolvers>
#include <unsupported/Eigen/IterativeSolvers>

#include "vector_multiplier.h"

using Eigen::SparseMatrix;

namespace LWS {
  namespace Product {

    // Matrix-free wrapper for a VectorMultiplier, for use with Eigen iterative solvers.
    template<typename T>
    class MatrixReplacement : public Eigen::EigenBase<MatrixReplacement<T>> {
    public:
      // Required typedefs, constants, and method:
      typedef double Scalar;
      typedef double RealScalar;
      typedef int StorageIndex;
      enum {
        ColsAtCompileTime = Eigen::Dynamic,
        MaxColsAtCompileTime = Eigen::Dynamic,
        IsRowMajor = false
      };
      Eigen::Index rows() const { return nRows; }
      Eigen::Index cols() const { return nRows; }
      template<typename Rhs>
      Eigen::Product<MatrixReplacement<T>,Rhs,Eigen::AliasFreeProduct> operator*(const Eigen::MatrixBase<Rhs>& x) const {
        return Eigen::Product<MatrixReplacement<T>,Rhs,Eigen::AliasFreeProduct>(*this, x.derived());
      }
      // Custom API:
      MatrixReplacement() {}
      MatrixReplacement(VectorMultiplier<T>* m, int n) { multiplier = m; nRows = n; }

      double nRows;
      VectorMultiplier<T>* multiplier;
    };
  }
}

template<typename T>
using _LWS_MR = LWS::Product::MatrixReplacement<T>;

namespace Eigen {
namespace internal {
  // MatrixReplacement looks-like a SparseMatrix, so let's inherits its traits:
  template<typename T>
  struct traits<_LWS_MR<T>> :  public Eigen::internal::traits<Eigen::SparseMatrix<double> >
  {};
}
}

// Implementation of MatrixReplacement * Eigen::DenseVector though a specialization of internal::generic_product_impl:
namespace Eigen {
namespace internal {
  template<typename Rhs, typename T>
  struct generic_product_impl<_LWS_MR<T>, Rhs, SparseShape, DenseShape, GemvProduct> // GEMV stands for matrix-vector
  : generic_product_impl_base<_LWS_MR<T>,Rhs,generic_product_impl<_LWS_MR<T>,Rhs> >
  {
    typedef typename Product<_LWS_MR<T>,Rhs>::Scalar Scalar;
    template<typename Dest>
    static void scaleAndAddTo(Dest& dst, const _LWS_MR<T>& lhs, const Rhs& rhs, const Scalar& alpha)
    {
        // This method should implement "dst += alpha * lhs * rhs" inplace,
        // however, for iterative solvers, alpha is always equal to 1, so let's not bother about it.
        assert(alpha==Scalar(1) && "scaling is not implemented");
        EIGEN_ONLY_USED_FOR_DEBUG(alpha);
        // Do the multiplication
        Eigen::VectorXd res;
        res.setZero(rhs.rows());
        lhs.multiplier->Multiply(rhs, res);
        dst += res;
    }
  };
}
}