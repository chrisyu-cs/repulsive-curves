#pragma once

#include "multigrid_operator.h"
#include "nullspace_projector.h"
#include <Eigen/Core>

namespace LWS {

    template<typename T, typename Mult>
    class MultigridDomain {
        public:
        typedef Mult MultType;
        virtual ~MultigridDomain() {}

        MultigridDomain<T, Mult>* Coarsen(MultigridOperator* prolongOp) const {
            return static_cast<T const&>(*this).Coarsen(prolongOp);
        }

        Mult* GetMultiplier() const {
            return static_cast<T const&>(*this).GetMultiplier();
        }

        Eigen::MatrixXd GetFullMatrix() const {
            return static_cast<T const&>(*this).GetFullMatrix();
        }

        Eigen::VectorXd DirectSolve(Eigen::VectorXd &b) const {
            return static_cast<T const&>(*this).DirectSolve(b);
        }

        int NumVertices() const {
            return static_cast<T const&>(*this).NumVertices();
        }
        
        int NumRows() const {
            return static_cast<T const&>(*this).NumRows();
        }

        MultigridOperator* MakeNewOperator() const {
            return static_cast<T const&>(*this).MakeNewOperator();
        }
        
        NullSpaceProjector* GetConstraintProjector() const {
            return static_cast<T const&>(*this).GetConstraintProjector();
        }
    };
}