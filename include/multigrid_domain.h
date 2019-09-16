#pragma once

#include "poly_curve.h"

namespace LWS {

    template<typename T>
    class MultigridDomain {
        public:
        MultigridDomain<T> Coarsen(MultigridOperator &sparsifyOp, MultigridOperator &prolongOp);
    };

}