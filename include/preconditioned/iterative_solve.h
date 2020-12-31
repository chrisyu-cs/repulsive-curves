#pragma once

#include "preconditioned/sparse_preconditioner.h"
#include "preconditioned/hs_mat_replacement.h"

namespace LWS
{
    namespace preconditioned
    {
        void SolveIterative(PolyCurveNetwork *curves, const Eigen::VectorXd &v, Eigen::VectorXd &dst, double sep);

    }
} // namespace LWS
