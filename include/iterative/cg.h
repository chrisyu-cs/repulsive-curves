#pragma once

#include "../spatial/block_cluster_tree.h"

namespace LWS {

    class ConjugateGradient {

        public:
        // Solve the system Ax = b, using whatever is in x as an initial guess.
        // Stores the solution in x once finished.
        static void CGSolve(BlockClusterTree &A, std::vector<double> &x, std::vector<double> &b);
        // One iteration of conjugate gradient.
        // Returns true if iterations have converged; false otherwise.
        static bool CGIteration(BlockClusterTree &A, std::vector<double> &x_k, std::vector<double> &r_k, std::vector<double> &p_k);

    };

}
