#pragma once

#include "../product/block_cluster_tree.h"

namespace LWS {

    class ConjugateGradient {

        public:
        // Solve the system Ax = b, using whatever is in x as an initial guess.
        // Stores the solution in x once finished.
        static void CGSolve(BlockClusterTree *tree, std::vector<double> &x, std::vector<double> &b);
        // Solve the 3 systems Ax = b obtained by separating x and b by components.
        // Solution is stored in x when finished.
        static void CGSolveComponents(BlockClusterTree *tree, std::vector<Vector3> &x, std::vector<Vector3> &b);
        // One iteration of conjugate gradient.
        // Returns true if iterations have converged; false otherwise.
        static bool CGIteration(BlockClusterTree *tree, std::vector<double> &x_k, std::vector<double> &r_k, std::vector<double> &p_k);

    };

}
