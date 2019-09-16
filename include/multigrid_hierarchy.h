#pragma once

#include "poly_curve.h"

#include "geometrycentral/utilities/vector3.h"
#include "Eigen/Sparse"

namespace LWS {
    class MultigridHierarchy {
        public:
        std::vector<PolyCurveGroup*> levels;
        std::vector<MultigridOperator> prolongationOps;
        std::vector<MultigridOperator> sparsifyOps;

        MultigridHierarchy(PolyCurveGroup* topLevel, size_t numLevels);
        ~MultigridHierarchy();

        void AddNextLevel();

        // Solve Gx = b, where G is the Sobolev Gram matrix of the top-level curve.
        Eigen::VectorXd VCycleSolve(Eigen::VectorXd b, double sepCoeff, double alpha, double beta, MultigridMode mode);

    };
}
