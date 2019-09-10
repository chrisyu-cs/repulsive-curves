#pragma once

#include "geometrycentral/utilities/vector3.h"
#include "Eigen/Sparse"
#include "poly_curve.h"

namespace LWS {
    class PolyCurveGroupHierarchy {
        public:
        std::vector<PolyCurveGroup*> levels;
        std::vector<MultigridOperator> prolongationOps;
        std::vector<MultigridOperator> sparsifyOps;

        PolyCurveGroupHierarchy(PolyCurveGroup* topLevel, size_t numLevels);
        ~PolyCurveGroupHierarchy();

        void AddNextLevel();

    };
}
