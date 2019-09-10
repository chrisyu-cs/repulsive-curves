#include "poly_curve_hierarchy.h"

namespace LWS {

    PolyCurveGroupHierarchy::PolyCurveGroupHierarchy(PolyCurveGroup* topLevel, size_t numLevels) {
        levels.push_back(topLevel);

        while (levels.size() < numLevels) {
            AddNextLevel();
        }
    }

    PolyCurveGroupHierarchy::~PolyCurveGroupHierarchy() {
        for (size_t i = 1; i < levels.size(); i++) {
            delete levels[i];
        }
    }

    void PolyCurveGroupHierarchy::AddNextLevel() {
        PolyCurveGroup* lastLevel = levels[levels.size() - 1];
        MultigridOperator prolongOp, sparsifyOp;
        PolyCurveGroup* nextLevel = lastLevel->Coarsen(prolongOp, sparsifyOp);
        prolongOp.lowerSize = nextLevel->NumVertices();
        prolongOp.upperSize = lastLevel->NumVertices();

        prolongationOps.push_back(prolongOp);
        sparsifyOps.push_back(prolongOp);
        levels.push_back(nextLevel);

        std::cout << prolongOp.matrices[prolongOp.matrices.size() - 1].M.toDense() << std::endl;
    }
}
