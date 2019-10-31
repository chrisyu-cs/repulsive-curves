#pragma once

#include <string>
#include <vector>

#include "geometrycentral/utilities/vector3.h"

namespace LWS {
    namespace CurveIO {
        using namespace geometrycentral;

        std::vector<std::string> split(const std::string& s, char delimiter);

        void readVerticesAndEdges(std::string fname, std::vector<Vector3> &all_positions,
        std::vector<std::array<size_t, 2>> &all_edges);
    }
}