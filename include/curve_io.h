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

        void writeOBJLineElements(std::string fname, const std::vector<Vector3> &all_positions,
        const std::vector<std::vector<size_t> > &components);
        // Writes a collection of curves to an OBJ file, where they are described as
        // line elements, i.e., each curve becomes a single line starting with 'l'.
        // Closed curves are not handled directly by OBJ, and should simply be passed
        // in with a repeated start/end vertex.
    }
}
