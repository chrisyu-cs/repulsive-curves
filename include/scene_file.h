#pragma once

#include <string>
#include <vector>
#include <sstream>
#include <algorithm>
#include <iterator>

#include "flow/gradient_constraint_enum.h"
#include "extra_potentials.h"

namespace LWS {

    struct ObstacleData {
        std::string filename;
        double weight;
    };

    class SceneData {
        public:
        std::string curve_filename;
        std::vector<ObstacleData> obstacles;
        std::vector<ConstraintType> constraints;
        std::vector<int> pinnedVertices;
        std::vector<int> pinnedTangents;
        std::vector<int> surfaceConstrainedVertices;
        double tpe_alpha;
        double tpe_beta;
        double tpe_weight;
        std::vector<CurvePotential> extraPotentials;
    };

    template <class Container>
    void splitString(const std::string& str, Container& cont, char delim = ' ')
    {
        std::stringstream ss(str);
        std::string token;
        while (std::getline(ss, token, delim)) {
            cont.push_back(token);
        }
    }

    std::string getDirectoryFromPath(std::string str);

    SceneData ParseSceneFile(std::string filename);

}