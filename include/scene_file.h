#pragma once

#include <string>
#include <vector>
#include <sstream>
#include <algorithm>
#include <iterator>

#include "flow/gradient_constraint_enum.h"
#include "extra_potentials.h"
#include "implicit_surface.h"

namespace LWS {

    struct ObstacleData {
        std::string filename;
        double weight;
    };

    enum class PotentialType {
        Length, Area, VectorField
    };

    struct PotentialData {
        PotentialType type;
        double weight;
        std::string extraInfo;
    };

    class SceneData {
        public:
        std::string curve_filename;
        std::vector<ObstacleData> obstacles;
        std::vector<std::string> surfacesToShow;
        std::vector<ConstraintType> constraints;
        std::vector<int> pinnedVertices;
        std::vector<int> pinnedTangents;
        std::vector<int> surfaceConstrainedVertices;
        bool pinSpecialVertices;
        bool pinSpecialTangents;
        bool constrainAllToSurface;
        int subdivideLimit;
        int iterationLimit;
        double tpe_alpha;
        double tpe_beta;
        double tpe_weight;
        std::vector<PotentialData> extraPotentials;
        bool useLengthScale;
        double edgeLengthScale;
        bool useTotalLengthScale;
        double totalLengthScale;
        ImplicitSurface* constraintSurface;
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