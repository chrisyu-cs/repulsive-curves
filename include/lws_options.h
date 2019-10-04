#pragma once

#include "geometrycentral/surface/vertex_position_geometry.h"

namespace LWS {
    class LWSOptions {
        public:
        static int iterationNum;
        static int frameNum;
        static bool outputFrames;
        static bool showWindow;
        static bool runLWSFlow;
        static bool runBoundaryFlow;
        static bool runTPE;
        static bool useSobalev;
        static bool useMultigrid;
    };
}