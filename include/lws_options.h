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
        static bool useSobolev;
        static bool useMultigrid;
        static bool useBarnesHut;
        static bool normalizeView;

        static double tpeAlpha;
        static double tpeBeta;
    };
}
