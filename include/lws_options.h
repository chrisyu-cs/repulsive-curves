#pragma once

#include "geometrycentral/surface/vertex_position_geometry.h"
#include "lws_flow.h"
#include "lws_energy.h"
#include "lws_cluster.h"
#include "boundary_loop.h"

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
    };
}