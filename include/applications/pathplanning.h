#pragma once

#include "poly_curve_network.h"

namespace LWS {
    namespace Applications {
        std::vector<std::vector<Vector3>> SamplePathsAlongY(PolyCurveNetwork* curves, int nSamples);
        void SampleAndWritePaths(PolyCurveNetwork* curves, int nSamples, std::string fname);
    }
}
