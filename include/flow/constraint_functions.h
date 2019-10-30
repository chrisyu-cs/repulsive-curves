#pragma once

#include <Eigen/Sparse>

namespace LWS {
    class PolyCurveNetwork;

    class ConstraintFunctions {
        public:
        static void AddBarycenterTriplets3X(PolyCurveNetwork* curves, std::vector<Eigen::Triplet<double>> &triplets, int rowStart);
        static void AddEdgeLengthTriplets(PolyCurveNetwork* curves, std::vector<Eigen::Triplet<double>> &triplets, int rowStart);
    };
}