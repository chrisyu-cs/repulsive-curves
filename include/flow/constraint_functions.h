#pragma once

#include <Eigen/Sparse>

namespace LWS {
    class PolyCurveNetwork;

    class ConstraintFunctions {
        public:
        static void AddBarycenterTriplets3X(PolyCurveNetwork* curves, std::vector<Eigen::Triplet<double>> &triplets, int rowStart);
        static void AddEdgeLengthTriplets(PolyCurveNetwork* curves, std::vector<Eigen::Triplet<double>> &triplets, int rowStart);
        static void AddPinTriplets(PolyCurveNetwork* curves, std::vector<Eigen::Triplet<double>> &triplets, int rowStart);

        static void NegativeBarycenterViolation(PolyCurveNetwork* curves, Eigen::VectorXd &b, Eigen::VectorXd &targets, int rowStart);
        static void NegativeEdgeLengthViolation(PolyCurveNetwork* curves, Eigen::VectorXd &b, Eigen::VectorXd &targets, int rowStart);
        static void NegativePinViolation(PolyCurveNetwork* curves, Eigen::VectorXd &b, Eigen::VectorXd &targets, int rowStart);
    };
}