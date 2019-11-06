#pragma once

#include "poly_curve_network.h"

namespace LWS {

    class Obstacle {
        public:
        virtual ~Obstacle() = 0;
        virtual void AddGradient(PolyCurveNetwork* curves, Eigen::MatrixXd &gradient, double alpha, double beta) = 0;
    };

}