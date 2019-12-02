#pragma once

#include "poly_curve_network.h"

namespace LWS {

    class Obstacle {
        public:
        virtual ~Obstacle() = 0;
        virtual void AddGradient(PolyCurveNetwork* curves, Eigen::MatrixXd &gradient) = 0;
        virtual double ComputeEnergy(PolyCurveNetwork* curves) = 0; 
    };

}