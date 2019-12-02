#pragma once

#include "obstacles/obstacle.h"

namespace LWS {
    class SphereObstacle : public Obstacle {
        public:
        Vector3 center;
        double radius;
        double p;
        SphereObstacle(Vector3 c, double r, double p_exp);
        virtual ~SphereObstacle();
        virtual void AddGradient(PolyCurveNetwork* curves, Eigen::MatrixXd &gradient);
        virtual double ComputeEnergy(PolyCurveNetwork* curves);

        inline Vector3 ClosestPoint(Vector3 input) {
            Vector3 dir = (input - center).normalize();
            return center + radius * dir;
        }
    };
}

