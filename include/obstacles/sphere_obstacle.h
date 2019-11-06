#pragma once

#include "obstacles/obstacle.h"

namespace LWS {
    class SphereObstacle : public Obstacle {
        public:
        Vector3 center;
        double radius;
        SphereObstacle(Vector3 c, double r);
        virtual ~SphereObstacle();
        virtual void AddGradient(PolyCurveNetwork* curves, Eigen::MatrixXd &gradient, double alpha, double beta);

        inline Vector3 ClosestPoint(Vector3 input) {
            Vector3 dir = (input - center).normalize();
            return center + radius * dir;
        }
    };
}

