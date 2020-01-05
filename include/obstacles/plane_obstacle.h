#pragma once

#include "obstacles/obstacle.h"

namespace LWS {
    class PlaneObstacle : public Obstacle {
        public:
        Vector3 center;
        Vector3 normal;
        double p;
        double weight;
        PlaneObstacle(Vector3 c, Vector3 n, double p_exp, double wt);
        virtual ~PlaneObstacle();
        virtual void AddGradient(PolyCurveNetwork* curves, Eigen::MatrixXd &gradient);
        virtual double ComputeEnergy(PolyCurveNetwork* curves);

        inline Vector3 ClosestPoint(Vector3 input) {
            Vector3 fromInput = center - input;
            Vector3 normalProj = normal * dot(normal, fromInput);
            return input + normalProj;
        }
    };
}
