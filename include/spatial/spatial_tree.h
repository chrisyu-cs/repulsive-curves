#pragma once

#include "vertex_body.h"
#include "../poly_curve.h"

#include "Eigen/Core"

namespace LWS {
    class SpatialTree {
        public:
        // Virtual destructor
        virtual ~SpatialTree() = 0;

        // Recursively recompute all centers of mass in this tree
        virtual void recomputeCentersOfMass(PolyCurveGroup* curves) = 0;
        
        // Compute the total energy contribution from a single vertex
        virtual void accumulateVertexEnergy(double &result, PointOnCurve &i_pt,
            PolyCurveGroup* curves, double alpha, double beta) = 0;

        // Compute the total TPE gradient at a single vertex and its neighbors
        virtual void accumulateTPEGradient(Eigen::MatrixXd &gradients, PointOnCurve &i_pt,
            PolyCurveGroup* curves, double alpha, double beta) = 0;
    };
}