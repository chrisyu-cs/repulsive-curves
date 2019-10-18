#pragma once

#include "vertex_body.h"
#include "poly_curve_network.h"

#include "Eigen/Core"

namespace LWS {
    class SpatialTree {
        public:
        // Virtual destructor
        virtual ~SpatialTree() = 0;

        virtual void refreshWeightsVector(PolyCurveNetwork* curves, BodyType bType) = 0;

        // Recursively recompute all centers of mass in this tree
        virtual void recomputeCentersOfMass(PolyCurveNetwork* curves) = 0;
        
        // Compute the total energy contribution from a single vertex
        virtual void accumulateVertexEnergy(double &result, CurveVertex* &i_pt,
            PolyCurveNetwork* curves, double alpha, double beta) = 0;

        // Compute the total TPE gradient at a single vertex and its neighbors
        virtual void accumulateTPEGradient(Eigen::MatrixXd &gradients, CurveVertex* &i_pt,
            PolyCurveNetwork* curves, double alpha, double beta) = 0;

        // Use the given spatial tree to compute the TPE gradient with Barnes-Hut.
        static void TPEGradientBarnesHut(PolyCurveNetwork* curveNetwork, SpatialTree *root,
        Eigen::MatrixXd &gradients, double alpha, double beta);

        // Use the given spatial tree to compute the TPE energy with Barnes-Hut.
        static double TPEnergyBH(PolyCurveNetwork* curveNetwork, SpatialTree *root, double alpha, double beta);
    };
}