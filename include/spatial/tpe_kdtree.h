#pragma once

#include "spatial_tree.h"

namespace LWS {

    class KDTreeNode3D /*: public SpatialTree */ {
        public:
        // Note that a KD tree is constructed by supplying all points at the start.
        KDTreeNode3D(std::vector<VertexBody6D> &points, int axis, PosTan mins, PosTan maxs);
        virtual ~KDTreeNode3D();

        double AxisSplittingPlane(std::vector<VertexBody6D> &points, int axis, PosTan mins, PosTan maxs);

        double totalMass;
        Vector3 centerOfMass;
        Vector3 averageTangent;
        // Recursively recompute all centers of mass in this tree
        virtual void recomputeCentersOfMass(PolyCurveGroup* curves);
        // Compute the total energy contribution from a single vertex
        virtual void accumulateVertexEnergy(double &result, PointOnCurve &i_pt, PolyCurveGroup* curves, double alpha, double beta);
        virtual void accumulateTPEGradient(std::vector<Vector3> &gradients, PointOnCurve &i_pt, 
            PolyCurveGroup* curves, double alpha, double beta);

        private:
        bool shouldUseCell(Vector3 vertPos);
        void combineValuesFromChildren();
        inline double nodeRadius();
        int splitAxis;
        double splitPoint;
        VertexBody6D body;
        bool isEmpty;
        bool isLeaf;
        KDTreeNode3D *children[2];
        PosTan minCoords;
        PosTan maxCoords;
        double thresholdTheta;
    };

    KDTreeNode3D* CreateKDTreeFromCurve(PolyCurveGroup *curves);

}