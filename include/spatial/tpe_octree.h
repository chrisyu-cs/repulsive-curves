#pragma once

#include "spatial_tree.h"

namespace LWS {

    class OctreeNode3D : public SpatialTree {
        public:
        OctreeNode3D(Vector3 center, double width);
        virtual ~OctreeNode3D();
        // Add the given body to the tree under this root
        void AddToSubtree(VertexBody b);
        // Recursively recompute all centers of mass in this tree
        virtual void recomputeCentersOfMass(PolyCurveGroup* curves);
        // Compute the total energy contribution from a single vertex
        virtual void accumulateVertexEnergy(double &result, PointOnCurve &i_pt, PolyCurveGroup* curves, double alpha, double beta);
        // Compute the total TPE gradient at a single vertex and its neighbors
        virtual void accumulateTPEGradient(std::vector<Vector3> &gradients, PointOnCurve &i_pt,
            PolyCurveGroup* curves, double alpha, double beta);
        // Aggregated values for this cell
        double totalMass;
        Vector3 centerOfMass;
        Vector3 averageTangent;

        private:
        // Check if this is empty
        inline bool isEmpty();
        // Find the octant in which the given position would belong
        int findOctant(Vector3 position);
        // Update the center of mass by adding a body of given mass at given pos
        void updateCenterOfMass(Vector3 pos, Vector3 tangent, double mass);
        // Make the given child if it doesn't already exist
        void createChildIfNonexistent(int childIndex);
        Vector3 getChildCenter(int octant);
        VertexBody body;
        bool isLeaf;
        OctreeNode3D *children[8];
        Vector3 nodeCenter;
        double nodeWidth;
        double theta;

    };

    OctreeNode3D* CreateOctreeFromCurve(PolyCurveGroup *curves);

}
