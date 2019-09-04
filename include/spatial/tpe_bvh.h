#pragma once

#include "spatial_tree.h"
#include "../tpe_energy_sc.h"
#include "geometrycentral/utilities/vector2.h"
#include "utils.h"

namespace LWS {

    struct BHPlotData {
        double theta;
        double error;
        double gradientNorm;
        double minWidth;
        double maxWidth;
    };

    class BVHNode3D : public SpatialTree {
        public:
        // Build a BVH of the given points
        BVHNode3D(std::vector<VertexBody6D> &points, int axis);
        virtual ~BVHNode3D();

        double totalMass;
        Vector3 centerOfMass;
        Vector3 averageTangent;

        void findCurveSegments(std::vector<VertexBody6D> &points, PolyCurveGroup* curves);

        // Recursively recompute all centers of mass in this tree
        virtual void recomputeCentersOfMass(PolyCurveGroup* curves);
        // Compute the total energy contribution from a single vertex
        virtual void accumulateVertexEnergy(double &result, PointOnCurve &i_pt, PolyCurveGroup* curves, double alpha, double beta);
        virtual void accumulateTPEGradient(Eigen::MatrixXd &gradients, PointOnCurve &i_pt, 
            PolyCurveGroup* curves, double alpha, double beta);
        int NumElements();
        
        virtual double bodyEnergyEvaluation(PointOnCurve &i_pt, double alpha, double beta);
        virtual Vector3 bodyForceEvaluation(PointOnCurve &i_pt, double alpha, double beta);

        Vector3 exactGradient(PointOnCurve basePoint, PolyCurveGroup* curves, double alpha, double beta);
        void testGradientSingle(std::vector<BHPlotData> &out, PointOnCurve basePoint,
            PolyCurveGroup* curves, double alpha, double beta);

        PosTan minBound();
        PosTan maxBound();
        Vector3 BoxCenter();
        std::vector<BVHNode3D*> children;

        void accumulateChildren(std::vector<VertexBody6D> &result);
        Vector2 viewspaceBounds(Vector3 point);

        private:
        int numElements;
        bool shouldUseCell(Vector3 vertPos);
        double AxisSplittingPlane(std::vector<VertexBody6D> &points, int axis);
        inline double nodeRadius();
        void setLeafData(PolyCurveGroup* curves);

        int splitAxis;
        double splitPoint;
        VertexBody6D body;
        bool isEmpty;
        bool isLeaf;
        PosTan minCoords;
        PosTan maxCoords;
        double thresholdTheta;
    };

    BVHNode3D* CreateBVHFromCurve(PolyCurveGroup *curves);
    BVHNode3D* CreateEdgeBVHFromCurve(PolyCurveGroup *curves);
}