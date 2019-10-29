#pragma once

#include "spatial_tree.h"
#include "../tpe_energy_sc.h"
#include "geometrycentral/utilities/vector2.h"
#include "utils.h"

#include <fstream>
#include <Eigen/Core>

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
        static int globalID;
        int thisNodeID;
        int numNodes;

        inline void recursivelyAssignIDs() {
            thisNodeID = globalID++;
            for (BVHNode3D* child : children) {
                child->recursivelyAssignIDs();
            }
        }

        inline void assignIDs() {
            globalID = 1;
            recursivelyAssignIDs();
        }

        inline void printIDs(std::ofstream &stream, int parentID = 0) {
            if (parentID > 0) {
                stream << thisNodeID << ", " << parentID << std::endl;
            }
            for (BVHNode3D* child : children) {
                child->printIDs(stream, thisNodeID);
            }
        }

        // Build a BVH of the given points
        BVHNode3D(std::vector<VertexBody6D> &points, int axis, BVHNode3D* root);
        virtual ~BVHNode3D();

        double totalMass;
        Vector3 centerOfMass;
        Vector3 averageTangent;
        std::vector<int> clusterIndices;
        BVHNode3D* bvhRoot;
        Eigen::VectorXd fullMasses;

        // Fields for use by matrix-vector products; not used
        // by any BVH functions.
        double V_I;
        double B_I;
        double aIJ_VJ;

        inline void zeroMVFields() {
            V_I = 0;
            B_I = 0;
            aIJ_VJ = 0;
        }

        inline void recursivelyZeroMVFields() {
            zeroMVFields();
            if (!isLeaf) {
                for (BVHNode3D* child : children) {
                    child->recursivelyZeroMVFields();
                }
            }
        }

        void findCurveSegments(std::vector<VertexBody6D> &points, PolyCurveNetwork* curves);

        // Copy the new weights from the curves
        virtual void refreshWeightsVector(PolyCurveNetwork* curves, BodyType bType);
        // Recursively recompute all centers of mass in this tree
        virtual void recomputeCentersOfMass(PolyCurveNetwork* curves);
        // Compute the total energy contribution from a single vertex
        virtual void accumulateVertexEnergy(double &result, CurveVertex* &i_pt, PolyCurveNetwork* curves, double alpha, double beta);
        virtual void accumulateTPEGradient(Eigen::MatrixXd &gradients, CurveVertex* &i_pt, 
            PolyCurveNetwork* curves, double alpha, double beta);
        int NumElements();
        
        virtual double bodyEnergyEvaluation(CurveVertex* &i_pt, double alpha, double beta);
        virtual Vector3 bodyForceEvaluation(CurveVertex* &i_pt, double alpha, double beta);

        Vector3 exactGradient(CurveVertex* basePoint, PolyCurveNetwork* curves, double alpha, double beta);

        PosTan minBound();
        PosTan maxBound();
        Vector3 BoxCenter();
        std::vector<BVHNode3D*> children;

        void accumulateChildren(std::vector<VertexBody6D> &result);
        Vector2 viewspaceBounds(Vector3 point);

        inline void fillClusterMassVector(Eigen::VectorXd &w) {
            w.setZero(clusterIndices.size());
            for (size_t i = 0; i < clusterIndices.size(); i++) {
                w(i) = bvhRoot->fullMasses(clusterIndices[i]);
            }
        }

        inline bool IsLeaf() {
            return isLeaf;
        }
        inline bool IsEmpty() {
            return isEmpty;
        }
        inline int VertexIndex() {
            return body.elementIndex;
        }

        private:
        int numElements;
        bool shouldUseCell(Vector3 vertPos);
        double AxisSplittingPlane(std::vector<VertexBody6D> &points, int axis);
        inline double nodeRadius();
        void setLeafData(PolyCurveNetwork* curves);

        int splitAxis;
        double splitPoint;
        VertexBody6D body;
        bool isEmpty;
        bool isLeaf;
        PosTan minCoords;
        PosTan maxCoords;
        double thresholdTheta;
    };

    BVHNode3D* CreateBVHFromCurve(PolyCurveNetwork *curves);
    BVHNode3D* CreateEdgeBVHFromCurve(PolyCurveNetwork *curves);
}