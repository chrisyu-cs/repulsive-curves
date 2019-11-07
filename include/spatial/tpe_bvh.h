#pragma once

#include "spatial_tree.h"
#include "../tpe_energy_sc.h"
#include "geometrycentral/utilities/vector2.h"
#include "geometrycentral/surface/halfedge_mesh.h"
#include "geometrycentral/surface/vertex_position_geometry.h"
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
        
        inline double nodeRatio(double d) {
            // Compute diagonal distance from corner to corner
            // double diag = norm(maxCoords.position - minCoords.position);
            Vector3 diag = maxCoords.position - minCoords.position;
            double maxCoord = fmax(diag.x, fmax(diag.y, diag.z));

            Vector3 tanDiag = maxCoords.tangent - minCoords.tangent;
            double maxTanCoord = fmax(tanDiag.x, fmax(tanDiag.y, tanDiag.z));

            // return (maxCoord * sqrt(3.0) / 2) / d;
            return fmax((maxCoord * sqrt(3.0) / 2) / d, maxTanCoord * sqrt(3.0) / 2);
        }

        inline bool testTangent() {
            Vector3 tanDiag = maxCoords.tangent - minCoords.tangent;
            double maxTanCoord = fmax(tanDiag.x, fmax(tanDiag.y, tanDiag.z));
            double r = maxTanCoord * sqrt(3.0) / 2;
            return r < thresholdTheta;
        }

        template<typename T>
        void setLeafData(T &curves);

        int splitAxis;
        double splitPoint;
        VertexBody6D body;
        bool isEmpty;
        bool isLeaf;
        PosTan minCoords;
        PosTan maxCoords;
        double thresholdTheta;
    };

    template<typename T>
    void BVHNode3D::setLeafData(T &curves) {
        std::cerr << "Type not supported" << std::endl;
        throw 1;
    }
    
    template<>
    inline void BVHNode3D::setLeafData(PolyCurveNetwork* &curves) {
        if (body.type == BodyType::Vertex) {
            CurveVertex* p = curves->GetVertex(body.elementIndex);
            body.mass = p->DualLength();
            body.pt.position = p->Position();
            body.pt.tangent = p->Tangent();
        }
        else if (body.type == BodyType::Edge) {
            CurveEdge* p1 = curves->GetEdge(body.elementIndex);

            // Mass of an edge is its length
            body.mass = p1->Length();
            // Use midpoint as center of mass
            body.pt.position = p1->Midpoint();
            // Tangent direction is normalized edge vector
            body.pt.tangent = p1->Tangent();
        }

        totalMass = body.mass;
        centerOfMass = body.pt.position;
        averageTangent = body.pt.tangent;

        minCoords = PosTan{body.pt.position, body.pt.tangent};
        maxCoords = minCoords;
    }

    BVHNode3D* CreateBVHFromCurve(PolyCurveNetwork *curves);
    BVHNode3D* CreateEdgeBVHFromCurve(PolyCurveNetwork *curves);
    BVHNode3D* CreateBVHFromMesh(std::shared_ptr<geometrycentral::surface::HalfedgeMesh> &mesh,
        std::shared_ptr<geometrycentral::surface::VertexPositionGeometry> &geom);
}