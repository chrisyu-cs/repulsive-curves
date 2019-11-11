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
        BVHNode3D(std::vector<VertexBody6D> &points, int axis, BVHNode3D* root, bool splitTangents);
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
        void refreshWeightsVector(PolyCurveNetwork* curves, BodyType bType);

        // Recursively recompute all centers of mass in this tree
        template<typename T>
        void recomputeCentersOfMass(T &curves);
        
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

        inline int countBadNodes() {
            if (isEmpty) return 0;
            else if (isLeaf) {
                if (testTangent()) return 0;
                else return 1;
            }
            else {
                int sum = 0;
                if (!testTangent()) sum++;
                for (BVHNode3D* child : children) {
                    sum += child->countBadNodes();
                } 
                return sum;
            }
        }

        inline double nodeRatio(double d) {
            // Compute diagonal distance from corner to corner
            // double diag = norm(maxCoords.position - minCoords.position);
            Vector3 diag = maxCoords.position - minCoords.position;
            double maxCoord = fmax(diag.x, fmax(diag.y, diag.z));
            // double spatialR = diag.norm() / 2;
            return diag.norm() / d;

            // Vector3 tanDiag = maxCoords.tangent - minCoords.tangent;
            // double maxTanCoord = fmax(tanDiag.x, fmax(tanDiag.y, tanDiag.z));
            // double tangentR = tanDiag.norm() / 2;
            // return fmax(spatialR / d, tangentR);
        }

        bool shouldUseCell(Vector3 vertPos);

        private:
        int numElements;
        double AxisSplittingPlane(std::vector<VertexBody6D> &points, int axis);
        
        inline bool testTangent() {
            Vector3 tanDiag = maxCoords.tangent - minCoords.tangent;
            double r = tanDiag.norm() / 2;
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

    template<>
    inline void BVHNode3D::setLeafData(std::pair<std::shared_ptr<geometrycentral::surface::HalfedgeMesh>,
    std::shared_ptr<geometrycentral::surface::VertexPositionGeometry>> &pair) {
        std::shared_ptr<geometrycentral::surface::HalfedgeMesh> mesh = pair.first;
        std::shared_ptr<geometrycentral::surface::VertexPositionGeometry> geom = pair.second;

        using namespace geometrycentral;
        using namespace surface;

        if (body.type == BodyType::Vertex) {
            Vertex v = mesh->vertex(body.elementIndex);
            body.mass = geom->vertexDualAreas[v];
            body.pt.position = geom->vertexPositions[v];
            body.pt.tangent = geom->vertexNormals[v];
        }
        else {
            std::cerr << "Element types besides vertex are not supported for meshes" << std::endl;
            throw 1;
        }

        totalMass = body.mass;
        centerOfMass = body.pt.position;
        averageTangent = body.pt.tangent;

        minCoords = PosTan{body.pt.position, body.pt.tangent};
        maxCoords = minCoords;
    }
    
    template<typename T>
    inline void BVHNode3D::recomputeCentersOfMass(T &curves) {
        if (isEmpty) {
            totalMass = 0;
            numElements = 0;
        }
        // For a leaf, just set centers and bounds from the one body
        else if (isLeaf) {
            setLeafData(curves);
            numElements = 1;
        }
        else {
            // Recursively compute bounds for all children
            for (size_t i = 0; i < children.size(); i++) {
                children[i]->recomputeCentersOfMass(curves);
            }

            minCoords = children[0]->minCoords;
            maxCoords = children[0]->maxCoords;

            totalMass = 0;
            centerOfMass = Vector3{0, 0, 0};
            averageTangent = Vector3{0, 0, 0};
            
            // Accumulate max/min over all nonempty children
            for (size_t i = 0; i < children.size(); i++) {
                if (!children[i]->isEmpty) {
                    minCoords = postan_min(children[i]->minCoords, minCoords);
                    maxCoords = postan_max(children[i]->maxCoords, maxCoords);

                    totalMass += children[i]->totalMass;
                    centerOfMass += children[i]->centerOfMass * children[i]->totalMass;
                    averageTangent += children[i]->averageTangent * children[i]->totalMass;
                }
            }

            centerOfMass /= totalMass;
            averageTangent /= totalMass;

            averageTangent = averageTangent.normalize();

            numElements = 0;
            for (size_t i = 0; i < children.size(); i++) {
                numElements += children[i]->numElements;
            }
        }
    }

    BVHNode3D* CreateBVHFromCurve(PolyCurveNetwork *curves);
    BVHNode3D* CreateEdgeBVHFromCurve(PolyCurveNetwork *curves);
    BVHNode3D* CreateBVHFromMesh(std::shared_ptr<geometrycentral::surface::HalfedgeMesh> &mesh,
        std::shared_ptr<geometrycentral::surface::VertexPositionGeometry> &geom);
}