#pragma once

#include "polyscope/polyscope.h"
#include "polyscope/surface_mesh.h"

#include "geometrycentral/surface/vertex_position_geometry.h"
#include "geometrycentral/surface/halfedge_mesh.h"
#include "geometrycentral/surface/halfedge_element_types.h"
#include "geometrycentral/numerical/suitesparse_utilities.h"

#include <vector>
#include <queue>
#include <unordered_set>

#include <Eigen/Core>

namespace LWS {

    using namespace geometrycentral;

    class LWSCluster {
        public:
        LWSCluster();
        LWSCluster(surface::HalfedgeMesh* mesh, surface::VertexPositionGeometry* geom);
        std::vector<surface::Vertex> vertices;
        int index;
        double a;
        double b;
        double c;
        void computeCoefficients(surface::VertexData<double> &As, surface::VertexData<double> &Hs, surface::VertexData<double> &Ks);
        surface::Vertex minEnergyVertex();
        double vertexEnergy(surface::Vertex v);
        double vertexEnergyPre(surface::Vertex v, surface::VertexData<double> &As, surface::VertexData<double> &Hs, surface::VertexData<double> &Ks);

        private:
        Eigen::Matrix3d clusterQuadric(surface::VertexData<double> &As, surface::VertexData<double> &Hs, surface::VertexData<double> &Ks);
        surface::HalfedgeMesh* mesh;
        surface::VertexPositionGeometry* geom;
    };

    typedef struct ClusterQueueEntryStruct {
        surface::Vertex vertex;
        LWSCluster* cluster;
        double priority;
    } ClusterQueueEntry;

    class CompareEntries {
        public:
        bool operator() (const ClusterQueueEntry &x, const ClusterQueueEntry &y);
    };

    class LWSClustering {
        public:
        LWSClustering();
        LWSClustering(polyscope::SurfaceMesh* surf, surface::HalfedgeMesh* mesh, surface::VertexPositionGeometry* geom);
        std::vector<LWSCluster*> clusters;
        void initializeClusters(int numClusters);
        void priorityFillClusters();
        void colorByClusters();
        surface::VertexData<int> vertexOwnership;
        private:
        surface::VertexData<double> areas;
        surface::VertexData<double> meanCurvatures;
        surface::VertexData<double> gaussCurvatures;
        void precomputeScalars(surface::HalfedgeMesh* mesh, surface::VertexPositionGeometry* geom);
        void printClusters();
        std::unordered_set<surface::Vertex> claimedVertices;
        std::priority_queue<ClusterQueueEntry, std::vector<ClusterQueueEntry>, CompareEntries>* queue;
        surface::HalfedgeMesh* mesh;
        surface::VertexPositionGeometry* geom;
        polyscope::SurfaceMesh* surface;
    };
}