#pragma once

#include "polyscope/polyscope.h"
#include "polyscope/surface_mesh.h"

#include <iostream>
#include "geometrycentral/surface/vertex_position_geometry.h"
#include "geometrycentral/surface/halfedge_mesh.h"
#include "geometrycentral/surface/polygon_soup_mesh.h"

#include "lws_options.h"
#include "tpe_energy_sc.h"
#include "tpe_flow_sc.h"
#include "Eigen/SparseLU"
#include "poly_curve_network.h"

namespace LWS {
    class LWSApp {
        public:
        static LWSApp* instance;
        static double LWSVertexEnergy(surface::VertexPositionGeometry* geom, surface::Vertex vert);
        static Vector3 LWSVertexGradient(surface::VertexPositionGeometry* geom, surface::Vertex base, surface::Vertex other);
        void customWindow();
        void initSolver();
        void processFileOBJ(std::string filename);
        void processLoopFile(std::string filename);
        void processSceneFile(std::string filename);

        void AddMeshObstacle(std::string objName, Vector3 center);
        void AddPlaneObstacle(Vector3 center, Vector3 normal);
        void AddSphereObstacle(Vector3 center, double radius);

        void DisplayWireSphere(Vector3 center, double radius, std::string name);
        void DisplayPlane(Vector3 center, Vector3 normal, std::string name);
        void DisplayCurves(PolyCurveNetwork* curves, std::string name);
        void DisplayCyclicList(std::vector<Vector3> &positions, std::string name);
        std::string surfaceName;
        PolyCurveNetwork* curves;
        TPEFlowSolverSC* tpeSolver;

        private:
        void centerLoopBarycenter(PolyCurveNetwork* curves);
        void UpdateCurvePositions();
        void outputFrame();
        
        std::unique_ptr<surface::HalfedgeMesh> mesh;
        std::unique_ptr<surface::VertexPositionGeometry> geom;

    };
}