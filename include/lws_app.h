#pragma once

#include "polyscope/polyscope.h"
#include "polyscope/surface_mesh.h"

#include <iostream>
#include "geometrycentral/surface/vertex_position_geometry.h"
#include "geometrycentral/surface/halfedge_mesh.h"
#include "geometrycentral/surface/polygon_soup_mesh.h"

#include "vertexderivatives.h"
#include "boundary_derivatives.h"
#include "lws_options.h"
#include "lws_flow.h"
#include "mesh_helpers.h"
#include "lws_cluster.h"
#include "tpe_energy_sc.h"
#include "tpe_flow_sc.h"
#include "Eigen/SparseLU"
#include "poly_curve.h"

namespace LWS {
    class LWSApp {
        public:
        static LWSApp* instance;
        static double LWSVertexEnergy(surface::VertexPositionGeometry* geom, surface::Vertex vert);
        static Vector3 LWSVertexGradient(surface::VertexPositionGeometry* geom, surface::Vertex base, surface::Vertex other);
        void customWindow();
        void initSolver();
        void processFileOBJ(std::string filename);

        private:
        void centerMeshBarycenter(polyscope::SurfaceMesh* surface);
        void centerLoopBarycenter(PolyCurveGroup* curves);
        void UpdateMeshPositions();
        void UpdateCurvePositions();
        void outputFrame();
        void plotBHError(double alpha, double beta);
        std::string getCurveName(int i);
        
        std::unique_ptr<surface::HalfedgeMesh> mesh;
        std::unique_ptr<surface::VertexPositionGeometry> geom;
        std::string surfaceName;

        GradientSolver* solver;
        LWSEnergyFunction* lwsFunc;
        LWSClustering* clustering;
        PolyCurveGroup* curves;
        
        TPEFlowSolverSC* tpeSolver;
    };
}