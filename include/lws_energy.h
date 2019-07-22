#pragma once

#include "polyscope/polyscope.h"

#include "geometrycentral/surface/vertex_position_geometry.h"
#include "geometrycentral/surface/halfedge_mesh.h"
#include "geometrycentral/surface/halfedge_element_types.h"
#include "geometrycentral/numerical/suitesparse_utilities.h"

#include "vertexderivatives.h"
#include "lws_cluster.h"

namespace LWS {

    class LWSEnergyFunction {
        public:
        LWSEnergyFunction();
        LWSEnergyFunction(double coeff_a, double coeff_b, double coeff_c);
        ~LWSEnergyFunction();
        void SetCoefficients(double coeff_a, double coeff_b, double coeff_c);
        // Coefficients for linear function
        double a;
        double b;
        double c;
        // Per-vertex energy and gradients
        double BaseEnergy(surface::VertexPositionGeometry* geom, surface::Vertex vert);
        void PrecomputeEnergies(surface::HalfedgeMesh* mesh, surface::VertexPositionGeometry* geom);
        double VertexEnergy(surface::VertexPositionGeometry* geom, surface::Vertex vert);
        Vector3 VertexGradient(surface::VertexPositionGeometry* geom, surface::Vertex base, surface::Vertex other);
        void SetClustering(LWSClustering* c);

        private:
        surface::VertexData<double>* vertexEnergies;
        LWSClustering* clustering;

    };

}