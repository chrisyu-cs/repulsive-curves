#pragma once

#include "geometrycentral/surface/vertex_position_geometry.h"
#include "geometrycentral/surface/halfedge_mesh.h"
#include "geometrycentral/surface/halfedge_element_types.h"

namespace LWS {
    using namespace geometrycentral;

    // Typedef for a (pointer to a) function that evaluates an energy
    // at a single vertex.
    typedef double (*EnergyFunction)(surface::VertexPositionGeometry*, surface::Vertex);

    // Basic energy functions
    double VertexArea(surface::VertexPositionGeometry* geom, surface::Vertex vert);
    double VertexMeanCurvature(surface::VertexPositionGeometry* geom, surface::Vertex vert);
    double VertexGaussCurvature(surface::VertexPositionGeometry* geom, surface::Vertex vert);
    // LWS tubular energy function
    double VertexTubularEnergy(surface::VertexPositionGeometry* geom, surface::Vertex vert);

    double MeshEnergy(surface::VertexPositionGeometry* geom, surface::HalfedgeMesh* mesh, EnergyFunction energy);

    // Typedef for a (pointer to a) function that evaluates the gradient of some quantity
    // at one vertex wrt another vertex.
    typedef Vector3 (*GradientFunction)(surface::VertexPositionGeometry*, surface::Vertex, surface::Vertex);

    // Compute the derivative of the dual area of a vertex (base) wrt. another vertex (other
    Vector3 GradientArea(surface::VertexPositionGeometry* geom, surface::Vertex base, surface::Vertex other);

    double DihedralAngle(surface::VertexPositionGeometry* geom, surface::Edge base);
    // Compute the derivative of the dihedral angle of an edge (base) wrt. a vertex (other)
    Vector3 GradientDihedralAngle(surface::VertexPositionGeometry* geom, surface::Edge base, surface::Vertex other);
    Vector3 GradientDihedralNumerical(surface::VertexPositionGeometry* geom, surface::Edge base, surface::Vertex other, double h);

    // Compute the derivative of mean curvature at a vertex (base) wrt. another vertex (other)
    Vector3 GradientMeanCurvature(surface::VertexPositionGeometry* geom, surface::Vertex base, surface::Vertex other);

    // Compute the derivative of Gauss curvature at a vertex (base) wrt. another vertex (other)
    Vector3 GradientGaussCurvature(surface::VertexPositionGeometry* geom, surface::Vertex base, surface::Vertex other);

    // Compute derivative of (1, -2, 1) LWS energy
    Vector3 GradientTubularEnergy(surface::VertexPositionGeometry* geom, surface::Vertex base, surface::Vertex other);

    // Numerical derivatives for checking
    Vector3 GradientNumerical(surface::VertexPositionGeometry* geom, EnergyFunction energy, surface::Vertex base, surface::Vertex other, double h);
    Vector3 GradientAreaNumerical(surface::VertexPositionGeometry* geom, surface::Vertex base, surface::Vertex other);
    Vector3 GradientMeanNumerical(surface::VertexPositionGeometry* geom, surface::Vertex base, surface::Vertex other);
    Vector3 GradientGaussNumerical(surface::VertexPositionGeometry* geom, surface::Vertex base, surface::Vertex other);
}

