#pragma once

#include "geometrycentral/surface/vertex_position_geometry.h"
#include "geometrycentral/surface/halfedge_mesh.h"
#include "geometrycentral/surface/halfedge_element_types.h"
#include "polyscope/surface_mesh.h"

#include "vertexderivatives.h"

namespace LWS {
    using namespace geometrycentral;

    struct HalfedgePair {
        surface::Halfedge next;
        surface::Halfedge prev;
        bool valid;
    };

    void TestBoundaryQuantities(surface::HalfedgeMesh* mesh, surface::VertexPositionGeometry* geom);
    HalfedgePair BoundaryEdges(surface::Vertex vert);

    double VertexBoundaryLength(surface::VertexPositionGeometry* geom, surface::Vertex vert);
    double VertexBoundaryCurvature(surface::VertexPositionGeometry* geom, surface::Vertex vert);
    double VertexBoundaryTorsion(surface::VertexPositionGeometry* geom, surface::Vertex vert);

    Vector3 GradientBLength(surface::VertexPositionGeometry* geom, surface::Vertex vert, surface::Vertex other);
    Vector3 GradientBCurvature(surface::VertexPositionGeometry* geom, surface::Vertex vert, surface::Vertex other);
    Vector3 GradientBTorsion(surface::VertexPositionGeometry* geom, surface::Vertex vert, surface::Vertex other);

    Vector3 GradientBLengthNum(surface::VertexPositionGeometry* geom, surface::Vertex vert, surface::Vertex other);
    Vector3 GradientBCurvatureNum(surface::VertexPositionGeometry* geom, surface::Vertex vert, surface::Vertex other);
    Vector3 GradientBTorsionNum(surface::VertexPositionGeometry* geom, surface::Vertex vert, surface::Vertex other);
}