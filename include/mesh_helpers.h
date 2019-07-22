#include "polyscope/polyscope.h"

#include "geometrycentral/surface/vertex_position_geometry.h"
#include "geometrycentral/surface/halfedge_mesh.h"
#include "geometrycentral/surface/polygon_soup_mesh.h"

namespace LWS {

    using namespace geometrycentral;

    Vector3 Barycenter(surface::HalfedgeMesh* mesh, surface::VertexPositionGeometry* geom);
    double NormOfVectors(std::vector<Vector3> &vecs);

    void PrecomputeAreas(surface::HalfedgeMesh* mesh, surface::VertexPositionGeometry* geom, surface::VertexData<double> &data);
    void PrecomputeMeanCurvatures(surface::HalfedgeMesh* mesh, surface::VertexPositionGeometry* geom, surface::VertexData<double> &data);
    void PrecomputeGaussCurvatures(surface::HalfedgeMesh* mesh, surface::VertexPositionGeometry* geom, surface::VertexData<double> &data);

}