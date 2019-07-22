#include "mesh_helpers.h"
#include "vertexderivatives.h"

namespace LWS {
    
    Vector3 Barycenter(surface::HalfedgeMesh* mesh, surface::VertexPositionGeometry* geom) {
        Vector3 center = Vector3{0, 0, 0};
        double totalArea = 0;

        for (surface::Vertex v : mesh->vertices()) {
            double area = VertexArea(geom, v);
            totalArea += area;
            center += geom->vertexPositions[v] * area;
        }
        center /= totalArea;
        return center;
    }

    void PrecomputeAreas(surface::HalfedgeMesh* mesh, surface::VertexPositionGeometry* geom, surface::VertexData<double> &data) {
        for (surface::Vertex v : mesh->vertices()) {
            data[v] = VertexArea(geom, v);
        }
    }

    void PrecomputeMeanCurvatures(surface::HalfedgeMesh* mesh, surface::VertexPositionGeometry* geom, surface::VertexData<double> &data) {
        for (surface::Vertex v : mesh->vertices()) {
            data[v] = VertexMeanCurvature(geom, v);
        }
    }

    void PrecomputeGaussCurvatures(surface::HalfedgeMesh* mesh, surface::VertexPositionGeometry* geom, surface::VertexData<double> &data) {
        for (surface::Vertex v : mesh->vertices()) {
            data[v] = VertexGaussCurvature(geom, v);
        }
    }

    double NormOfVectors(std::vector<Vector3> &vecs) {
        double sum = 0;
        for (size_t i = 0; i < vecs.size(); i++) {
            Vector3 entry = vecs[i];
            sum += dot(entry, entry);
        }
        return sqrt(sum);
    }
}