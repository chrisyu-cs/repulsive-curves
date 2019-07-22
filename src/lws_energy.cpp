#include "lws_energy.h"

using namespace geometrycentral::surface;

namespace LWS {

    // Default constructor goes to tubular energy
    LWSEnergyFunction::LWSEnergyFunction() {
        a = 1;
        b = -2;
        c = 1;
        vertexEnergies = 0;
        clustering = 0;
    }

    LWSEnergyFunction::~LWSEnergyFunction() {
        // We don't allocate anything, so nothing to do
    }

    LWSEnergyFunction::LWSEnergyFunction(double coeff_a, double coeff_b, double coeff_c) {
        SetCoefficients(coeff_a, coeff_b, coeff_c);
        vertexEnergies = 0;
        clustering = 0;
    }

    void LWSEnergyFunction::SetClustering(LWSClustering* c) {
        clustering = c;
    }

    void LWSEnergyFunction::SetCoefficients(double coeff_a, double coeff_b, double coeff_c) {
        a = coeff_a;
        b = coeff_b;
        c = coeff_c;
    }

    double LWSEnergyFunction::BaseEnergy(VertexPositionGeometry* geom, Vertex vert) {
        // If a clustering is assigned, use the coefficients from this vertex's cluster
        if (clustering) {
            LWSCluster* c = clustering->clusters[clustering->vertexOwnership[vert]];
            return c->a * VertexArea(geom, vert) + c->b * VertexMeanCurvature(geom, vert) + c->c * VertexGaussCurvature(geom, vert);
        }
        else {
            return a * VertexArea(geom, vert) + b * VertexMeanCurvature(geom, vert) + c * VertexGaussCurvature(geom, vert);
        }
    }

    void LWSEnergyFunction::PrecomputeEnergies(HalfedgeMesh* mesh, VertexPositionGeometry* geom) {
        if (!vertexEnergies) {
            vertexEnergies = new VertexData<double>(*mesh);
        }

        for (Vertex v : mesh->vertices()) {
            (*vertexEnergies)[v] = BaseEnergy(geom, v);
        }
    }

    double LWSEnergyFunction::VertexEnergy(VertexPositionGeometry* geom, Vertex vert) {
        double e = BaseEnergy(geom, vert);
        return e * e;
    }

    Vector3 LWSEnergyFunction::VertexGradient(VertexPositionGeometry* geom, Vertex base, Vertex other) {
        Vector3 gradA = GradientArea(geom, base, other);
        Vector3 gradH = GradientMeanCurvature(geom, base, other);
        Vector3 gradK = GradientGaussCurvature(geom, base, other);

        double energy = (*vertexEnergies)[base];

        if (clustering) {
            LWSCluster* c = clustering->clusters[clustering->vertexOwnership[base]];
            return energy * (c->a * gradA + c->b * gradH + c->c * gradK);
        }
        else {
            return energy * (a * gradA + b * gradH + c * gradK);
        }
    }

}