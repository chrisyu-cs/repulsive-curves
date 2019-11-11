#pragma once

#include "obstacles/obstacle.h"
#include "geometrycentral/surface/vertex_position_geometry.h"
#include "geometrycentral/surface/halfedge_mesh.h"
#include "spatial/tpe_bvh.h"

namespace LWS {
    using namespace geometrycentral;
    using namespace surface;

    class MeshObstacle : public Obstacle {
        public:
        std::shared_ptr<HalfedgeMesh> mesh;
        std::shared_ptr<VertexPositionGeometry> geometry;

        MeshObstacle(std::shared_ptr<HalfedgeMesh> m, std::shared_ptr<VertexPositionGeometry> geom);
        virtual ~MeshObstacle();
        virtual void AddGradient(PolyCurveNetwork* curves, Eigen::MatrixXd &gradient, double alpha, double beta);

        private:
        BVHNode3D* bvh;
        Vector3 AccumulateForce(BVHNode3D* node, Vector3 point, double alpha, double beta);
    };
}
