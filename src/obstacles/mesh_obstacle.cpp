#include "obstacles/mesh_obstacle.h"
#include "tpe_energy_sc.h"

namespace LWS {
    MeshObstacle::MeshObstacle(std::shared_ptr<HalfedgeMesh> m, std::shared_ptr<VertexPositionGeometry> geom) {
        mesh = m;
        geometry = geom;
        bvh = CreateBVHFromMesh(m, geom);
    }

    MeshObstacle::~MeshObstacle() {
        if (bvh) {
            delete bvh;
        }
    }

    void MeshObstacle::AddGradient(PolyCurveNetwork* curves, Eigen::MatrixXd &gradient, double alpha, double beta) {
        int nVerts = curves->NumVertices();
        for (int i = 0; i < nVerts; i++) {
            Vector3 pos = curves->GetVertex(i)->Position();
            Vector3 force = AccumulateForce(bvh, pos, alpha, beta);
            AddToRow(gradient, i, force);
        }
    }

    inline Vector3 bodyForce(BVHNode3D* bvh, Vector3 point, double alpha, double beta) {
        Vector3 center = bvh->centerOfMass;
        double mass = bvh->totalMass;
        
        Vector3 toPoint = center - point;
        double distance = toPoint.norm();
        toPoint /= distance;
        Vector3 grad_i = toPoint * 1.0 / pow(distance, beta - alpha + 1);
        
        return mass * grad_i;
    }

    Vector3 MeshObstacle::AccumulateForce(BVHNode3D* node, Vector3 point, double alpha, double beta) {
        if (node->IsEmpty()) {
            return Vector3{0, 0, 0};
        }
        else if (node->IsLeaf()) {
            return bodyForce(node, point, alpha, beta);
        }
        else {
            if (node->shouldUseCell(point)) {
                return bodyForce(node, point, alpha, beta);
            }
            else {
                Vector3 total{0, 0, 0};
                for (BVHNode3D* child : node->children) {
                    total += AccumulateForce(child, point, alpha, beta);
                }
                return total;
            }
        }
    }


}