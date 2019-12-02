#include "obstacles/mesh_obstacle.h"
#include "tpe_energy_sc.h"

namespace LWS {
    MeshObstacle::MeshObstacle(std::shared_ptr<HalfedgeMesh> m, std::shared_ptr<VertexPositionGeometry> geom, double p_exp, double w) {
        mesh = m;
        geometry = geom;
        p = p_exp;
        weight = w;
        bvh = CreateBVHFromMesh(m, geom);
    }

    MeshObstacle::~MeshObstacle() {
        if (bvh) {
            delete bvh;
        }
    }

    void MeshObstacle::AddGradient(PolyCurveNetwork* curves, Eigen::MatrixXd &gradient) {
        int nVerts = curves->NumVertices();
        for (int i = 0; i < nVerts; i++) {
            Vector3 pos = curves->GetVertex(i)->Position();
            Vector3 force = AccumulateForce(bvh, pos);
            AddToRow(gradient, i, force * weight);
        }
    }

    double MeshObstacle::ComputeEnergy(PolyCurveNetwork* curves) {
        int nVerts = curves->NumVertices();
        double sumE = 0;
        for (int i = 0; i < nVerts; i++) {
            Vector3 pos = curves->GetVertex(i)->Position();
            sumE += AccumulateEnergy(bvh, pos);
        }
        return sumE;
    }

    inline double bodyEnergy(BVHNode3D* bvh, Vector3 point, double p) {
        Vector3 center = bvh->centerOfMass;
        double mass = bvh->totalMass;
        double distance = (center - point).norm();
        return 1.0 / pow(distance, p);
    }

    double MeshObstacle::AccumulateEnergy(BVHNode3D* node, Vector3 point) {
        if (node->IsEmpty()) {
            return 0;
        }
        else if (node->IsLeaf()) {
            return bodyEnergy(node, point, p);
        }
        else {
            if (node->shouldUseCell(point)) {
                return bodyEnergy(node, point, p);
            }
            else {
                double total = 0;
                for (BVHNode3D* child : node->children) {
                    total += AccumulateEnergy(child, point);
                }
                return total;
            }
        }
    }

    inline Vector3 bodyForce(BVHNode3D* bvh, Vector3 point, double p) {
        Vector3 center = bvh->centerOfMass;
        double mass = bvh->totalMass;
        
        Vector3 toPoint = center - point;
        double distance = toPoint.norm();
        toPoint /= distance;
        // Derivative of 1 / x^p = -p / x^(p+1)
        Vector3 grad_i = toPoint * p / pow(distance, p + 1);
        
        return mass * grad_i;
    }

    Vector3 MeshObstacle::AccumulateForce(BVHNode3D* node, Vector3 point) {
        if (node->IsEmpty()) {
            return Vector3{0, 0, 0};
        }
        else if (node->IsLeaf()) {
            return bodyForce(node, point, p);
        }
        else {
            if (node->shouldUseCell(point)) {
                return bodyForce(node, point, p);
            }
            else {
                Vector3 total{0, 0, 0};
                for (BVHNode3D* child : node->children) {
                    total += AccumulateForce(child, point);
                }
                return total;
            }
        }
    }


}