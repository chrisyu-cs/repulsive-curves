#include "obstacles/sphere_obstacle.h"
#include "tpe_energy_sc.h"

namespace LWS {
    
    SphereObstacle::SphereObstacle(Vector3 c, double r, double p_exp) {
        center = c;
        radius = r;
        p = p_exp;
    }

    SphereObstacle::~SphereObstacle() {}

    double SphereObstacle::ComputeEnergy(PolyCurveNetwork* curves) {
        int nVerts = curves->NumVertices();
        double sumE = 0;

        for (int i = 0; i < nVerts; i++) {
            CurveVertex* v_i = curves->GetVertex(i);
            Vector3 p_i = v_i->Position();
            // If we're very close to the center of the sphere, gradient is 0
            if ((p_i - center).norm() < 1e-6) continue;
            // Find the closest point on the plane
            Vector3 nearest = ClosestPoint(p_i);
            // Simulate an energy contribution of 1 / r^(b - a)
            Vector3 toPoint = nearest - p_i;
            double distance = (nearest - p_i).norm();
            sumE += 1.0 / pow(distance, p);
        }

        return sumE;
    }

    void SphereObstacle::AddGradient(PolyCurveNetwork* curves, Eigen::MatrixXd &gradient) {
        int nVerts = curves->NumVertices();

        for (int i = 0; i < nVerts; i++) {
            CurveVertex* v_i = curves->GetVertex(i);
            Vector3 p_i = v_i->Position();
            // If we're very close to the center of the sphere, gradient is 0
            if ((p_i - center).norm() < 1e-6) continue;
            // Find the closest point on the plane
            Vector3 nearest = ClosestPoint(p_i);
            // Simulate an energy contribution of 1 / r^(b - a)
            Vector3 toPoint = nearest - p_i;
            double distance = (nearest - p_i).norm();
            toPoint /= distance;
            Vector3 grad_i = toPoint * p / pow(distance, p + 1);

            // Add an antipodal component
            // double oppDistance = 2 * radius - distance;
            // grad_i += -toPoint * 1.0 / pow(oppDistance, p); 

            AddToRow(gradient, v_i->GlobalIndex(), grad_i);
        }
    }
}
