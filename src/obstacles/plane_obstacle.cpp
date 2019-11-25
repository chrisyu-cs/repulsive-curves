#include "obstacles/plane_obstacle.h"
#include "tpe_energy_sc.h"

namespace LWS {
    
    PlaneObstacle::PlaneObstacle(Vector3 c, Vector3 n, double p_exp) {
        center = c;
        normal = n.normalize();
        p = p_exp;
    }

    PlaneObstacle::~PlaneObstacle() {}

    void PlaneObstacle::AddGradient(PolyCurveNetwork* curves, Eigen::MatrixXd &gradient) {
        int nVerts = curves->NumVertices();

        for (int i = 0; i < nVerts; i++) {
            CurveVertex* v_i = curves->GetVertex(i);
            Vector3 p_i = v_i->Position();
            // Find the closest point on the plane
            Vector3 nearest = ClosestPoint(p_i);
            // Simulate an energy contribution of 1 / r^(b - a)
            Vector3 toPoint = nearest - p_i;
            double distance = (nearest - p_i).norm();
            toPoint /= distance;
            Vector3 grad_i = toPoint * 1.0 / pow(distance, p);

            AddToRow(gradient, v_i->GlobalIndex(), grad_i);
        }
    }
}