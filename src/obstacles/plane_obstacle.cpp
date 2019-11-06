#include "obstacles/plane_obstacle.h"
#include "tpe_energy_sc.h"

namespace LWS {
    
    PlaneObstacle::PlaneObstacle(Vector3 c, Vector3 n) {
        center = c;
        normal = n.normalize();
    }

    PlaneObstacle::~PlaneObstacle() {}

    void PlaneObstacle::AddGradient(PolyCurveNetwork* curves, Eigen::MatrixXd &gradient, double alpha, double beta) {
        int nVerts = curves->NumVertices();

        for (int i = 0; i < nVerts; i++) {
            CurveVertex* v_i = curves->GetVertex(i);
            Vector3 p_i = v_i->Position();
            Vector3 nearest = ClosestPoint(p_i);
            TangentMassPoint fake_y{Vector3{0, 0, 0}, 1, nearest, 0, 0};

            Vector3 grad_i = TPESC::tpe_grad(v_i, fake_y, alpha, beta, v_i);
            AddToRow(gradient, v_i->GlobalIndex(), grad_i);
        }
    }
}