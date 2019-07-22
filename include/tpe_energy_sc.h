#pragma once

#include "boundary_derivatives.h"
#include "vert_jacobian.h"
#include "poly_curve.h"

namespace LWS {

    enum class TPEPointType {
        Cluster, Point, Edge
    };

    struct TangentMassPoint {
        Vector3 tangent;
        double mass;
        Vector3 point;
        PointOnCurve curvePt;
        PointOnCurve curvePt2;

        TPEPointType type();
    };

    class TPESC {
        public:
        static double tpe_Kf(PointOnCurve i, PointOnCurve j, double alpha, double beta);
        static double tpe_Kf_pts(Vector3 p_x, Vector3 p_y, Vector3 tangent_x, double alpha, double beta);
        
        static double tpe_pair(PointOnCurve i, PointOnCurve j, double alpha, double beta);
        static double tpe_pair_pts(Vector3 p_x, Vector3 p_y, Vector3 tangent_x, double l_x, double l_y, double alpha, double beta);
        static double tpe_total(PolyCurve* loop, double alpha, double beta);
        static double tpe_total(PolyCurveGroup* curves, double alpha, double beta);

        static Vector3 tpe_grad_Kf(PointOnCurve i, PointOnCurve j, double alpha, double beta, PointOnCurve wrt);
        static Vector3 tpe_grad_Kf(TangentMassPoint i, PointOnCurve j, double alpha, double beta, PointOnCurve wrt);
        static Vector3 tpe_grad_Kf(PointOnCurve i, TangentMassPoint j, double alpha, double beta, PointOnCurve wrt);

        static Vector3 tpe_grad(PointOnCurve x, PointOnCurve y, double alpha, double beta, PointOnCurve wrt);
        static Vector3 tpe_grad(TangentMassPoint x, PointOnCurve y, double alpha, double beta, PointOnCurve wrt);
        static Vector3 tpe_grad(PointOnCurve x, TangentMassPoint y, double alpha, double beta, PointOnCurve wrt);

        // Projection of (f(i) - f(j)) onto the plane normal to the tangent at i
        static Vector3 proj_normal_plane(PointOnCurve i, PointOnCurve j);
        
        // Gradient of |proj_N|^alpha
        static Vector3 grad_norm_proj_alpha(PointOnCurve i, PointOnCurve j, double alpha, double beta, PointOnCurve wrt);
        static Vector3 grad_norm_proj_alpha(TangentMassPoint i, PointOnCurve j, double alpha, double beta, PointOnCurve wrt);
        static Vector3 grad_norm_proj_alpha(PointOnCurve i, TangentMassPoint j, double alpha, double beta, PointOnCurve wrt);
        static Vector3 grad_norm_proj_num(PointOnCurve i, PointOnCurve j, double alpha, double beta, PointOnCurve wrt, double h);

        // Jacobian of proj_N = ((f(x) - f(y)) - <f(x) - f(y), T> * T)
        static VertJacobian grad_tangent_proj(PointOnCurve i, PointOnCurve j, PointOnCurve wrt);
        static VertJacobian grad_tangent_proj(TangentMassPoint i, PointOnCurve j, PointOnCurve wrt);
        static VertJacobian grad_tangent_proj(PointOnCurve i, TangentMassPoint j, PointOnCurve wrt);
        static VertJacobian grad_tangent_proj_num(PointOnCurve i, PointOnCurve j, PointOnCurve wrt, double h);

        // Gradient of vertex dual length wrt a vertex
        static Vector3 length_wrt_vert(PointOnCurve lengthVert, PointOnCurve wrt);

        // Jacobian of the edge "tangent" (i.e. the edge vector) wrt a vertex
        static VertJacobian edge_tangent_wrt_vert(PointOnCurve tangentVert, PointOnCurve wrtVert);
        // Jacobian of the vertex tangent (average of surrounding edges) wrt a vertex
        static VertJacobian vertex_tangent_wrt_vert(PointOnCurve tangentVert, PointOnCurve wrtVert);
    };
}
