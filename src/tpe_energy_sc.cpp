#include "tpe_energy_sc.h"

namespace LWS {

    TPEPointType TangentMassPoint::type() {
        if (curvePt) {
            return TPEPointType::Point;
        }
        else {
            return TPEPointType::Cluster;
        }
    }


    void TPESC::FillGradientSingle(PolyCurveNetwork* curveNetwork, Eigen::MatrixXd &gradients, int i, int j, double alpha, double beta) {
        if (i == j) return;
        CurveVertex* i_pt = curveNetwork->GetVertex(i);
        CurveVertex* j_pt = curveNetwork->GetVertex(j);

        // Add i and neighbors of i
        std::vector<CurveVertex*> i_pts;
        i_pts.push_back(i_pt);
        for (int e = 0; e < i_pt->numEdges(); e++) {
            i_pts.push_back(i_pt->edge(e)->Opposite(i_pt));
        }

        // Add j and neighbors of j
        std::vector<CurveVertex*> j_pts;
        j_pts.push_back(j_pt);
        for (int e = 0; e < j_pt->numEdges(); e++) {
            j_pts.push_back(j_pt->edge(e)->Opposite(j_pt));
        }

        // Differentiate wrt neighbors of i
        for (CurveVertex* i_n : i_pts) {
            AddToRow(gradients, i_n->GlobalIndex(), TPESC::tpe_grad(i_pt, j_pt, alpha, beta, i_n));
        }
        // Differentiate wrt neighbors of j
        for (CurveVertex* j_n : j_pts) {
            bool noOverlap = true;
            // Only compute this derivative if j_n is not already included in one of the previous pairs
            for (CurveVertex* i_n : i_pts) {
                if (i_n == j_n) noOverlap = false;
            }
            if (noOverlap) {
                AddToRow(gradients, j_n->GlobalIndex(), TPESC::tpe_grad(i_pt, j_pt, alpha, beta, j_n));
            }
        }
    }

    void TPESC::FillGradientVectorDirect(PolyCurveNetwork* curveNetwork, Eigen::MatrixXd &gradients, double alpha, double beta) {
        int nVerts = curveNetwork->NumVertices();
        // Fill with zeros, so that the constraint entries are 0
        gradients.setZero();
        // Fill vertex entries with accumulated gradients
        for (int i = 0; i < nVerts; i++) {
            for (int j = 0; j < nVerts; j++) {
                if (i == j) continue;
                TPESC::FillGradientSingle(curveNetwork, gradients, i, j, alpha, beta);
            }
        }
    }

    double TPESC::tpe_Kf(CurveVertex* i, CurveVertex* j, double alpha, double beta) {
        if (i == j) return 0;
        else if (i->numEdges() > 2) return 0;

        Vector3 disp = i->Position() - j->Position();
        Vector3 T_i = i->Tangent();

        Vector3 normal_proj = disp - dot(disp, T_i) * T_i;
        double numer = pow(norm(normal_proj), alpha);
        double denom = pow(norm(disp), beta);
        return numer / denom;
    }

    double TPESC::tpe_pair(CurveVertex* i, CurveVertex* j, double alpha, double beta) {
        double kfxy = tpe_Kf(i, j, alpha, beta);
        double l_x = i->DualLength();
        double l_y = j->DualLength();
        return l_x * l_y * kfxy;
    }

    double TPESC::tpe_pair_pts(Vector3 p_x, Vector3 p_y, Vector3 tangent_x, double l_x, double l_y, double alpha, double beta) {
        double kfxy = tpe_Kf_pts(p_x, p_y, tangent_x, alpha, beta);
        return l_x * l_y * kfxy;
    }

    double TPESC::tpe_total(PolyCurveNetwork* curves, double alpha, double beta) {
        int nVerts = curves->NumVertices();
        double sumEnergy = 0;

        for (int i = 0; i < nVerts; i++) {
            for (int j = 0; j < nVerts; j++) {
                CurveVertex* pt_i = curves->GetVertex(i);
                CurveVertex* pt_j = curves->GetVertex(j);
                sumEnergy += tpe_pair(pt_i, pt_j, alpha, beta);
            }
        }
        return sumEnergy;
    }

    Vector3 TPESC::proj_normal_plane(CurveVertex* i, CurveVertex* j) {
        Vector3 disp = i->Position() - j->Position();
        Vector3 T_i = i->Tangent();
        return disp - dot(disp, T_i) * T_i;
    }
    
    Vector3 TPESC::tpe_grad_Kf(CurveVertex* i, CurveVertex* j, double alpha, double beta, CurveVertex* wrt) {
        if (i == j) return Vector3{0, 0, 0};
        else if (i->numEdges() > 2) return Vector3{0, 0, 0};

        // Get positions and displacement vectors
        Vector3 disp = i->Position() - j->Position();
        Vector3 T_i = i->Tangent();
        // Normalized displacement direction
        Vector3 unit_disp = disp;
        unit_disp = unit_disp.normalize();

        // Evaluate projection onto normal plane, v - <v,T> * T
        Vector3 normal_proj = disp - dot(disp, T_i) * T_i;
        // Numerator of energy is norm of projection ^ alpha
        double A = pow(norm(normal_proj), alpha);
        // Denominator of energy is distance between points ^ beta
        double B = pow(norm(disp), beta);

        // Derivative of numerator
        Vector3 deriv_A = grad_norm_proj_alpha(i, j, alpha, beta, wrt);
        // Derivative of denominator
        Vector3 deriv_B{0, 0, 0};
        if (wrt == i) {
            deriv_B = beta * pow(norm(disp), beta - 1) * unit_disp;
        }
        else if (wrt == j) {
            deriv_B = -beta * pow(norm(disp), beta - 1) * unit_disp;
        }
        // Quotient rule for A / B
        Vector3 total = (deriv_A * B - A * deriv_B) / (B * B);
        return total;
    }
    
    Vector3 TPESC::tpe_grad_Kf(TangentMassPoint i, CurveVertex* j, double alpha, double beta, CurveVertex* wrt) {
        // Get positions and displacement vectors
        Vector3 disp = i.point - j->Position();
        Vector3 T_i = i.tangent;
        // Normalized displacement direction
        Vector3 unit_disp = disp;
        unit_disp = unit_disp.normalize();

        // Evaluate projection onto normal plane, v - <v,T> * T
        Vector3 normal_proj = disp - dot(disp, T_i) * T_i;
        // Numerator of energy is norm of projection ^ alpha
        double A = pow(norm(normal_proj), alpha);
        // Denominator of energy is distance between points ^ beta
        double B = pow(norm(disp), beta);

        // Derivative of numerator
        Vector3 deriv_A = grad_norm_proj_alpha(i, j, alpha, beta, wrt);

        // Derivative of denominator
        Vector3 deriv_B{0, 0, 0};
        if (i.type() == TPEPointType::Point && wrt == i.curvePt) {
            deriv_B = beta * pow(norm(disp), beta - 1) * unit_disp;
        }
        else if (i.type() == TPEPointType::Edge && (wrt == i.curvePt || wrt == i.curvePt2)) {
            deriv_B = (beta * pow(norm(disp), beta - 1) * unit_disp) / 2;
        }
        if (wrt == j) {
            deriv_B += -beta * pow(norm(disp), beta - 1) * unit_disp;
        }
        // Quotient rule for A / B
        return (deriv_A * B - A * deriv_B) / (B * B);
    }
    
    Vector3 TPESC::tpe_grad_Kf(CurveVertex* i, TangentMassPoint j, double alpha, double beta, CurveVertex* wrt) {
        // Get positions and displacement vectors
        Vector3 disp = i->Position() - j.point;
        Vector3 T_i = i->Tangent();
        // Normalized displacement direction
        Vector3 unit_disp = disp;
        unit_disp = unit_disp.normalize();

        // Evaluate projection onto normal plane, v - <v,T> * T
        Vector3 normal_proj = disp - dot(disp, T_i) * T_i;
        // Numerator of energy is norm of projection ^ alpha
        double A = pow(norm(normal_proj), alpha);
        // Denominator of energy is distance between points ^ beta
        double B = pow(norm(disp), beta);

        // Derivative of numerator
        Vector3 deriv_A = grad_norm_proj_alpha(i, j, alpha, beta, wrt);
        // Derivative of denominator
        Vector3 deriv_B{0, 0, 0};
        if (wrt == i) {
            deriv_B = beta * pow(norm(disp), beta - 1) * unit_disp;
        }
        if (j.type() == TPEPointType::Point && wrt == j.curvePt) {
            deriv_B += -beta * pow(norm(disp), beta - 1) * unit_disp;
        }
        else if (j.type() == TPEPointType::Edge && (wrt == j.curvePt || wrt == j.curvePt2)) {
            deriv_B += (-beta * pow(norm(disp), beta - 1) * unit_disp) / 2;
        }
        // Quotient rule for A / B
        return (deriv_A * B - A * deriv_B) / (B * B);
    }

    Vector3 TPESC::tpe_grad(CurveVertex* x, CurveVertex* y, double alpha, double beta, CurveVertex* wrt) {
        // Computes the gradient of the kernel (K_f(x, y) dx dy) with respect to
        // the position of the vertex "wrt".
        if (x->numEdges() > 2) return Vector3{0, 0, 0};
        
        // First get the gradient of K_f(x, y)
        Vector3 grad_Kf = tpe_grad_Kf(x, y, alpha, beta, wrt);
        double Kf = tpe_Kf(x, y, alpha, beta);

        double l_x = x->DualLength();
        double l_y = y->DualLength();

        // d/dy of area(x)
        Vector3 grad_lx = length_wrt_vert(x, wrt);
        // d/dy of area(y)
        Vector3 grad_ly = length_wrt_vert(y, wrt);
        // Evaluate the product rule for dx*dy
        Vector3 prod_rule = grad_lx * l_y + l_x * grad_ly;
        // Evaluate the product rule for k dx dy
        Vector3 total = grad_Kf * l_x * l_y + Kf * prod_rule;
        return total;
    }

    Vector3 TPESC::tpe_grad(TangentMassPoint x, CurveVertex* y, double alpha, double beta, CurveVertex* wrt) {
        // Computes the gradient of the kernel (K_f(x, y) dx dy) with respect to
        // the position of the vertex "wrt".
        
        // First get the gradient of K_f(x, y)
        Vector3 grad_Kf = tpe_grad_Kf(x, y, alpha, beta, wrt);
        double Kf = tpe_Kf_pts(x.point, y->Position(), x.tangent, alpha, beta);

        double l_x = x.mass;
        double l_y = y->DualLength();

        // Area gradient for x depends on whether the mass point is distant
        Vector3 grad_lx{0, 0, 0};
        if (x.type() == TPEPointType::Point) {
            grad_lx = length_wrt_vert(x.curvePt, wrt);
        }
        else if (x.type() == TPEPointType::Edge) {
            if (wrt == x.curvePt) {
                grad_lx = (x.curvePt->Position() - x.curvePt2->Position());
                grad_lx = grad_lx.normalize();
            }
            else if (wrt == x.curvePt2) {
                grad_lx = (x.curvePt2->Position() - x.curvePt->Position());
                grad_lx = grad_lx.normalize();
            }
        }
        // d/dy of area(y)
        Vector3 grad_ly = length_wrt_vert(y, wrt);
        // Evaluate the product rule for dx*dy
        Vector3 prod_rule = grad_lx * l_y + l_x * grad_ly;

        // Evaluate the product rule for k dx dy
        return grad_Kf * l_x * l_y + Kf * prod_rule;
    }

    Vector3 TPESC::tpe_grad(CurveVertex* x, TangentMassPoint y, double alpha, double beta, CurveVertex* wrt) {
        // Here, y is a mass point (e.g. from Barnes-Hut), so gradients of y are assumed to be zero.
        // First get the gradient of K_f(x, y)
        Vector3 grad_Kf = tpe_grad_Kf(x, y, alpha, beta, wrt);
        double Kf = tpe_Kf_pts(x->Position(), y.point, x->Tangent(), alpha, beta);

        double l_x = x->DualLength();
        double l_y = y.mass;

        // d/dy of area(x)
        Vector3 grad_lx = length_wrt_vert(x, wrt);
        // Area gradient for y depends on whether the mass point is distant
        Vector3 grad_ly{0, 0, 0};
        if (y.type() == TPEPointType::Point) {
            grad_ly = length_wrt_vert(y.curvePt, wrt);
        }
        else if (y.type() == TPEPointType::Edge) {
            if (wrt == y.curvePt) {
                grad_ly = (y.curvePt->Position() - y.curvePt2->Position());
                grad_ly = grad_ly.normalize();
            }
            else if (wrt == y.curvePt2) {
                grad_ly = (y.curvePt2->Position() - y.curvePt->Position());
                grad_ly = grad_ly.normalize();
            }
        }

        // Evaluate the product rule for dx*dy
        Vector3 prod_rule = grad_lx * l_y + l_x * grad_ly;

        // Evaluate the product rule for k dx dy
        return grad_Kf * l_x * l_y + Kf * prod_rule;
    }

    Vector3 TPESC::grad_norm_proj_alpha(CurveVertex* i, CurveVertex* j, double alpha, double beta, CurveVertex* wrt) {
        Vector3 disp = i->Position() - j->Position();
        Vector3 T_i = i->Tangent();
        // Projection onto normal plane
        Vector3 normal_proj = disp - dot(disp, T_i) * T_i;
        double proj_len = norm(normal_proj);

        // If the displacement is actually exactly perpendicular to the tangent,
        // then the contribution is exactly 0.
        if (proj_len < 1e-12) return Vector3{0, 0, 0};

        // Derivative of |f(x) - ...|^alpha = alpha * |f(x) - ...|^(alpha - 1)
        double alpha_deriv = alpha * pow(proj_len, alpha - 1);
        // Normalized vector of projection onto normal plane
        Vector3 proj_normalized = normal_proj / proj_len;

        Vector3 zero{0, 0, 0};
        VertJacobian deriv_disp{zero, zero, zero};
        if (wrt == i) {
            deriv_disp.directional_x = Vector3{1, 0, 0};
            deriv_disp.directional_y = Vector3{0, 1, 0};
            deriv_disp.directional_z = Vector3{0, 0, 1};
        }
        else if (wrt == j) {
            deriv_disp.directional_x = Vector3{-1, 0, 0};
            deriv_disp.directional_y = Vector3{0, -1, 0};
            deriv_disp.directional_z = Vector3{0, 0, -1};
        }

        // Derivative of <f(x) - f(y), T> * T
        VertJacobian deriv_T_inner = grad_tangent_proj(i, j, wrt);
        VertJacobian deriv_N_proj = deriv_disp - deriv_T_inner;

        Vector3 total = alpha_deriv * deriv_N_proj.LeftMultiply(proj_normalized);
        return total;
    }

    Vector3 TPESC::grad_norm_proj_alpha(TangentMassPoint i, CurveVertex* j, double alpha, double beta, CurveVertex* wrt) {
        Vector3 disp = i.point - j->Position();
        Vector3 T_i = i.tangent;
        // Projection onto normal plane
        Vector3 normal_proj = disp - dot(disp, T_i) * T_i;
        double proj_len = norm(normal_proj);

        // Derivative of |f(x) - ...|^alpha = alpha * |f(x) - ...|^(alpha - 1)
        double alpha_deriv = alpha * pow(proj_len, alpha - 1);
        // Normalized vector of projection onto normal plane
        Vector3 proj_normalized = normal_proj / proj_len;

        Vector3 zero{0, 0, 0};
        VertJacobian deriv_disp{zero, zero, zero};
        if (i.type() == TPEPointType::Point && wrt == i.curvePt) {
            deriv_disp.directional_x = Vector3{1, 0, 0};
            deriv_disp.directional_y = Vector3{0, 1, 0};
            deriv_disp.directional_z = Vector3{0, 0, 1};
        }
        if (i.type() == TPEPointType::Edge && (wrt == i.curvePt || wrt == i.curvePt2)) {
            deriv_disp.directional_x = Vector3{0.5, 0, 0};
            deriv_disp.directional_y = Vector3{0, 0.5, 0};
            deriv_disp.directional_z = Vector3{0, 0, 0.5};
        }
        else if (wrt == j) {
            deriv_disp.directional_x = Vector3{-1, 0, 0};
            deriv_disp.directional_y = Vector3{0, -1, 0};
            deriv_disp.directional_z = Vector3{0, 0, -1};
        }

        // Derivative of <f(x) - f(y), T> * T
        VertJacobian deriv_T_inner = grad_tangent_proj(i, j, wrt);
        VertJacobian deriv_N_proj = deriv_disp - deriv_T_inner;

        return alpha_deriv * deriv_N_proj.LeftMultiply(proj_normalized);
    }

    Vector3 TPESC::grad_norm_proj_alpha(CurveVertex* i, TangentMassPoint j, double alpha, double beta, CurveVertex* wrt) {
        Vector3 disp = i->Position() - j.point;
        Vector3 T_i = i->Tangent();
        // Projection onto normal plane
        Vector3 normal_proj = disp - dot(disp, T_i) * T_i;
        double proj_len = norm(normal_proj);

        // Derivative of |f(x) - ...|^alpha = alpha * |f(x) - ...|^(alpha - 1)
        double alpha_deriv = alpha * pow(proj_len, alpha - 1);
        // Normalized vector of projection onto normal plane
        Vector3 proj_normalized = normal_proj / proj_len;

        Vector3 zero{0, 0, 0};
        VertJacobian deriv_disp{zero, zero, zero};
        if (wrt == i) {
            deriv_disp.directional_x = Vector3{1, 0, 0};
            deriv_disp.directional_y = Vector3{0, 1, 0};
            deriv_disp.directional_z = Vector3{0, 0, 1};
        }
        else if (j.type() == TPEPointType::Point && wrt == j.curvePt) {
            deriv_disp.directional_x = Vector3{-1, 0, 0};
            deriv_disp.directional_y = Vector3{0, -1, 0};
            deriv_disp.directional_z = Vector3{0, 0, -1};
        }
        else if (j.type() == TPEPointType::Edge && (wrt == j.curvePt || wrt == j.curvePt2)) {
            deriv_disp.directional_x = Vector3{-0.5, 0, 0};
            deriv_disp.directional_y = Vector3{0, -0.5, 0};
            deriv_disp.directional_z = Vector3{0, 0, -0.5};
        }

        // Derivative of <f(x) - f(y), T> * T
        VertJacobian deriv_T_inner = grad_tangent_proj(i, j, wrt);
        VertJacobian deriv_N_proj = deriv_disp - deriv_T_inner;

        return alpha_deriv * deriv_N_proj.LeftMultiply(proj_normalized);
    }

    Vector3 TPESC::grad_norm_proj_num(CurveVertex* i, CurveVertex* j, double alpha, double beta, CurveVertex* wrt, double h) {
        Vector3 origPos = wrt->Position();
        double orig = pow(norm(proj_normal_plane(i, j)), alpha);

        wrt->SetPosition(origPos + Vector3{h, 0, 0});
        double xVal = pow(norm(proj_normal_plane(i, j)), alpha);

        wrt->SetPosition(origPos + Vector3{0, h, 0});
        double yVal = pow(norm(proj_normal_plane(i, j)), alpha);

        wrt->SetPosition(origPos + Vector3{0, 0, h});
        double zVal = pow(norm(proj_normal_plane(i, j)), alpha);

        wrt->SetPosition(origPos);

        double xDeriv = (xVal - orig) / h;
        double yDeriv = (yVal - orig) / h;
        double zDeriv = (zVal - orig) / h;

        return Vector3{xDeriv, yDeriv, zDeriv};
    }

    VertJacobian TPESC::grad_tangent_proj(CurveVertex* i, CurveVertex* j, CurveVertex* wrt) {
        // Differentiate the inner product
        Vector3 disp = i->Position() - j->Position();
        Vector3 T_i = i->Tangent();
        double disp_dot_T = dot(disp, T_i);

        Vector3 inner_deriv_A_B{0, 0, 0};
        if (wrt == i) {
            inner_deriv_A_B = T_i;
        }
        else if (wrt == j) {
            inner_deriv_A_B = -T_i;
        }
        VertJacobian deriv_T = vertex_tangent_wrt_vert(i, wrt);
        Vector3 inner_A_deriv_B = deriv_T.LeftMultiply(disp);
        Vector3 deriv_inner = inner_deriv_A_B + inner_A_deriv_B;

        // Now use product rule for <f(x) - f(y), T)> * T
        VertJacobian deriv_A_B = outer_product_to_jacobian(T_i, deriv_inner);
        VertJacobian A_deriv_B = disp_dot_T * deriv_T;

        return deriv_A_B + A_deriv_B;
    }

    VertJacobian TPESC::grad_tangent_proj(TangentMassPoint i, CurveVertex* j, CurveVertex* wrt) {
        // Differentiate the inner product
        Vector3 disp = i.point - j->Position();
        Vector3 T_i = i.tangent;
        double disp_dot_T = dot(disp, T_i);

        Vector3 inner_deriv_A_B{0, 0, 0};
        if (wrt == j) {
            inner_deriv_A_B = -T_i;
        }
        VertJacobian deriv_T = vertex_tangent_wrt_vert(i.curvePt, wrt);
        Vector3 inner_A_deriv_B = deriv_T.LeftMultiply(disp);
        Vector3 deriv_inner = inner_deriv_A_B + inner_A_deriv_B;

        // Now use product rule for <f(x) - f(y), T)> * T
        VertJacobian deriv_A_B = outer_product_to_jacobian(T_i, deriv_inner);
        VertJacobian A_deriv_B = disp_dot_T * deriv_T;

        return deriv_A_B + A_deriv_B;
    }

    VertJacobian TPESC::grad_tangent_proj(CurveVertex* i, TangentMassPoint j, CurveVertex* wrt) {
        // Differentiate the inner product
        Vector3 disp = i->Position() - j.point;
        Vector3 T_i = i->Tangent();
        double disp_dot_T = dot(disp, T_i);

        Vector3 inner_deriv_A_B{0, 0, 0};
        if (wrt == i) {
            inner_deriv_A_B = T_i;
        }
        VertJacobian deriv_T = vertex_tangent_wrt_vert(i, wrt);
        Vector3 inner_A_deriv_B = deriv_T.LeftMultiply(disp);
        Vector3 deriv_inner = inner_deriv_A_B + inner_A_deriv_B;

        // Now use product rule for <f(x) - f(y), T)> * T
        VertJacobian deriv_A_B = outer_product_to_jacobian(T_i, deriv_inner);
        VertJacobian A_deriv_B = disp_dot_T * deriv_T;

        return deriv_A_B + A_deriv_B;
    }

    VertJacobian TPESC::grad_tangent_proj_num(CurveVertex* i, CurveVertex* j, CurveVertex* wrt, double h) {
        Vector3 origPos = wrt->Position();
        Vector3 origTangent = dot(i->Position() - j->Position(), i->Tangent()) * i->Tangent();

        wrt->SetPosition(origPos + Vector3{h, 0, 0});
        Vector3 xTangent = dot(i->Position() - j->Position(), i->Tangent()) * i->Tangent();

        wrt->SetPosition(origPos + Vector3{0, h, 0});
        Vector3 yTangent = dot(i->Position() - j->Position(), i->Tangent()) * i->Tangent();

        wrt->SetPosition(origPos + Vector3{0, 0, h});
        Vector3 zTangent = dot(i->Position() - j->Position(), i->Tangent()) * i->Tangent();

        wrt->SetPosition(origPos);

        Vector3 xDeriv = (xTangent - origTangent) / h;
        Vector3 yDeriv = (yTangent - origTangent) / h;
        Vector3 zDeriv = (zTangent - origTangent) / h;

        return VertJacobian{xDeriv, yDeriv, zDeriv};
    }
    

    VertJacobian TPESC::edge_tangent_wrt_vert(CurveEdge* edge, CurveVertex* wrtVert) {
        // get positions
        CurveVertex* prevVert = edge->prevVert;
        CurveVertex* nextVert = edge->nextVert;
        Vector3 v_h = prevVert->Position();
        Vector3 v_i = nextVert->Position();

        if (wrtVert != prevVert && wrtVert != nextVert) {
            Vector3 zero{0, 0, 0};
            return VertJacobian{zero, zero, zero};
        }

        Vector3 v_tangent = v_i - v_h;
        double v_norm = norm(v_tangent);
        Vector3 v_normalized = v_tangent / v_norm;

        VertJacobian I{Vector3{1, 0, 0}, Vector3{0, 1, 0}, Vector3{0, 0, 1}};

        VertJacobian deriv_A_B = I * v_norm;
        VertJacobian A_deriv_B = outer_product_to_jacobian(v_tangent, v_normalized);
        VertJacobian deriv = (deriv_A_B - A_deriv_B) * (1.0 / (v_norm * v_norm));

        // If we're differentiating the tail vertex, the derivative is negative
        if (wrtVert == prevVert) return -1 * deriv;
        // Otherwise we're differentiating the head vertex, and the derivative is positive
        else return deriv;
    }

    VertJacobian TPESC::vertex_tangent_wrt_vert(CurveVertex* tangentVert, CurveVertex* wrtVert) {
        if (!tangentVert || !wrtVert) {
            return VertJacobian{Vector3{0, 0, 0}, Vector3{0, 0, 0}, Vector3{0, 0, 0}};
        }

        if (tangentVert->numEdges() != 2) {
            if (tangentVert->numEdges() == 1) {
                CurveEdge* edge = tangentVert->edge(0);
                Vector3 tangent = edge->Tangent();
                // Derivative of T
                return edge_tangent_wrt_vert(edge, wrtVert);
            }

            else {
                return VertJacobian{Vector3{0, 0, 0}, Vector3{0, 0, 0}, Vector3{0, 0, 0}};
            }
        }

        CurveEdge* prevEdge = tangentVert->edge(0);
        CurveEdge* nextEdge = tangentVert->edge(1);

        Vector3 prevTangent = prevEdge->Tangent();
        Vector3 nextTangent = nextEdge->Tangent();

        Vector3 sumTangents = prevTangent + nextTangent;
        double normSum = norm(sumTangents);
        Vector3 vertTangent = sumTangents;
        vertTangent = vertTangent.normalize();

        // Quotient rule on (T1 + T2) / |T1 + T2| 
        VertJacobian derivSumTs = edge_tangent_wrt_vert(prevEdge, wrtVert)
            + edge_tangent_wrt_vert(nextEdge, wrtVert);

        Vector3 derivNorm = derivSumTs.LeftMultiply(vertTangent);
        VertJacobian deriv_A_B = derivSumTs * normSum;
        VertJacobian A_deriv_B = outer_product_to_jacobian(sumTangents, derivNorm);

        return (deriv_A_B - A_deriv_B) * (1.0 / (normSum * normSum)); 
    }


    Vector3 TPESC::length_wrt_vert(CurveVertex* lengthVert, CurveVertex* wrt) {
        // If differentiating wrt self, need to consider both side
        if (lengthVert == wrt) {
            Vector3 sumDirections{0, 0, 0};
            Vector3 center = lengthVert->Position();
            for (int e = 0; e < lengthVert->numEdges(); e++) {
                Vector3 other = lengthVert->edge(e)->Opposite(lengthVert)->Position();
                Vector3 outward = other - center;
                outward = outward.normalize();
                sumDirections += outward;
            }
            // Gradient is half the sum of unit vectors along outgoing edges
            return -0.5 * sumDirections;
        }
        // Otherwise, only consider the one edge between other and vert
        else {
            for (int e = 0; e < lengthVert->numEdges(); e++) {
                CurveVertex* other = lengthVert->edge(e)->Opposite(lengthVert);
                if (other == wrt) {
                    Vector3 outward = other->Position() - lengthVert->Position();
                    outward = outward.normalize();
                    return outward / 2;
                }
            }
            return Vector3{0, 0, 0};
        }
    }

}