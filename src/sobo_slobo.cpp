#include "sobo_slobo.h"

namespace LWS {
    void SobolevCurves::LocalMatrix(EdgePositionPair e1, EdgePositionPair e2, double t1, double t2,
        double alpha, double beta, double out[4][4]) {
        
        Vector3 f1 = (1 - t1) * e1.f1 + t1 * e1.f2;
        Vector3 f2 = (1 - t2) * e2.f1 + t2 * e2.f2;

        // Edge lengths
        double l1 = norm(e1.f1 - e1.f2);
        double l2 = norm(e2.f1 - e2.f2);
        double area = l1 * l2;

        // Squared distance between the interpolated points
        double r = norm(f1 - f2);

        // Assign the local matrix
        out[0][0] = 1.0 / (l1 * l1);
        out[0][1] = -1.0 / (l1 * l1);
        out[0][2] = -1.0 / (l1 * l2);
        out[0][3] = 1.0 / (l1 * l2);

        out[1][0] = -1.0 / (l1 * l1);
        out[1][1] = 1.0 / (l1 * l1);
        out[1][2] = 1.0 / (l1 * l2);
        out[1][3] = -1.0 / (l1 * l2);
        
        out[2][0] = -1.0 / (l1 * l2);
        out[2][1] = 1.0 / (l1 * l2);
        out[2][2] = 1.0 / (l2 * l2);
        out[2][3] = -1.0 / (l2 * l2);
        
        out[3][0] = 1.0 / (l1 * l2);
        out[3][1] = -1.0 / (l1 * l2);
        out[3][2] = -1.0 / (l2 * l2);
        out[3][3] = 1.0 / (l2 * l2);
        
        //double denominator = pow(r2, (beta - alpha) / 2);
        // This exponent comes from http://brickisland.net/winegarden/2018/12/10/sobolev-slobodeckij-spaces/
        // double pow_s = 1. - 1. / alpha;
        double pow_s = (beta - alpha);
        double denominator = pow(r, pow_s);

        for (int i = 0; i < 4; i++) {
            for (int j = 0; j < 4; j++) {
                out[i][j] /= denominator;
                out[i][j] *= area;
            }
        }
    }

    void SobolevCurves::KfMatrix(EdgePositionPair e1, EdgePositionPair e2, double t1, double t2,
        double alpha, double beta, double out[4][4]) {

        // Interpolated position on edges
        Vector3 f1 = (1 - t1) * e1.f1 + t1 * e1.f2;
        Vector3 f2 = (1 - t2) * e2.f1 + t2 * e2.f2;
        // Tangent of edge 1
        Vector3 tangent_x = e1.f2 - e1.f1;
        // Edge lengths
        double l1 = norm(e1.f1 - e1.f2);
        double l2 = norm(e2.f1 - e2.f2);
        double area = l1 * l2;

        // Squared distance between the interpolated points
        double r2 = norm2(f1 - f2);

        // Assign the local matrix
        out[0][0] = (1 - t1) * (1 - t1);
        out[0][1] = (1 - t1) * t1;
        out[0][2] = (1 - t1) * (-1 + t2);
        out[0][3] = -(1 - t1) * t2;

        out[1][0] = (1 - t1) * t1;
        out[1][1] = t1 * t1;
        out[1][2] = t1 * (-1 + t2);
        out[1][3] = -t1 * t2;
        
        out[2][0] = (1 - t1) * (-1 + t2);
        out[2][1] = t1 * (-1 + t2);
        out[2][2] = (-1 + t2) * (-1 + t2);
        out[2][3] = -(-1 + t2) * t2;
        
        out[3][0] = -(1 - t1) * t2;
        out[3][1] = -t1 * t2;
        out[3][2] = -(-1 + t2) * t2;
        out[3][3] = t2 * t2;

        double denominator = sqrt(r2);
        double kfxy = TPESC::tpe_Kf_pts(f1, f2, tangent_x, alpha, beta);

        for (int i = 0; i < 4; i++) {
            for (int j = 0; j < 4; j++) {
                out[i][j] = (area * kfxy * out[i][j]) / denominator;
            }
        }
    }

    void SobolevCurves::AddToFirst(double acc[4][4], double add[4][4], double weight) {
        for (int i = 0; i < 4; i++){
            for (int j = 0; j < 4; j++) {
                acc[i][j] += add[i][j] * weight;
            }
        }
    }

    void SobolevCurves::IntegrateLocalMatrix(EdgePositionPair e1, EdgePositionPair e2,
        double alpha, double beta, double out[4][4]) {

        // Quadrature points and weights
        double points[1] = {0.5};
        double weights[1] = {1};
        //double points[3] = {0.1127015, 0.5, 0.8872985};
        //double weights[3] = {0.2777778, 0.4444444, 0.2777778};
        //double points[3] = {0.1666667, 0.5, 0.8333333};
        //double weights[3] = {0.3333333, 0.3333333, 0.3333333};
        double temp_out[4][4];

        // Zero out the output array
        for (int i = 0; i < 4; i++) {
            for (int j = 0; j < 4; j++) {
                out[i][j] = 0;
            }
        }

        // Loop over all pairs of quadrature points
        for (int x = 0; x < 1; x++) {
            for (int y = 0; y < 1; y++) {
                double weight = weights[x] * weights[y];
                // Compute contribution to local matrix from each quadrature point
                LocalMatrix(e1, e2, points[x], points[y], alpha, beta, temp_out);
                // Add to result matrix, weighted by quadrature weight
                AddToFirst(out, temp_out, weight);
                // Compute contribution to local matrix from each quadrature point
                // KfMatrix(e1, e2, points[x], points[y], alpha, beta, temp_out);
                // AddToFirst(out, temp_out, weight);
            }
        } 
    }

    Vector3 HatGradientOnEdge(PointOnCurve edgeStart, PointOnCurve vertex) {
        PointOnCurve edgeEnd = edgeStart.Next();

        if (vertex == edgeStart) {
            Vector3 towardsVertex = edgeStart.Position() - edgeEnd.Position();
            double length = norm(towardsVertex);
            // 1 over length times normalized edge vector towards the vertex
            return towardsVertex / (length * length);
        }
        else if (vertex == edgeEnd) {
            Vector3 towardsVertex = edgeEnd.Position() - edgeStart.Position();
            double length = norm(towardsVertex);
            return towardsVertex / (length * length);
        }
        else {
            return Vector3{0, 0, 0};
        }
    }

    void SobolevCurves::AddEdgePairContribution(PolyCurveGroup* loop, double alpha, double beta,
    PointOnCurve s, PointOnCurve t, Eigen::MatrixXd &A) {

        PointOnCurve endpoints[4] = {s, s.Next(), t, t.Next()};
        double len1 = norm(s.Position() - s.Next().Position());
        double len2 = norm(t.Position() - t.Next().Position());
        Vector3 mid1 = (s.Position() + s.Next().Position()) / 2;
        Vector3 mid2 = (t.Position() + t.Next().Position()) / 2;
        double denom = pow(norm(mid1 - mid2), beta - alpha);

        for (PointOnCurve u : endpoints) {
            for (PointOnCurve v : endpoints) {
                Vector3 u_hat_s = HatGradientOnEdge(s, u);
                Vector3 u_hat_t = HatGradientOnEdge(t, u);
                Vector3 v_hat_s = HatGradientOnEdge(s, v);
                Vector3 v_hat_t = HatGradientOnEdge(t, v);

                double numer = dot(u_hat_s - u_hat_t, v_hat_s - v_hat_t);
                int index_u = loop->GlobalIndex(u);
                int index_v = loop->GlobalIndex(v);

                A(index_u, index_v) += (numer / denom) * len1 * len2;
            }
        }
    }

    double HatMidpoint(PointOnCurve edgeStart, PointOnCurve vertex) {
        PointOnCurve edgeEnd = edgeStart.Next();

        if (vertex == edgeStart) {
            return 0.5;
        }
        else if (vertex == edgeEnd) {
            return 0.5;
        }
        else {
            return 0;
        }
    }

    void SobolevCurves::AddEdgePairContributionLow(PolyCurveGroup* loop, double alpha, double beta,
    PointOnCurve s, PointOnCurve t, Eigen::MatrixXd &A) {

        PointOnCurve endpoints[4] = {s, s.Next(), t, t.Next()};

        double len1 = norm(s.Position() - s.Next().Position());
        double len2 = norm(t.Position() - t.Next().Position());
        Vector3 mid_s = (s.Position() + s.Next().Position()) / 2;
        Vector3 mid_t = (t.Position() + t.Next().Position()) / 2;
        double denom = pow(norm(mid_s - mid_t), 2);

        Vector3 tangent_s = (s.Next().Position() - s.Position()).normalize();

        double kf_st = TPESC::tpe_Kf_pts(mid_s, mid_t, tangent_s, alpha, beta);

        for (PointOnCurve u : endpoints) {
            Vector3 tangent_u = (u.Next().Position() - u.Position()).normalize();
            
            for (PointOnCurve v : endpoints) {
                double u_s = HatMidpoint(s, u);
                double u_t = HatMidpoint(t, u);
                double v_s = HatMidpoint(s, v);
                double v_t = HatMidpoint(t, v);

                double numer = (u_s - u_t) * (v_s - v_t);
                int index_u = loop->GlobalIndex(u);
                int index_v = loop->GlobalIndex(v);

                A(index_u, index_v) += (numer / denom) * kf_st * len1 * len2;
            }
        }
    }

    void SobolevCurves::FillGlobalMatrix(PolyCurveGroup* curves, double alpha, double beta, Eigen::MatrixXd &A) {
        double out[4][4];
        int nVerts = curves->NumVertices();

        int numTerms = 0;

        for (int i = 0; i < nVerts; i++) {
            PointOnCurve pc_i = curves->GetCurvePoint(i);
            
            for (int j = 0; j < nVerts; j++) {
                PointOnCurve pc_j = curves->GetCurvePoint(j);
                // if (pc_i == pc_j) continue;
                if (pc_i == pc_j || pc_i.Next() == pc_j || pc_i == pc_j.Next() || pc_i.Next() == pc_j.Next()) continue;

                AddEdgePairContribution(curves, alpha, beta, pc_i, pc_j, A);
                // AddEdgePairContributionLow(curves, alpha, beta, pc_i, pc_j, A);
            }
        }
    }

    double SobolevCurves::SobolevDot(std::vector<Vector3> &as, std::vector<Vector3> &bs, Eigen::MatrixXd &J) {
        Eigen::VectorXd A, B;
        A.setZero(as.size());
        B.setZero(bs.size());

        // X component
        for (size_t i = 0; i < as.size(); i++) {
            A(i) = as[i].x;
            B(i) = bs[i].x;
        }
        double dot_x = A.transpose() * J * B;

        // Y component
        for (size_t i = 0; i < as.size(); i++) {
            A(i) = as[i].y;
            B(i) = bs[i].y;
        }
        double dot_y = A.transpose() * J * B;

        // Z component
        for (size_t i = 0; i < as.size(); i++) {
            A(i) = as[i].z;
            B(i) = bs[i].z;
        }
        double dot_z = A.transpose() * J * B;

        return dot_x + dot_y + dot_z;
    }

    void SobolevCurves::SobolevNormalize(std::vector<Vector3> &as, Eigen::MatrixXd &J) {
        double sobo_norm = sqrt(SobolevDot(as, as, J));
        for (size_t i = 0; i < as.size(); i++) {
            as[i] /= sobo_norm;
        }
    }

    void SobolevCurves::SobolevOrthoProjection(std::vector<Vector3> &as, std::vector<Vector3> &bs, Eigen::MatrixXd &J) {
        double d1 = SobolevDot(as, bs, J);
        double d2 = SobolevDot(bs, bs, J);
        // Take the projection, and divide out the norm of b
        for (size_t i = 0; i < as.size(); i++) {
            as[i] -= (d1 / d2) * bs[i];
        }
    }

    void SobolevCurves::DfMatrix(PolyCurveGroup* loop, Eigen::MatrixXd &out) {
        int nVerts = loop->NumVertices();
        out.setZero(nVerts * 3, nVerts);
        
        for (int i = 0; i < nVerts; i++) {
            PointOnCurve p1 = loop->GetCurvePoint(i);
            int i1 = loop->GlobalIndex(p1);
            PointOnCurve p2 = p1.Next();
            int i2 = loop->GlobalIndex(p2);

            // Positive weight on the second vertex, negative on the first
            Vector3 e12 = (p2.Position() - p1.Position());
            double length = norm(e12);

            Vector3 weight = e12 / (length * length);
            out(3 * i1, i1) = -weight.x;
            out(3 * i1 + 1, i1) = -weight.y;
            out(3 * i1 + 2, i1) = -weight.z;

            out(3 * i1, i2) = weight.x;
            out(3 * i1 + 1, i2) = weight.y;
            out(3 * i1 + 2, i2) = weight.z;
        }
    }

    void SobolevCurves::FillAfMatrix(PolyCurveGroup* loop, double alpha, double beta, Eigen::MatrixXd &A) {
        int nVerts = loop->NumVertices();
        double pow_s = (beta - alpha);

        for (int i = 0; i < nVerts; i++) {
            PointOnCurve p_i = loop->GetCurvePoint(i);
            Vector3 mid_i = (p_i.Position() + p_i.Next().Position()) / 2;
            for (int j = 0; j < nVerts; j++) {
                PointOnCurve p_j = loop->GetCurvePoint(j);
                if (p_i == p_j || p_i.Next() == p_j || p_i == p_j.Next() || p_i.Next() == p_j.Next()) continue;
                Vector3 mid_j = (p_j.Position() + p_j.Next().Position()) / 2;

                A(i, j) = (p_i.DualLength() * p_j.DualLength()) / pow(norm(mid_i - mid_j), pow_s);
            }
        }
    }

    void SobolevCurves::ApplyDf(PolyCurveGroup* loop, Eigen::VectorXd &as, Eigen::MatrixXd &out) {
        for (PolyCurve *c : loop->curves) {
            int nVerts = c->NumVertices();
            for (int i = 0; i < nVerts; i++) {
                PointOnCurve p1 = loop->GetCurvePoint(i);
                int i1 = loop->GlobalIndex(p1);
                PointOnCurve p2 = p1.Next();
                int i2 = loop->GlobalIndex(p2);

                // Positive weight on the second vertex, negative on the first
                double d12 = as(i2) - as(i1);
                Vector3 e12 = (p2.Position() - p1.Position());
                double length = norm(e12);

                Vector3 grad12 = (d12 * e12) / (length * length);

                SetRow(out, i1, grad12);
            }
        }
    } 

    void SobolevCurves::ApplyDfTranspose(PolyCurveGroup* loop, Eigen::MatrixXd &es, Eigen::VectorXd &out) {
        for (PolyCurve *c : loop->curves) {
            int nVerts = c->NumVertices();
            for (int i = 0; i < nVerts; i++) {
                // This is the first vertex of the next edge, so negative weight
                PointOnCurve p1 = loop->GetCurvePoint(i);
                // This is the second vertex of the previous edge, so positive weight
                PointOnCurve prev = p1.Prev();
                int i_prev = loop->GlobalIndex(prev);
                int i_next = loop->GlobalIndex(p1);

                PointOnCurve next = p1.Next();

                // Compute the two weights
                Vector3 v_prev = (p1.Position() - prev.Position());
                double len_prev = norm(v_prev);
                double w_prev = 1.0 / len_prev;
                v_prev /= len_prev;
                
                // However, the "forward" edge had a negative weight, so it is also negative here
                Vector3 v_next = (next.Position() - p1.Position());
                double len_next = norm(v_next);
                double w_next = -1.0 / len_next;
                v_next /= len_next;

                double result = w_prev * dot(v_prev, SelectRow(es, i_prev)) + w_next * dot(v_next, SelectRow(es, i_next));
                out(i_next) += result;
            }
        }
    }

    void SobolevCurves::ApplyGramMatrix(PolyCurveGroup* loop, Eigen::VectorXd &as, double alpha, double beta, Eigen::VectorXd &out) {
        int n = as.size();
        Eigen::MatrixXd Df_v(n, 3);
        Df_v.setZero();

        ApplyDf(loop, as, Df_v);

        Eigen::MatrixXd Af;
        Af.setZero(n, n);

        FillAfMatrix(loop, alpha, beta, Af);

        Eigen::MatrixXd Af_Df_v(n, 3);
        Af_Df_v.setZero();
        Af_Df_v = Af * Df_v;

        Eigen::VectorXd ones;
        ones.setOnes(n);
        Eigen::VectorXd Af_1 = Af * ones;

        /*
        for (int i = 0; i < n; i++) {
            Df_v[i] = Af_1(i) * Df_v(i) - Af_Df_v(i);
        }
        */
        Df_v = Af_1.asDiagonal() * Df_v - Af_Df_v;

        ApplyDfTranspose(loop, Df_v, out);
    }

    void SobolevCurves::MultiplyComponents(Eigen::MatrixXd &A, std::vector<Vector3> &x, std::vector<Vector3> &out) {
        int n = x.size();
        Eigen::VectorXd xVec;
        xVec.setZero(n);

        // First multiply x coordinates
        for (int i = 0; i < n; i++) {
            xVec(i) = x[i].x;
        }
        Eigen::VectorXd res = A * xVec;
        for (int i = 0; i < n; i++) {
            out[i].x = res(i);
        }

        // Multiply y coordinates
        for (int i = 0; i < n; i++) {
            xVec(i) = x[i].y;
        }
        res = A * xVec;
        for (int i = 0; i < n; i++) {
            out[i].y = res(i);
        }

        // Multiply z coordinates
        for (int i = 0; i < n; i++) {
            xVec(i) = x[i].z;
        }
        res = A * xVec;
        for (int i = 0; i < n; i++) {
            out[i].z = res(i);
        }
    }
}
