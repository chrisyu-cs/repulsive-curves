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

    Vector3 HatGradientOnEdge(CurveEdge* edge, CurveVertex* vertex) {
        CurveVertex* edgeStart = edge->prevVert;
        CurveVertex* edgeEnd = edge->nextVert;

        if (vertex == edgeStart) {
            Vector3 towardsVertex = edgeStart->Position() - edgeEnd->Position();
            double length = norm(towardsVertex);
            // 1 over length times normalized edge vector towards the vertex
            return towardsVertex / (length * length);
        }
        else if (vertex == edgeEnd) {
            Vector3 towardsVertex = edgeEnd->Position() - edgeStart->Position();
            double length = norm(towardsVertex);
            return towardsVertex / (length * length);
        }
        else {
            return Vector3{0, 0, 0};
        }
    }

    void SobolevCurves::AddEdgePairContribution(PolyCurveNetwork* loop, double alpha, double beta,
    CurveEdge* s, CurveEdge* t, Eigen::MatrixXd &A) {
        CurveVertex* endpoints[4] = {s->prevVert, s->nextVert, t->prevVert, t->nextVert};
        double len1 = s->Length();
        double len2 = t->Length();
        Vector3 mid1 = s->Midpoint();
        Vector3 mid2 = t->Midpoint();
        double denom = pow(norm(mid1 - mid2), beta - alpha);

        for (CurveVertex* u : endpoints) {
            for (CurveVertex* v : endpoints) {
                Vector3 u_hat_s = HatGradientOnEdge(s, u);
                Vector3 u_hat_t = HatGradientOnEdge(t, u);
                Vector3 v_hat_s = HatGradientOnEdge(s, v);
                Vector3 v_hat_t = HatGradientOnEdge(t, v);

                double numer = dot(u_hat_s - u_hat_t, v_hat_s - v_hat_t);
                int index_u = u->GlobalIndex();
                int index_v = v->GlobalIndex();

                A(index_u, index_v) += (numer / denom) * len1 * len2;
            }
        }
    }

    double HatMidpoint(CurveEdge* edge, CurveVertex* vertex) {
        CurveVertex* edgeStart = edge->prevVert;
        CurveVertex* edgeEnd = edge->nextVert;

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

    void SobolevCurves::AddEdgePairContributionLow(PolyCurveNetwork* loop, double alpha, double beta,
    CurveEdge* s, CurveEdge* t, Eigen::MatrixXd &A) {

        CurveVertex* endpoints[4] = {s->prevVert, s->nextVert, t->prevVert, t->nextVert};

        double len1 = s->Length();
        double len2 = t->Length();
        Vector3 mid_s = s->Midpoint();
        Vector3 mid_t = t->Midpoint();
        double denom = pow(norm(mid_s - mid_t), 2);

        Vector3 tangent_s = s->Tangent();

        double kf_st = TPESC::tpe_Kf_pts(mid_s, mid_t, tangent_s, alpha, beta);

        for (CurveVertex* u : endpoints) {
            for (CurveVertex* v : endpoints) {
                double u_s = HatMidpoint(s, u);
                double u_t = HatMidpoint(t, u);
                double v_s = HatMidpoint(s, v);
                double v_t = HatMidpoint(t, v);

                double numer = (u_s - u_t) * (v_s - v_t);
                int index_u = u->GlobalIndex();
                int index_v = v->GlobalIndex();

                A(index_u, index_v) += (numer / denom) * kf_st * len1 * len2;
            }
        }
    }

    void SobolevCurves::SobolevGramMatrix(PolyCurveNetwork* curves, double alpha, double beta, Eigen::MatrixXd &A, double diagEps) {
        double out[4][4];
        int nEdges = curves->NumEdges();

        int numTerms = 0;

        for (int i = 0; i < nEdges; i++) {
            CurveEdge* pc_i = curves->GetEdge(i);
            
            for (int j = 0; j < nEdges; j++) {
                CurveEdge* pc_j = curves->GetEdge(j);
                // if (pc_i == pc_j) continue;
                if (pc_i == pc_j || pc_i->IsNeighbors(pc_j)) continue;

                AddEdgePairContribution(curves, alpha, beta, pc_i, pc_j, A);
                // AddEdgePairContributionLow(curves, alpha, beta, pc_i, pc_j, A);
            }
            A(i, i) += diagEps * curves->GetVertex(i)->DualLength();
        }
    }

    void SobolevCurves::SobolevGramMatrix3X(PolyCurveNetwork* curves, double alpha, double beta, Eigen::MatrixXd &A, double diagEps) {
        Eigen::MatrixXd topLeft;
        int nVerts = curves->NumVertices();
        topLeft.setZero(nVerts, nVerts);

        SobolevGramMatrix(curves, alpha, beta, topLeft, diagEps);

        for (int i = 0; i < nVerts; i++) {
            for (int j = 0; j < nVerts; j++) {
                A(3 * i, 3 * j) = topLeft(i, j);
                A(3 * i + 1, 3 * j + 1) = topLeft(i, j);
                A(3 * i + 2, 3 * j + 2) = topLeft(i, j);
            }
        }
    }

    void SobolevCurves::SobolevPlusBarycenter(PolyCurveNetwork* loop, double alpha, double beta, Eigen::MatrixXd &A, double diagEps) {
        int nVerts = loop->NumVertices();
        // Fill the top-left block with the gram matrix
        SobolevCurves::SobolevGramMatrix(loop, alpha, beta, A, diagEps);

        double sumLength = loop->TotalLength();

        // Fill the bottom row with weights for the constraint
        for (int i = 0; i < nVerts; i++) {
            double areaWeight = loop->GetVertex(i)->DualLength() / sumLength;
            // Fill in bottom row and rightmost column
            A(i, nVerts) = areaWeight;
            A(nVerts, i) = areaWeight;
        }
    }

    void SobolevCurves::SobolevPlusBarycenter3X(PolyCurveNetwork* loop, double alpha, double beta, Eigen::MatrixXd &A, double diagEps) {
        int nVerts = loop->NumVertices();
        Eigen::MatrixXd topLeft;
        topLeft.setZero(nVerts + 1, nVerts + 1);
        SobolevPlusBarycenter(loop, alpha, beta, topLeft, diagEps);

        for (int i = 0; i < nVerts + 1; i++) {
            for (int j = 0; j < nVerts + 1; j++) {
                A(3 * i, 3 * j) = topLeft(i, j);
                A(3 * i + 1, 3 * j + 1) = topLeft(i, j);
                A(3 * i + 2, 3 * j + 2) = topLeft(i, j);
            }
        }
    }

    void SobolevCurves::AddEdgeLengthConstraints(PolyCurveNetwork* curves, Eigen::MatrixXd &A, int baseIndex) {
        int nEdges = curves->NumEdges();
        for (int i = 0; i < nEdges; i++) {
            CurveEdge* e_i = curves->GetEdge(i);
            CurveVertex* pt1 = e_i->prevVert;
            CurveVertex* pt2 = e_i->nextVert;
            // This is the gradient of edge length wrt pt1; the gradient wrt pt2 is just negative of this.
            Vector3 grad1 = pt1->Position() - pt2->Position();
            grad1 = grad1.normalize();

            int j1 = pt1->GlobalIndex();
            int j2 = pt2->GlobalIndex();
            int curRow = baseIndex + i;

            // Write the three gradient entries for pt1 into the row and column
            A(curRow, 3 * j1    ) = grad1.x;
            A(curRow, 3 * j1 + 1) = grad1.y;
            A(curRow, 3 * j1 + 2) = grad1.z;
            A(3 * j1,     curRow) = grad1.x;
            A(3 * j1 + 1, curRow) = grad1.y;
            A(3 * j1 + 2, curRow) = grad1.z;

            // Similarly write the three gradient entries for pt2 into the same row and column
            A(curRow, 3 * j2    ) = -grad1.x;
            A(curRow, 3 * j2 + 1) = -grad1.y;
            A(curRow, 3 * j2 + 2) = -grad1.z;
            A(3 * j2,     curRow) = -grad1.x;
            A(3 * j2 + 1, curRow) = -grad1.y;
            A(3 * j2 + 2, curRow) = -grad1.z;
        }
    }

    double SobolevCurves::SobolevDot(Eigen::MatrixXd &as, Eigen::MatrixXd &bs, Eigen::MatrixXd &J) {
        double dot_x = as.col(0).transpose() * J * bs.col(0);
        double dot_y = as.col(1).transpose() * J * bs.col(1);
        double dot_z = as.col(2).transpose() * J * bs.col(2);

        return dot_x + dot_y + dot_z;
    }

    void SobolevCurves::SobolevNormalize(Eigen::MatrixXd &as, Eigen::MatrixXd &J) {
        double sobo_norm = sqrt(SobolevDot(as, as, J));
        for (int i = 0; i < as.rows(); i++) {
            as(i, 0) /= sobo_norm;
            as(i, 1) /= sobo_norm;
            as(i, 2) /= sobo_norm;
        }
    }

    void SobolevCurves::SobolevOrthoProjection(Eigen::MatrixXd &as, Eigen::MatrixXd &bs, Eigen::MatrixXd &J) {
        double d1 = SobolevDot(as, bs, J);
        double d2 = SobolevDot(bs, bs, J);
        // Take the projection, and divide out the norm of b
        for (int i = 0; i < as.rows(); i++) {
            AddToRow(as, i, -(d1 / d2) * SelectRow(bs, i));
        }
    }

    void SobolevCurves::DfMatrix(PolyCurveNetwork* loop, Eigen::MatrixXd &out) {
        int nVerts = loop->NumVertices();
        int nEdges = loop->NumEdges();
        // Maps from 1 scalar per vertex to a 3-vector per edge
        out.setZero(nEdges * 3, nVerts);
        
        for (int i = 0; i < nEdges; i++) {
            CurveEdge* e_i = loop->GetEdge(i);
            CurveVertex* p1 = e_i->prevVert;
            int i1 = p1->GlobalIndex();
            CurveVertex* p2 = e_i->nextVert;
            int i2 = p2->GlobalIndex();

            // Positive weight on the second vertex, negative on the first
            Vector3 e12 = (p2->Position() - p1->Position());
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
