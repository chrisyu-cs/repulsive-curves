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
                KfMatrix(e1, e2, points[x], points[y], alpha, beta, temp_out);
                AddToFirst(out, temp_out, weight);
            }
        } 
    }

    void SobolevCurves::FillGlobalMatrix(PolyCurve* loop, double alpha, double beta, Eigen::MatrixXd &A) {
        double out[4][4];
        int nVerts = loop->NumVertices();

        for (int i = 0; i < nVerts - 1; i++) {
            int i_next = loop->NextVertex(i);
            Vector3 p_i = loop->Position(i);
            Vector3 p_i_next = loop->Position(i_next);
            EdgePositionPair e1{p_i, p_i_next};

            for (int j = i_next; j < nVerts; j++) {
                int j_next = loop->NextVertex(j);
                Vector3 p_j = loop->Position(j);
                Vector3 p_j_next = loop->Position(j_next);
                EdgePositionPair e2{p_j, p_j_next};

                IntegrateLocalMatrix(e1, e2, alpha, beta, out);
                
                int indices[] = {i, i_next, j, j_next};
                for (int r = 0; r < 4; r++) {
                    for (int c = 0; c < 4; c++) {
                        A(indices[r], indices[c]) += out[r][c];
                    }
                }
            }
        }
    }

    void SobolevCurves::FillGlobalMatrix(PolyCurveGroup* curves, double alpha, double beta, Eigen::MatrixXd &A) {
        double out[4][4];
        int nVerts = curves->NumVertices();

        int numTerms = 0;

        for (int i = 0; i < nVerts; i++) {
            PointOnCurve pc_i = curves->GetCurvePoint(i);
            
            Vector3 p_i = pc_i.Position();
            Vector3 p_i_next = pc_i.Next().Position();
            EdgePositionPair e1{p_i, p_i_next};

            for (int j = i + 1; j < nVerts; j++) {
                PointOnCurve pc_j = curves->GetCurvePoint(j);
                if (i == j) continue;

                Vector3 p_j = pc_j.Position();
                Vector3 p_j_next = pc_j.Next().Position();
                EdgePositionPair e2{p_j, p_j_next};

                numTerms++;
                IntegrateLocalMatrix(e1, e2, alpha, beta, out);

                int indices[] = {i, curves->NextIndexInCurve(i), j, curves->NextIndexInCurve(j)};
                for (int r = 0; r < 4; r++) {
                    for (int c = 0; c < 4; c++) {
                        A(indices[r], indices[c]) += out[r][c];
                    }
                }
            }
        }
    }

    void SobolevCurves::ApproximateGlobalMatrix(PolyCurveGroup* loop, double alpha, double beta, Eigen::MatrixXd &A) {
        std::cout << "  (TODO: Approximate global matrix)" << std::endl;
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

    void SobolevCurves::ApplyDf(PolyCurveGroup* loop, std::vector<double> &as, std::vector<Vector3> &out) {
        for (PolyCurve *c : loop->curves) {
            int nVerts = c->NumVertices();
            for (int i = 0; i < nVerts; i++) {
                PointOnCurve p1 = loop->GetCurvePoint(i);
                int i1 = loop->GlobalIndex(p1);
                PointOnCurve p2 = p1.Next();
                int i2 = loop->GlobalIndex(p2);

                // Positive weight on the second vertex, negative on the first
                double d12 = as[i2] - as[i1];
                Vector3 e12 = (p2.Position() - p1.Position());
                double length = norm(e12);

                Vector3 grad12 = (d12 * e12) / (length * length);

                out[i1] = grad12;
            }
        }
    } 

    void SobolevCurves::ApplyDfTranspose(PolyCurveGroup* loop, std::vector<Vector3> &es, std::vector<double> &out) {
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
                double w_prev = len_prev;
                v_prev /= len_prev;
                // However, the "forward" edge had a negative weight, so it is also negative here
                Vector3 v_next = (next.Position() - p1.Position());
                double len_next = norm(v_next);
                double w_next = -len_next;
                v_next /= len_next;

                double result = w_prev * dot(v_prev, es[i_prev]) + w_next * dot(v_next, es[i_next]);
                out[i_next] += result;
            }
        }
    }
}
