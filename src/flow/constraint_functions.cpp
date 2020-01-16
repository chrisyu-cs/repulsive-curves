#include "flow/constraint_functions.h"
#include "poly_curve_network.h"

#include "tpe_energy_sc.h"

namespace LWS {
    
    void ConstraintFunctions::NegativeBarycenterViolation(PolyCurveNetwork* curves,
    Eigen::VectorXd &b, Eigen::VectorXd &targets, int rowStart) {
        Vector3 barycenter = curves->Barycenter();

        // Fill in difference from barycenter to target point
        b(rowStart + 0) = targets(rowStart + 0) - barycenter.x;
        b(rowStart + 1) = targets(rowStart + 1) - barycenter.y;
        b(rowStart + 2) = targets(rowStart + 2) - barycenter.z;
    }

    void ConstraintFunctions::NegativeTotalLengthViolation(PolyCurveNetwork* curves,
    Eigen::VectorXd &b, Eigen::VectorXd &targets, int rowStart) {
        double totalLength = curves->TotalLength();
        b(rowStart) = targets(rowStart) - totalLength;
    }

    void ConstraintFunctions::NegativeEdgeLengthViolation(PolyCurveNetwork* curves,
    Eigen::VectorXd &b, Eigen::VectorXd &targets, int rowStart) {
        int nEdges = curves->NumEdges();

        // For each edge, fill in its deviation from target length
        for (int i = 0; i < nEdges; i++) {
            CurveEdge* e_i = curves->GetEdge(i);
            double curLen = e_i->Length();
            int id = e_i->id;
            double negError = targets(rowStart + id) - curLen;
            b(rowStart + id) = negError;
        }
    }

    void ConstraintFunctions::NegativePinViolation(PolyCurveNetwork* curves,
    Eigen::VectorXd &b, Eigen::VectorXd &targets, int rowStart) {
        int nPins = curves->NumPins();

        for (int i = 0; i < nPins; i++) {
            CurveVertex* v_i = curves->GetPinnedVertex(i);
            Vector3 p_i = v_i->Position();
            b(rowStart + 3 * i    ) = targets(rowStart + 3 * i    ) - p_i.x;
            b(rowStart + 3 * i + 1) = targets(rowStart + 3 * i + 1) - p_i.y;
            b(rowStart + 3 * i + 2) = targets(rowStart + 3 * i + 2) - p_i.z;
        }
    }

    void ConstraintFunctions::NegativeTangentViolation(PolyCurveNetwork* curves,
    Eigen::VectorXd &b, Eigen::VectorXd &targets, int rowStart) {
        int nTanPins = curves->NumTangentPins();
        
        for (int i = 0; i < nTanPins; i++) {
            CurveVertex* v_i = curves->GetPinnedTangent(i);
            Vector3 t_i = v_i->Tangent();
            b(rowStart + 3 * i    ) = targets(rowStart + 3 * i    ) - t_i.x;
            b(rowStart + 3 * i + 1) = targets(rowStart + 3 * i + 1) - t_i.y;
            b(rowStart + 3 * i + 2) = targets(rowStart + 3 * i + 2) - t_i.z;
        }
    }

    void ConstraintFunctions::NegativeSurfaceViolation(PolyCurveNetwork* curves, Eigen::VectorXd &b, Eigen::VectorXd &targets, int rowStart) {
        int nPins = curves->NumPinnedToSurface();

        for (int i = 0; i < nPins; i++) {
            CurveVertex* v = curves->GetPinnedToSurface(i);
            int row = rowStart + v->GlobalIndex();
            double distance = curves->constraintSurface->SignedDistance(v->Position());
            b(row) = -distance;
        }
    }

    void ConstraintFunctions::AddBarycenterTriplets3X(PolyCurveNetwork* curves,
    std::vector<Eigen::Triplet<double>> &triplets, int rowStart) {
        int nVerts = curves->NumVertices();
        double totalLength = curves->TotalLength();
        // Fill a single row with normalized vertex weights
        for (int i = 0; i < nVerts; i++) {
            double wt = curves->GetVertex(i)->DualLength() / totalLength;
            triplets.push_back(Eigen::Triplet<double>(rowStart + 0, 3 * i, wt));
            triplets.push_back(Eigen::Triplet<double>(rowStart + 1, 3 * i + 1, wt));
            triplets.push_back(Eigen::Triplet<double>(rowStart + 2, 3 * i + 2, wt));
        }
    }

    void ConstraintFunctions::AddTotalLengthTriplets(PolyCurveNetwork* curves,
    std::vector<Eigen::Triplet<double>> &triplets, int rowStart) {
        int nVerts = curves->NumVertices();
        for (int i = 0; i < nVerts; i++) {
            CurveVertex* v_i = curves->GetVertex(i);
            int nNeighbors = v_i->numEdges();
            Vector3 gradient{0, 0, 0};

            for (int j = 0; j < nNeighbors; j++) {
                CurveVertex* v_j = v_i->edge(j)->Opposite(v_i);
                Vector3 grad_j = v_i->Position() - v_j->Position();
                grad_j = grad_j.normalize();

                gradient += grad_j;
            }

            int index = v_i->GlobalIndex();
            triplets.push_back(Eigen::Triplet<double>(rowStart, 3 * index,     gradient.x));
            triplets.push_back(Eigen::Triplet<double>(rowStart, 3 * index + 1, gradient.y));
            triplets.push_back(Eigen::Triplet<double>(rowStart, 3 * index + 2, gradient.z));
        }
    }

    void ConstraintFunctions::AddEdgeLengthTriplets(PolyCurveNetwork* curves,
    std::vector<Eigen::Triplet<double>> &triplets, int rowStart) {
        int nEdges = curves->NumEdges();
        
        // Add the edge length rows
        for (int i = 0; i < nEdges; i++) {
            CurveEdge* edge = curves->GetEdge(i);
            CurveVertex* pt1 = edge->prevVert;
            CurveVertex* pt2 = edge->nextVert;

            // This is the gradient of edge length wrt pt1; the gradient wrt pt2 is just negative of this.
            Vector3 grad1 = pt1->Position() - pt2->Position();
            grad1 = grad1.normalize();

            int j1 = pt1->GlobalIndex();
            int j2 = pt2->GlobalIndex();
            int start = i + rowStart;

            // Write the three gradient entries for pt1 into the row
            triplets.push_back(Eigen::Triplet<double>(start, 3 * j1,     grad1.x));
            triplets.push_back(Eigen::Triplet<double>(start, 3 * j1 + 1, grad1.y));
            triplets.push_back(Eigen::Triplet<double>(start, 3 * j1 + 2, grad1.z));

            // Similarly write the three gradient entries for pt2 into the same row
            triplets.push_back(Eigen::Triplet<double>(start, 3 * j2,     -grad1.x));
            triplets.push_back(Eigen::Triplet<double>(start, 3 * j2 + 1, -grad1.y));
            triplets.push_back(Eigen::Triplet<double>(start, 3 * j2 + 2, -grad1.z));
        }
    }

    void ConstraintFunctions::AddPinTriplets(PolyCurveNetwork* curves,
    std::vector<Eigen::Triplet<double>> &triplets, int rowStart) {
        int nPins = curves->NumPins();

        for (int i = 0; i < nPins; i++) {
            int id = curves->GetPinnedVertex(i)->id;
            triplets.push_back(Eigen::Triplet<double>(rowStart + 3 * i,     3 * id,     1));
            triplets.push_back(Eigen::Triplet<double>(rowStart + 3 * i + 1, 3 * id + 1, 1));
            triplets.push_back(Eigen::Triplet<double>(rowStart + 3 * i + 2, 3 * id + 2, 1));
        }
    }

    void AddJacobianTriplets(std::vector<Eigen::Triplet<double>> &triplets, VertJacobian &J, int row, int col) {
        triplets.push_back(Eigen::Triplet<double>(row,     col,     J.directional_x.x));
        triplets.push_back(Eigen::Triplet<double>(row,     col + 1, J.directional_y.x));
        triplets.push_back(Eigen::Triplet<double>(row,     col + 2, J.directional_z.x));
        
        triplets.push_back(Eigen::Triplet<double>(row + 1, col,     J.directional_x.y));
        triplets.push_back(Eigen::Triplet<double>(row + 1, col + 1, J.directional_y.y));
        triplets.push_back(Eigen::Triplet<double>(row + 1, col + 2, J.directional_z.y));

        triplets.push_back(Eigen::Triplet<double>(row + 2, col,     J.directional_x.z));
        triplets.push_back(Eigen::Triplet<double>(row + 2, col + 1, J.directional_y.z));
        triplets.push_back(Eigen::Triplet<double>(row + 2, col + 2, J.directional_z.z));
    }

    void ConstraintFunctions::AddTangentTriplets(PolyCurveNetwork* curves,
    std::vector<Eigen::Triplet<double>> &triplets, int rowStart) {
        int nTanPins = curves->NumTangentPins();

        for (int i = 0; i < nTanPins; i++) {
            CurveVertex* v = curves->GetPinnedTangent(i);
            int nNeighbors = v->numEdges();
            for (int j = 0; j < nNeighbors; j++) {
                CurveVertex* v_neighbor = v->edge(j)->Opposite(v);
                int id = v_neighbor->id;
                VertJacobian derivTangent = TPESC::vertex_tangent_wrt_vert(v, v_neighbor);
                AddJacobianTriplets(triplets, derivTangent, rowStart + 3 * i, 3 * id);
            }
            VertJacobian derivTangentCenter = -1 * TPESC::vertex_tangent_wrt_vert(v, v);
            AddJacobianTriplets(triplets, derivTangentCenter, rowStart + 3 * i, 3 * v->id);

        }
    }

    void ConstraintFunctions::AddSurfaceTriplets(PolyCurveNetwork* curves, std::vector<Eigen::Triplet<double>> &triplets, int rowStart) {
        int nPins = curves->NumPinnedToSurface();

        for (int i = 0; i < nPins; i++) {
            CurveVertex* v = curves->GetPinnedToSurface(i);
            Vector3 gradient = curves->constraintSurface->GradientOfDistance(v->Position());
            int id = v->GlobalIndex();
            int row = rowStart + i;

            triplets.push_back(Eigen::Triplet<double>(row, 3 * id, gradient.x));
            triplets.push_back(Eigen::Triplet<double>(row, 3 * id + 1, gradient.y));
            triplets.push_back(Eigen::Triplet<double>(row, 3 * id + 2, gradient.z));
        }
    }

    void ConstraintFunctions::SetBarycenterTargets(PolyCurveNetwork* curves, Eigen::VectorXd &targets, int rowStart) {
        Vector3 bcenter = curves->Barycenter();
        // Set the three barycenter target entries to 0
        targets(rowStart + 0) = bcenter.x;
        targets(rowStart + 1) = bcenter.y;
        targets(rowStart + 2) = bcenter.z;
    }

    void ConstraintFunctions::SetTotalLengthTargets(PolyCurveNetwork* curves, Eigen::VectorXd &targets, int rowStart) {
        targets(rowStart) = curves->TotalLength();
    }

    void ConstraintFunctions::SetEdgeLengthTargets(PolyCurveNetwork* curves, Eigen::VectorXd &targets, int rowStart) {
        int nEdges = curves->NumEdges();
        for (int i = 0; i < nEdges; i++) {
            targets(rowStart + i) = curves->GetEdge(i)->Length();
        }
    }

    void ConstraintFunctions::SetPinTargets(PolyCurveNetwork* curves, Eigen::VectorXd &targets, int rowStart) {
        int nPins = curves->NumPins();
        for (int i = 0; i < nPins; i++) {
            Vector3 p_i = curves->GetPinnedVertex(i)->Position();
            targets(rowStart + 3 * i    ) = p_i.x;
            targets(rowStart + 3 * i + 1) = p_i.y;
            targets(rowStart + 3 * i + 2) = p_i.z;
        }
    }

    void ConstraintFunctions::SetTangentTargets(PolyCurveNetwork* curves, Eigen::VectorXd &targets, int rowStart) {
        int nTanPins = curves->NumTangentPins();
        for (int i = 0; i < nTanPins; i++) {
            Vector3 p_i = curves->GetPinnedTangent(i)->Tangent();
            targets(rowStart + 3 * i    ) = p_i.x;
            targets(rowStart + 3 * i + 1) = p_i.y;
            targets(rowStart + 3 * i + 2) = p_i.z;
        }
    }

    void ConstraintFunctions::SetSurfaceTargets(PolyCurveNetwork* curves, Eigen::VectorXd &targets, int rowStart) {
        int nPins = curves->NumPinnedToSurface();
        for (int i = 0; i < nPins; i++) {
            targets(rowStart + i) = 0;
        }
    }

}