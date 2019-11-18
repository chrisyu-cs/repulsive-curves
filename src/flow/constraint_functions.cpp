#include "flow/constraint_functions.h"
#include "poly_curve_network.h"

namespace LWS {
    
    void ConstraintFunctions::NegativeBarycenterViolation(PolyCurveNetwork* curves,
    Eigen::VectorXd &b, Eigen::VectorXd &targets, int rowStart) {
        Vector3 barycenter = curves->Barycenter();

        // Fill in difference from barycenter to target point
        b(rowStart + 0) = targets(rowStart + 0) - barycenter.x;
        b(rowStart + 1) = targets(rowStart + 1) - barycenter.y;
        b(rowStart + 2) = targets(rowStart + 2) - barycenter.z;
    }

    void ConstraintFunctions::NegativeEdgeLengthViolation(PolyCurveNetwork* curves,
    Eigen::VectorXd &b, Eigen::VectorXd &targets, int rowStart) {
        int nEdges = curves->NumEdges();

        // For each edge, fill in its deviation from target length
        for (int i = 0; i < nEdges; i++) {
            CurveEdge* e_i = curves->GetEdge(i);
            double curLen = e_i->Length();
            int id = e_i->id;
            double negError = targets[rowStart + id] - curLen;
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

    void ConstraintFunctions::SetBarycenterTargets(PolyCurveNetwork* curves, Eigen::VectorXd &targets, int rowStart) {
        // Set the three barycenter target entries to 0
        targets(rowStart + 0) = 0;
        targets(rowStart + 1) = 0;
        targets(rowStart + 2) = 0;
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

}