#include "flow/constraint_functions.h"
#include "poly_curve_network.h"

namespace LWS {
    
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

}