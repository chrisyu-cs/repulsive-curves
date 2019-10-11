#include "poly_curve_network.h"

#include <queue>

namespace LWS {

    bool CurveVertex::operator ==(const CurveVertex &other) {
        return (other.curve == curve) && (other.id == id);
    }

    bool CurveVertex::operator !=(const CurveVertex &other) {
        return !((other.curve == curve) && (other.id == id));
    }

    PolyCurveNetwork::PolyCurveNetwork(std::vector<Vector3> &ps, std::vector<std::array<size_t, 2>> &es) {
        nVerts = ps.size();
        positions.setZero(nVerts, 3);
        for (int i = 0; i < nVerts; i++) {
            SetRow(positions, i, ps[i]);
        }
        adjacency = std::vector<std::vector<CurveEdge*>>(nVerts);

        InitStructs(es);
        FindComponents();
    }

    PolyCurveNetwork::PolyCurveNetwork(Eigen::MatrixXd &ps, std::vector<std::array<size_t, 2>> &es) {
        positions = ps;
        adjacency = std::vector<std::vector<CurveEdge*>>(nVerts);
        
        InitStructs(es);
        FindComponents();
    }

    void PolyCurveNetwork::InitStructs(std::vector<std::array<size_t, 2>> &es) {

        // Create all vertex structs
        for (int i = 0; i < nVerts; i++) {
            vertices.push_back(new CurveVertex());
            vertices[i]->id = i;
            vertices[i]->curve = this;
        }
        // Create all edge structs
        for (size_t i = 0; i < es.size(); i++) {
            edges.push_back(new CurveEdge());
            // Assume that the edge is oriented towards the second vertex,
            // and away from the first vertex
            edges[i]->prevVert = vertices[es[i][0]];
            edges[i]->nextVert = vertices[es[i][1]];
            adjacency[es[i][0]].push_back(edges[i]);
            adjacency[es[i][1]].push_back(edges[i]);
            edges[i]->id = i;
        }
    }

    void PolyCurveNetwork::FindComponents() {
        std::vector<bool> foundVerts(nVerts);
        std::vector<bool> foundEdges(edges.size());

        for (int i = 0; i < nVerts; i++) {
            foundVerts[i] = false;
        }

        for (size_t i = 0; i < edges.size(); i++) {
            foundEdges[i] = false;
        }

        for (int i = 0; i < nVerts; i++) {
            if (foundVerts[i]) continue;
            else {
                int cNum = verticesByComponent.size();
                verticesByComponent.push_back(std::vector<CurveVertex*>());
                edgesByComponent.push_back(std::vector<CurveEdge*>());

                CurveVertex* start = vertices[i];
                std::queue<CurveVertex*> frontier;
                frontier.push(start);
                // Do the BFS
                while (!frontier.empty()) {
                    CurveVertex* next = frontier.front();
                    frontier.pop();
                    if (foundVerts[next->id]) continue;
                    // Group vertices by component
                    foundVerts[next->id] = true;
                    verticesByComponent[cNum].push_back(next);
                    // Make the start vertex an endpoint, if we find it
                    if (next->numEdges() == 1) {
                        start = next;
                    }
                    // Add neighbors
                    for (int e = 0; e < next->numEdges(); e++) {
                        CurveEdge* nEdge = next->edge(e);
                        // Group edges by component
                        if (!foundEdges[nEdge->id]) {
                            foundEdges[nEdge->id] = true;
                            edgesByComponent[cNum].push_back(nEdge);
                        }
                        // Add neighbors as well
                        CurveVertex* neighbor = nEdge->Opposite(next);
                        if (foundVerts[neighbor->id]) continue;
                        frontier.push(neighbor);
                    }
                }

                std::cout << "Found component starting at " << start->id << ", with "
                    << verticesByComponent[cNum].size() << " vertices and "
                    << edgesByComponent[cNum].size() << " edges" << std::endl;
            }
        }
    }

    void PolyCurveNetwork::BoundingCube(Vector3 &center, double &width) {
        Vector3 min_coords;
        Vector3 max_coords;

        // Find min and max coordinates
        for (int i = 0; i < nVerts; i++) {
            Vector3 coords = SelectRow(positions, i);
            min_coords = vector_min(min_coords, coords);
            max_coords = vector_max(max_coords, coords);
        }

        // Midpoint of cell is just the average of the two
        center = (min_coords + max_coords) / 2;
        // Cell width is the largest coordinate difference
        Vector3 diff = vector_abs(max_coords - min_coords);
        width = fmax(diff.x, fmax(diff.y, diff.z)) * 1.05;
    }

    Vector3 PolyCurveNetwork::Barycenter() {
        Vector3 center{0, 0, 0};
        for (int i = 0; i < nVerts; i++) {
            center += Position(i);
        }
        return center / nVerts;
    }

    double PolyCurveNetwork::TotalLength() {
        double length = 0;
        for (size_t i = 0; i < edges.size(); i++) {
            length += edges[i]->Length();
        }
        return length;
    }

    PolyCurveNetwork* PolyCurveNetwork::Coarsen(MultigridOperator &op, bool doEdgeMatrix) {
        Eigen::SparseMatrix<double> prolongMatrix, edgeMatrix;
        // Determine which vertices should be kept and which shouldn't
        std::vector<bool> shouldKeep(nVerts);
        std::vector<bool> explored(nVerts);
        std::queue<std::pair<CurveVertex*, bool>> frontier;
        // Map coarsened indices back to finer ones in this curve
        std::vector<int> fineToCoarse(nVerts);
        int coarseCount = 0;

        for (int i = 0; i < NumComponents(); i++) {
            frontier.push(std::pair<CurveVertex*, bool>(verticesByComponent[i][0], true));
        }

        while (!frontier.empty()) {
            std::pair<CurveVertex*, bool> next_pair = frontier.front();
            CurveVertex* next = next_pair.first;
            bool tentativeKeep = next_pair.second;
            frontier.pop();
            if (!explored[next->id]) {
                explored[next->id] = true;
                // Keep it either if we previous said we should, or if
                // the vertex is an "interesting" one (i.e. degree neq 2)
                bool keep = tentativeKeep || (next->numEdges() != 2);
                shouldKeep[next->id] = keep;
                if (keep) {
                    int coarseID = coarseCount;
                    coarseCount++;
                    fineToCoarse[next->id] = coarseID;
                }
                for (int e = 0; e < next->numEdges(); e++) {
                    CurveVertex* neighbor = next->edge(e)->Opposite(next);
                    if (!explored[neighbor->id]) {
                        frontier.push(std::pair<CurveVertex*, bool>(neighbor, !keep));
                    }
                }
            }
        }

        std::vector<Eigen::Triplet<double>> triplets;

        // Assemble the prolongation operator
        for (int i = 0; i < nVerts; i++) {
            if (shouldKeep[i]) {
                // If we kept the vertex, then we need to look up the index of
                // the corresponding vertex in the coarsened network
                int coarseI = fineToCoarse[i];
                triplets.push_back(Eigen::Triplet<double>(i, coarseI, 1));
            }
            else {
                // Otherwise, if we deleted the vertex, we need to find the index
                // of the two neighboring coarsened vertices
                CurveVertex* v_i = vertices[i];
                if (v_i->numEdges() != 2) {
                    std::cerr << "Deleted an edge with vertex valuence != 2" << std::endl;
                    throw 1;
                }
                int coarseI1 = v_i->edge(0)->Opposite(v_i)->id;
                int coarseI2 = v_i->edge(1)->Opposite(v_i)->id;
                triplets.push_back(Eigen::Triplet<double>(i, coarseI1, 0.5));
                triplets.push_back(Eigen::Triplet<double>(i, coarseI2, 0.5));
            }
        }

        prolongMatrix.resize(vertices.size(), coarseCount);
        prolongMatrix.setFromTriplets(triplets.begin(), triplets.end());
        op.matrices.push_back(IndexedMatrix{prolongMatrix, 0, 0});

        if (doEdgeMatrix) {

        }
        
        // Assemble the edge prolongation operator
        std::vector<Eigen::Triplet<double>> edgeTriplets;
        std::vector<std::array<size_t, 2>> coarseEdges;

        for (int i = 0; i < nVerts; i++) {
            if (shouldKeep[i]) {
                // Look for neighbors that have both endpoints kept
                CurveVertex* v_i = vertices[i];
                for (int e = 0; e < v_i->numEdges(); i++) {
                    CurveEdge* e_i = v_i->edge(e);
                    CurveVertex* v_neighbor = e_i->Opposite(v_i);
                    // If both endpoints of the edge were kept,
                    // then map the old one to the new one with
                    // weight 1
                    if (shouldKeep[v_neighbor->id]) {
                        int coarseID = coarseEdges.size();
                        size_t coarsePrev = fineToCoarse[e_i->prevVert->id];
                        size_t coarseNext = fineToCoarse[e_i->nextVert->id];
                        coarseEdges.push_back({coarsePrev, coarseNext});
                        int fineID = e_i->id;
                        if (doEdgeMatrix) {
                            edgeTriplets.push_back(Eigen::Triplet<double>(fineID, coarseID, 1));
                        }
                    }
                }
            }
            else {
                // This means we deleted a vertex, so both its neighbors were kept,
                // and one longer edge runs between them
                CurveVertex* v_i = vertices[i];
                size_t coarsePrev, coarseNext;
                if (v_i == v_i->edge(0)->nextVert) {
                    coarsePrev = v_i->edge(0)->prevVert->id;
                    coarseNext = v_i->edge(1)->nextVert->id;
                }
                else {
                    coarsePrev = v_i->edge(0)->nextVert->id;
                    coarseNext = v_i->edge(1)->prevVert->id;
                }
                int coarseID = coarseEdges.size();
                coarseEdges.push_back({coarsePrev, coarseNext});
                if (doEdgeMatrix) {
                    edgeTriplets.push_back(Eigen::Triplet<double>(v_i->edge(0)->id, coarseID, 0.5));
                    edgeTriplets.push_back(Eigen::Triplet<double>(v_i->edge(1)->id, coarseID, 0.5));
                }
            }
        }

        if (doEdgeMatrix) {
            edgeMatrix.resize(edges.size(), coarseEdges.size());
            edgeMatrix.setFromTriplets(edgeTriplets.begin(), edgeTriplets.end());
            op.edgeMatrices.push_back(IndexedMatrix{edgeMatrix, 0, 0});
        }

        // Get the coarsened positions
        Eigen::MatrixXd coarsePosMat = ApplyPinv(prolongMatrix, positions);
        PolyCurveNetwork* p = new PolyCurveNetwork(coarsePosMat, coarseEdges);
        return p;
    }
}