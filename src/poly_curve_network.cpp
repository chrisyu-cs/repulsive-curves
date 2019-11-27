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
        nVerts = ps.rows();
        positions = ps;
        adjacency = std::vector<std::vector<CurveEdge*>>(nVerts);
        
        InitStructs(es);
        FindComponents();
    }

    PolyCurveNetwork::~PolyCurveNetwork() {
        CleanUpStructs();
    }

    void PolyCurveNetwork::InitStructs(std::vector<std::array<size_t, 2>> &es) {
        constraintProjector = 0;

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

    void PolyCurveNetwork::CleanUpStructs() {
        // Delete all vertex structs
        for (int i = 0; i < nVerts; i++) {
            delete vertices[i];
        }
        // Delete all edge structs
        int nEdges = NumEdges();
        for (int i = 0; i < nEdges; i++) {
            delete edges[i];
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
                    next->component = cNum;
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
                            nEdge->component = cNum;
                        }
                        // Add neighbors as well
                        CurveVertex* neighbor = nEdge->Opposite(next);
                        if (foundVerts[neighbor->id]) continue;
                        frontier.push(neighbor);
                    }
                }
            }
        }
    }

    void PolyCurveNetwork::PinAllSpecialVertices(bool includeTangents) {
        for (CurveVertex* v : vertices) {
            if (v->numEdges() != 2) {
                PinVertex(v->id);
                if (includeTangents) {
                    PinTangent(v->id);
                }
            }
        }
    }

    void PolyCurveNetwork::PinVertex(int i) {
        // Add if not already pinned
        if (!pinnedSet.count(i)) {
            pinnedVertices.push_back(i);
            pinnedSet.insert(i);
        }
    }

    void PolyCurveNetwork::PinTangent(int i) {
        if (!pinnedTangentSet.count(i)) {
            pinnedTangents.push_back(i);
            pinnedTangentSet.insert(i);
        }
    }

    void PolyCurveNetwork::PrintPins() {
        for (size_t i = 0; i < pinnedVertices.size(); i++) {
            CurveVertex* v = GetPinnedVertex(i);
            std::cout << "Pin " << i << " = " << v->id << " (position " << v->Position() << ")" << std::endl;
        }
        for (size_t i = 0; i < pinnedTangents.size(); i++) {
            CurveVertex* v = GetPinnedTangent(i);
            std::cout << "Pin " << i << " = " << v->id << " (tangent " << v->Tangent() << ")" << std::endl;
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
        double totalMass = 0;
        for (int i = 0; i < nVerts; i++) {
            double m = vertices[i]->DualLength();
            center += vertices[i]->Position() * m;
            totalMass += m;
        }
        return center / totalMass;
    }

    double PolyCurveNetwork::TotalLength() {
        double length = 0;
        for (size_t i = 0; i < edges.size(); i++) {
            length += edges[i]->Length();
        }
        return length;
    }

    Vector3 PolyCurveNetwork::AreaVector() {
        int nVerts = NumVertices();
        Vector3 sum{0, 0, 0};

        for (int i = 0; i < nVerts; i++) {
            CurveVertex* v_i = GetVertex(i);
            CurveVertex* v_next = v_i->edge(0)->Opposite(v_i);
            sum += cross(v_i->Position(), v_next->Position());
        }

        return sum / 6.0;
    }

    PolyCurveNetwork* PolyCurveNetwork::Coarsen(MultigridOperator &op, bool doEdgeMatrix) {
        Eigen::SparseMatrix<double> prolongMatrix, edgeMatrix;
        int nEdges = NumEdges();
        // Determine which vertices should be kept and which shouldn't
        std::vector<bool> shouldKeep(nVerts);
        std::vector<bool> explored(nVerts);
        std::vector<bool> seenEdges(nEdges);
        std::queue<std::pair<CurveVertex*, bool>> frontier;
        // Map fine indices to coarse ones in the next curve
        std::vector<int> fineToCoarse(nVerts);
        int coarseCount = 0;

        for (int i = 0; i < nVerts; i++) {
            shouldKeep[i] = false;
            explored[i] = false;
            fineToCoarse[i] = -1;
        }
        for (int i = 0; i < nEdges; i++) {
            seenEdges[i] = false;
        }

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
                // Keep it either if we previous said we should, if the vertex is pinned, or if
                // the vertex is an "interesting" one (i.e. degree neq 2)
                bool keep = tentativeKeep || (next->numEdges() != 2) || isPinned(next->id) || isTangentPinned(next->id);
                // Also check if any of the neighbors are already marked to be deleted
                int nEdges_i = next->numEdges();
                for (int i = 0; i < nEdges_i; i++) {
                    CurveVertex* neighbor = next->edge(i)->Opposite(next);
                    // If any neighbor is already marked to be deleted, then we need to
                    // keep this one
                    if (explored[neighbor->id] && !shouldKeep[neighbor->id]) {
                        keep = true;
                    }
                }

                shouldKeep[next->id] = keep;
                if (keep) {
                    int coarseID = coarseCount;
                    coarseCount++;
                    fineToCoarse[next->id] = coarseID;
                }
                else {
                    fineToCoarse[next->id] = -42;
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
        for (CurveVertex* v_i : vertices) {
            if (shouldKeep[v_i->id]) {
                // If we kept the vertex, then we need to look up the index of
                // the corresponding vertex in the coarsened network
                int coarseI = fineToCoarse[v_i->id];
                triplets.push_back(Eigen::Triplet<double>(v_i->id, coarseI, 1));
            }
            else {
                // Otherwise, if we deleted the vertex, we need to find the index
                // of the two neighboring coarsened vertices
                if (v_i->numEdges() != 2) {
                    std::cerr << "Deleted an edge with vertex valence != 2" << std::endl;
                    throw 1;
                }
                int n0 = v_i->edge(0)->Opposite(v_i)->id;
                int n1 = v_i->edge(1)->Opposite(v_i)->id;

                int coarseI1 = fineToCoarse[n0];
                int coarseI2 = fineToCoarse[n1];

                triplets.push_back(Eigen::Triplet<double>(v_i->id, coarseI1, 0.5));
                triplets.push_back(Eigen::Triplet<double>(v_i->id, coarseI2, 0.5));
            }
        }

        prolongMatrix.resize(vertices.size(), coarseCount);
        prolongMatrix.setFromTriplets(triplets.begin(), triplets.end());
        op.matrices.push_back(IndexedMatrix{prolongMatrix, 0, 0});
        
        // Assemble the edge prolongation operator
        std::vector<Eigen::Triplet<double>> edgeTriplets;
        std::vector<std::array<size_t, 2>> coarseEdges;

        for (CurveVertex* v_i : vertices) {
            if (shouldKeep[v_i->id]) {
                // Look for neighbors that have both endpoints kept
                for (int e = 0; e < v_i->numEdges(); e++) {
                    CurveEdge* e_i = v_i->edge(e);
                    if (seenEdges[e_i->id]) continue;
                    seenEdges[e_i->id] = true;
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
                size_t coarsePrev, coarseNext;
                int finePrev, fineNext;

                if (v_i == v_i->edge(0)->nextVert && v_i == v_i->edge(1)->prevVert) {
                    coarsePrev = fineToCoarse[v_i->edge(0)->prevVert->id];
                    coarseNext = fineToCoarse[v_i->edge(1)->nextVert->id];
                }
                else if (v_i == v_i->edge(0)->prevVert && v_i == v_i->edge(1)->nextVert) {
                    coarsePrev = fineToCoarse[v_i->edge(1)->prevVert->id];
                    coarseNext = fineToCoarse[v_i->edge(0)->nextVert->id];
                }
                else {
                    std::cerr << "Couldn't orient edge pair" << std::endl;
                    throw 1;
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

        std::vector<int> vertCounts(coarseCount);

        for (auto &e : coarseEdges) {
            vertCounts[e[0]]++;
            vertCounts[e[1]]++;
        }

        // Get the coarsened positions
        Eigen::MatrixXd coarsePosMat = ApplyPinv(prolongMatrix, positions);
        PolyCurveNetwork* p = new PolyCurveNetwork(coarsePosMat, coarseEdges);
        op.lowerSize = p->NumVertices();
        op.upperSize = NumVertices();

        // Pin corresponding vertices in the coarsened curve
        for (int pin : pinnedVertices) {
            p->PinVertex(fineToCoarse[pin]);
        }
        for (int pin : pinnedTangents) {
            p->PinTangent(fineToCoarse[pin]);
        }

        for (ConstraintType type : appliedConstraints) {
            p->appliedConstraints.push_back(type);
        }

        return p;
    }
}