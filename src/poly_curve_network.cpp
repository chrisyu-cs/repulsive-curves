#include "poly_curve_network.h"

namespace LWS {

    bool CurvePoint::operator ==(const CurvePoint &other) {
        return (other.curve == curve) && (other.pIndex == pIndex);
    }

    bool CurvePoint::operator !=(const CurvePoint &other) {
        return !((other.curve == curve) && (other.pIndex == pIndex));
    }

    Vector3 CurvePoint::Position() {
        return curve->Position(pIndex);
    }

    void CurvePoint::SetPosition(Vector3 pos) {
        curve->SetPosition(pIndex, pos);
    }

    Vector3 CurvePoint::Tangent() {
        return curve->VertexTangent(pIndex);
    }

    double CurvePoint::DualLength() {
        return curve->DualLength(pIndex);
    }

    CurvePoint CurvePoint::Prev() {
        return CurvePoint{curve->PrevVertex(pIndex), curve};
    }

    CurvePoint CurvePoint::Next() {
        return CurvePoint{curve->NextVertex(pIndex), curve};
    }


    PolyCurveNetwork::PolyCurveNetwork(std::vector<Vector3> &ps, std::vector<std::array<size_t, 2>> &es) {
        nVerts = ps.size();
        positions.setZero(nVerts, 3);
        for (int i = 0; i < nVerts; i++) {
            SetRow(positions, i, ps[i]);
        }
        
        // Create all vertex structs
        for (int i = 0; i < nVerts; i++) {
            vertices.push_back(new CurveVertex());
            vertices[i]->id = i;
        }
        // Create all edge structs
        for (size_t i = 0; i < es.size(); i++) {
            edges.push_back(new CurveEdge());
            // Assume that the edge is oriented towards the second vertex,
            // and away from the first vertex
            edges[i]->prevVert = vertices[es[i][0]];
            edges[i]->nextVert = vertices[es[i][1]];
            vertices[es[i][0]]->nextEdge = edges[i];
            vertices[es[i][1]]->prevEdge = edges[i];
            edges[i]->id = i;
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

    PointOnCurve PolyCurveNetwork::GetCurvePoint(int v) {

    }
}