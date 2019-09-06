#include "poly_curve.h"
#include "utils.h"

namespace LWS {

    bool PointOnCurve::operator ==(const PointOnCurve &other) {
        return (other.curve->index == curve->index) && (other.pIndex == pIndex);
    }

    bool PointOnCurve::operator !=(const PointOnCurve &other) {
        return !((other.curve->index == curve->index) && (other.pIndex == pIndex));
    }

    Vector3 PointOnCurve::Position() {
        return curve->Position(pIndex);
    }

    void PointOnCurve::SetPosition(Vector3 pos) {
        curve->SetPosition(pIndex, pos);
    }

    Vector3 PointOnCurve::Tangent() {
        return curve->VertexTangent(pIndex);
    }

    double PointOnCurve::DualLength() {
        return curve->DualLength(pIndex);
    }

    PointOnCurve PointOnCurve::Prev() {
        return PointOnCurve{curve->PrevVertex(pIndex), curve};
    }

    PointOnCurve PointOnCurve::Next() {
        return PointOnCurve{curve->NextVertex(pIndex), curve};
    }

    PolyCurve::PolyCurve(int numVerts) : positions(numVerts) {
        index = 0;
    }

    PolyCurve::PolyCurve(std::vector<Vector3> ps) : positions(ps) {
        index = 0;
    }

    int PolyCurve::NumVertices() {
        return positions.size();
    }

    int PolyCurve::NextVertex(int v) {
        if (v == NumVertices() - 1) {
            return 0;
        }
        else return v + 1;
    }

    int PolyCurve::PrevVertex(int v) {
        if (v == 0) {
            return NumVertices() - 1;
        }
        else return v - 1;
    }

    PointOnCurve PolyCurve::GetCurvePoint(int v) {
        return PointOnCurve{v, this};
    }
    
    Vector3 PolyCurve::Position(int v) {
        return positions[v];
    }

    void PolyCurve::SetPosition(int v, Vector3 pos) {
        positions[v] = pos;
    }

    Vector3 PolyCurve::BisectorNormal(int index) {
        int prev = PrevVertex(index);
        int next = NextVertex(index);

        Vector3 v1 = Position(prev);
        Vector3 v2 = Position(index);
        Vector3 v3 = Position(next);

        Vector3 e1 = (v2 - v1);
        Vector3 e2 = (v3 - v2);

        Vector3 up{0, 1, 0};

        Vector3 n1 = cross(up, e1);
        n1 = n1.normalize();
        Vector3 n2 = cross(up, e2);
        n2 = n2.normalize();

        // Return the angle bisector normal
        Vector3 normal = n1 + n2;
        normal = normal.normalize();
        return normal;
    }

    Vector3 PolyCurve::LengthWeightedNormal(int index) {
        int prev = PrevVertex(index);
        int next = NextVertex(index);

        Vector3 v1 = Position(prev);
        Vector3 v2 = Position(index);
        Vector3 v3 = Position(next);

        Vector3 e1 = (v2 - v1);
        Vector3 e2 = (v3 - v2);
        double w1 = norm(e1);
        double w2 = norm(e2);

        Vector3 up{0, 1, 0};

        Vector3 n1 = cross(up, e1);
        n1 = n1.normalize();
        Vector3 n2 = cross(up, e2);
        n2 = n2.normalize();

        // Weight normals by length
        Vector3 normal = (w1 * n1 + w2 * n2) / (w1 + w2);
        normal = normal.normalize();
        return normal;
    }

    Vector3 PolyCurve::EdgeNormal(int index) {
        int next = NextVertex(index);

        Vector3 edgeVec = Position(next) - Position(index);
        Vector3 up{0, 1, 0};
        Vector3 normal = cross(up, edgeVec);
        normal = normal.normalize();
        return normal;
    }

    Vector3 PolyCurve::VertexTangent(int index) {
        Vector3 prev = Position(PrevVertex(index));
        Vector3 next = Position(NextVertex(index));
        Vector3 center = Position(index);

        Vector3 prevT = center - prev;
        Vector3 nextT = next - center;
        prevT = prevT.normalize();
        nextT = nextT.normalize();

        Vector3 tangent = prevT + nextT;
        tangent = tangent.normalize();
        return tangent;
    }

    double PolyCurve::TotalLength() {
        double length = 0;
        for (int i = 0; i < NumVertices(); i++) {
            int i_next = NextVertex(i);
            Vector3 v_i = Position(i);
            Vector3 v_n = Position(i_next);
            length += norm(v_i - v_n);
        }
        return length;
    }

    double PolyCurve::DualLength(int i) {
        int i_next = NextVertex(i);
        int i_prev = PrevVertex(i);

        Vector3 v_i = Position(i);
        Vector3 v_p = Position(i_prev);
        Vector3 v_n = Position(i_next);

        double l1 = norm(v_i - v_p);
        double l2 = norm(v_i - v_n);
        return (l1 + l2) / 2;
    }

    Vector3 PolyCurve::Barycenter() {
        Vector3 center{0, 0, 0};
        int nVerts = NumVertices();

        for (int i = 0; i < NumVertices(); i++) {
            center += Position(i);
        }
        return center / nVerts;
    }

    Vector3 PolyCurve::TotalLengthGradient(int i) {
        // Moving vertex i only affects the lengths of the two edges
        // on either side of i.
        int prev = PrevVertex(i);
        int next = NextVertex(i);
    
        // Directions are away from the neighboring vertices.
        Vector3 prevEdge = Position(i) - Position(prev);
        Vector3 nextEdge = Position(i) - Position(next);
        // Magnitudes are both 1, since length changes at 1 unit / unit.
        prevEdge = prevEdge.normalize();
        nextEdge = nextEdge.normalize();

        return prevEdge + nextEdge;
    }

    PolyCurve* PolyCurve::Coarsen() {
        int nVerts = positions.size();
        bool isOddNumber = (nVerts % 2) == 1;

        int coarseVerts = (isOddNumber) ? nVerts / 2 + 1 : nVerts / 2;

        std::vector<Eigen::Triplet<double>> triplets;

        std::cout << "Coarsening " << nVerts << " -> " << coarseVerts << std::endl;

        for (int i = 0; i < coarseVerts; i++) {
            int oldI = 2 * i;
            if (isOddNumber && i == 0) {
                triplets.push_back(Eigen::Triplet<double>(i, oldI, 0.75));
                triplets.push_back(Eigen::Triplet<double>(i, oldI + 1, 0.25));
            }
            else if (isOddNumber && i == coarseVerts - 1) {
                triplets.push_back(Eigen::Triplet<double>(i, oldI - 1, 0.25));
                triplets.push_back(Eigen::Triplet<double>(i, oldI, 0.75));
            }
            else {
                triplets.push_back(Eigen::Triplet<double>(i, (oldI - 1 + nVerts) % nVerts, 0.25));
                triplets.push_back(Eigen::Triplet<double>(i, oldI, 0.5));
                triplets.push_back(Eigen::Triplet<double>(i, (oldI + 1) % nVerts, 0.25));
            }
        }

        std::cout << "Made " << triplets.size() << " triplets" << std::endl;

        Eigen::SparseMatrix<double> J(coarseVerts, nVerts);
        J.setFromTriplets(triplets.begin(), triplets.end());

        Eigen::MatrixXd positionMatrix(nVerts, 3);
        for (int i = 0; i < nVerts; i++) {
            SetRow(positionMatrix, i, positions[i]);
        }
        Eigen::MatrixXd coarsePosMat = J * positionMatrix;

        std::vector<Vector3> coarsePositions(coarsePosMat.rows());
        PolyCurve* p = new PolyCurve(coarsePosMat.rows());
        for (int i = 0; i < coarseVerts; i++) {
            p->positions[i] = SelectRow(coarsePosMat, i);
        }
        return p;
    }

    void PolyCurveGroup::AddCurve(PolyCurve* c) {
        c->index = curves.size();
        c->offset = NumVertices();
        std::cout << "Added curve (offset " << c->offset << ")" << std::endl;
        curves.push_back(c);
    }

    int PolyCurveGroup::NumVertices() {
        int verts = 0;
        for (PolyCurve* pc : curves) {
            verts += pc->NumVertices();
        }
        return verts;
    }

    PointOnCurve PolyCurveGroup::GetCurvePoint(int v) {
        // Traverse each curve sequentially
        for (size_t i = 0; i < curves.size(); i++) {
            // If the requested index is less than the length of the
            // current curve, then the point we want is on that curve
            if (v < curves[i]->NumVertices()) {
                return curves[i]->GetCurvePoint(v);
            }
            // Otherwise, we move onto the next curve, decrementing
            // the requested index
            else {
                v -= curves[i]->NumVertices();
            }
        }
        return curves[curves.size() - 1]->GetCurvePoint(v);
    }

    void PolyCurveGroup::BoundingCube(Vector3 &center, double &width) {
        PointOnCurve first = GetCurvePoint(0);
        Vector3 max_vec = first.Position();
        Vector3 min_vec = first.Position();

        // Traverse each curve sequentially
        for (size_t i = 0; i < curves.size(); i++) {
            size_t nVerts = curves[i]->NumVertices();
            for (size_t j = 0; j < nVerts; j++) {
                max_vec = vector_max(max_vec, curves[i]->positions[j]);
                min_vec = vector_min(min_vec, curves[i]->positions[j]);
            }
        }

        // Midpoint of cell is just the average of the two
        center = (max_vec + min_vec) / 2;
        // Cell width is the largest coordinate difference
        Vector3 diff = vector_abs(max_vec - min_vec);
        width = fmax(diff.x, fmax(diff.y, diff.z)) * 1.05;
        
        for (size_t i = 0; i < curves.size(); i++) {
            size_t nVerts = curves[i]->NumVertices();
            for (size_t j = 0; j < nVerts; j++) {
                Vector3 pos = curves[i]->positions[j];
                Vector3 cdist = vector_abs(center - pos);
                if (cdist.x > width / 2 || cdist.y > width / 2 || cdist.z > width / 2) {
                    std::cout << "Bounding box is not bounding" << std::endl;
                }
            }
        }
    }

    Vector3 PolyCurveGroup::Barycenter() {
        Vector3 center{0, 0, 0};
        int nVerts = 0;
        for (size_t i = 0; i < curves.size(); i++) {
            int n = curves[i]->NumVertices();
            center += curves[i]->Barycenter() * n;
            nVerts += n;
        }
        return (center / nVerts);
    }

    double PolyCurveGroup::TotalLength() {
        double length = 0;
        for (size_t i = 0; i < curves.size(); i++) {
            length += curves[i]->TotalLength();
        }
        return length;
    }

    int PolyCurveGroup::NextIndexInCurve(int v) {
        PointOnCurve p = GetCurvePoint(v);
        if (p.pIndex == p.curve->NumVertices() - 1) {
            return v + 1 - p.curve->NumVertices();
        }
        else {
            return v + 1;
        }
    }

    int PolyCurveGroup::GlobalIndex(PointOnCurve c) {
        return c.curve->offset + c.pIndex;
    }

    int PointOnCurve::GlobalIndex() {
        return curve->offset + pIndex;
    }
}
