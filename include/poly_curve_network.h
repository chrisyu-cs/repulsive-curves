#pragma once

#include <vector>
#include "geometrycentral/utilities/vector3.h"
#include "Eigen/Core"
#include "utils.h"
#include "multigrid/multigrid_operator.h"
#include <unordered_set>
#include "flow/gradient_constraint_enum.h"
#include "implicit_surface.h"

namespace LWS
{
    using namespace geometrycentral;
    class PolyCurveNetwork;
    struct CurveVertex;

    struct CurveEdge {
        CurveVertex* nextVert;
        CurveVertex* prevVert;
        int id;
        int component;

        inline CurveVertex* Opposite(CurveVertex* v);
        inline Vector3 Vector();
        inline Vector3 Tangent();
        inline double Length();
        inline Vector3 Midpoint();
        inline Vector3 Normal(int index);
        inline int GlobalIndex();
        inline bool IsNeighbors(CurveEdge* other);
    };

    struct CurveVertex {
        PolyCurveNetwork* curve;
        int id;
        int component;

        bool operator ==(const CurveVertex &other);
        bool operator !=(const CurveVertex &other);

        inline int numEdges();
        inline CurveEdge* edge(int i);
        inline CurveEdge* edge(size_t i);
        inline Vector3 Position();
        inline void SetPosition(Vector3 pos);
        inline Vector3 Tangent();
        inline double DualLength();
        inline int GlobalIndex();
        inline Vector3 TotalLengthGradient();
        inline bool IsNeighbors(CurveVertex* other);
    };

    class PolyCurveNetwork {
        friend struct CurveVertex;

        public:
        PolyCurveNetwork(std::vector<Vector3> &ps, std::vector<std::array<size_t, 2>> &es);
        PolyCurveNetwork(Eigen::MatrixXd &ps, std::vector<std::array<size_t, 2>> &es);
        ~PolyCurveNetwork();
        
        void InitStructs(std::vector<std::array<size_t, 2>> &es);
        void FindComponents();
        void PinVertex(int i);
        void PinTangent(int i);
        void PinToSurface(int i);
        void PinAllSpecialVertices(bool includeTangents);
        void PinAllEndpoints(bool includeTangents);
        void PrintPins();

        bool pinnedAllToSurface;

        inline Vector3 Position(CurveVertex* v) {
            return SelectRow(positions, v->id);
        }

        inline void SetPosition(CurveVertex* v, Vector3 newPos) {
            SetPosition(v->id, newPos);
        }

        inline int NumVertices() {
            return nVerts;
        }

        inline int NumEdges() {
            return edges.size();
        }

        inline int NumPins() {
            return pinnedVertices.size();
        }

        inline int NumTangentPins() {
            return pinnedTangents.size();
        }

        inline int NumPinnedToSurface() {
            return pinnedToSurface.size();
        }

        inline int NumComponents() {
            return verticesByComponent.size();
        }

        inline int NumVerticesInComponent(int i) {
            return verticesByComponent[i].size();
        }

        inline CurveVertex* GetVertexInComponent(int c, int i) {
            return verticesByComponent[c][i];
        }

        inline CurveVertex* GetVertex(int i) {
            return vertices[i];
        }

        inline CurveVertex* GetVertex(size_t i) {
            return vertices[i];
        }

        inline CurveVertex* GetPinnedVertex(int i) {
            return vertices[pinnedVertices[i]];
        }

        inline CurveVertex* GetPinnedTangent(int i) {
            return vertices[pinnedTangents[i]];
        }

        inline CurveVertex* GetPinnedToSurface(int i) {
            return vertices[pinnedToSurface[i]];
        }

        inline bool isPinned(int i) {
            return pinnedSet.count(i) > 0;
        }

        inline bool isTangentPinned(int i) {
            return pinnedTangentSet.count(i) > 0;
        }

        inline CurveEdge* GetEdge(int i) {
            return edges[i];
        }

        inline CurveEdge* GetEdge(size_t i) {
            return edges[i];
        }

        inline double BendingAngle(int i) {
            CurveVertex* v_i = GetVertex(i);
            if (v_i->numEdges() != 2) return 0;

            CurveVertex* v_prev = v_i->edge(0)->Opposite(v_i);
            CurveVertex* v_next = v_i->edge(1)->Opposite(v_i);

            Vector3 e_prev = (v_i->Position() - v_prev->Position()).normalize();
            Vector3 e_next = (v_next->Position() - v_i->Position()).normalize();

            double sin_angle = norm(cross(e_prev, e_next));
            double cos_angle = dot(e_prev, e_next);
            return atan2(sin_angle, cos_angle);
        }

        void BoundingCube(Vector3 &center, double &width);
        Vector3 Barycenter();
        double TotalLength();
        Vector3 AreaVector();
        PolyCurveNetwork* Subdivide();
        PolyCurveNetwork* Coarsen(MultigridOperator* op, bool doEdgeMatrix = false);

        NullSpaceProjector* constraintProjector;
        template<typename T>
        void AddConstraintProjector(GradientConstraints<T> &constraints) {
            if (constraintProjector) {
                delete constraintProjector;
            }
            constraintProjector = new NullSpaceProjector(constraints);
        }

        Eigen::MatrixXd positions;
        std::vector<ConstraintType> appliedConstraints;
        ImplicitSurface* constraintSurface;

        private:
        int nVerts;
        std::vector<int> pinnedToSurface;
        std::vector<int> pinnedVertices;
        std::unordered_set<int> pinnedSet;
        std::vector<int> pinnedTangents;
        std::unordered_set<int> pinnedTangentSet;
        std::vector<CurveVertex*> vertices;
        std::vector<CurveEdge*> edges;
        std::vector<std::vector<CurveEdge*>> adjacency;

        std::vector<std::vector<CurveVertex*>> verticesByComponent;
        std::vector<std::vector<CurveEdge*>> edgesByComponent;

        void CleanUpStructs();

        inline Vector3 Position(int i) {
            return SelectRow(positions, i);
        }

        inline void SetPosition(int i, Vector3 newPos) {
            SetRow(positions, i, newPos);
        }
    };

    inline CurveVertex* CurveEdge::Opposite(CurveVertex* v) {
        if (nextVert == v) return prevVert;
        else if (prevVert == v) return nextVert;
        else return 0;
    }

    inline Vector3 CurveEdge::Vector() {
        return (nextVert->Position() - prevVert->Position());
    }

    inline Vector3 CurveEdge::Tangent() {
        Vector3 dir = (nextVert->Position() - prevVert->Position());
        return dir.normalize();
    }

    inline double CurveEdge::Length() {
        return norm(nextVert->Position() - prevVert->Position());
    }

    inline Vector3 CurveEdge::Midpoint() {
        return (nextVert->Position() + prevVert->Position()) / 2;
    }

    inline int CurveEdge::GlobalIndex() {
        return id;
    }

    inline bool CurveEdge::IsNeighbors(CurveEdge* other) {
        bool shared = (nextVert == other->nextVert) ||
            (nextVert == other->prevVert) ||
            (prevVert == other->nextVert) ||
            (prevVert == other->prevVert);
        return shared && (other != this);
    }

    inline int CurveVertex::numEdges() {
        return curve->adjacency[id].size();
    }

    inline CurveEdge* CurveVertex::edge(int i) {
        return curve->adjacency[id][i];
    }

    inline CurveEdge* CurveVertex::edge(size_t i) {
        return curve->adjacency[id][i];
    }

    inline Vector3 CurveVertex::Position() {
        return SelectRow(curve->positions, id);
    }

    inline void CurveVertex::SetPosition(Vector3 pos) {
        curve->SetPosition(id, pos);
    }

    inline Vector3 CurveVertex::Tangent() {
        Vector3 tangent{0, 0, 0};
        for (int i = 0; i < numEdges(); i++) {
            tangent += edge(i)->Tangent();
        }
        tangent = tangent.normalize();
        return tangent;
    }
    
    inline double CurveVertex::DualLength() {
        double length = 0;
        for (int i = 0; i < numEdges(); i++) {
            length += edge(i)->Length();
        }
        return length / numEdges();
    }

    inline int CurveVertex::GlobalIndex() {
        return id;
    }

    inline Vector3 CurveVertex::TotalLengthGradient() {
        Vector3 total{0, 0, 0};
        // Moving vertex i only affects the lengths of the two edges
        // on either side of i.
        for (int e = 0; e < numEdges(); e++) {
            CurveEdge* e_i = edge(e);
            // Directions are away from the neighboring vertices.
            Vector3 inward = Position() - e_i->Opposite(this)->Position();
            // Magnitudes are 1, since length changes at 1 unit / unit.
            inward.normalize();
            total += inward;
        }
        return total;
    }

    inline bool CurveVertex::IsNeighbors(CurveVertex* other) {
        for (int i = 0; i < numEdges(); i++) {
            CurveEdge* e_i = edge(i);
            if (e_i->Opposite(this) == other) return true;
        }
        return false;
    }

} // namespace LWS

