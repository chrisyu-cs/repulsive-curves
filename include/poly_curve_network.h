#pragma once

#include <vector>
#include "geometrycentral/utilities/vector3.h"
#include "Eigen/Core"
#include "utils.h"
#include "multigrid/multigrid_operator.h"

namespace LWS
{
    using namespace geometrycentral;
    class PolyCurveNetwork;
    struct CurveVertex;

    struct CurveEdge {
        CurveVertex* nextVert;
        CurveVertex* prevVert;
        int id;

        inline CurveVertex* Opposite(CurveVertex* v);
        inline Vector3 Vector();
        inline Vector3 Tangent();
        inline double Length();
        inline Vector3 Normal(int index);
    };

    struct CurveVertex {
        PolyCurveNetwork* curve;
        int id;

        bool operator ==(const CurveVertex &other);
        bool operator !=(const CurveVertex &other);

        inline size_t numEdges();
        inline CurveEdge* edge(size_t i);
        inline Vector3 Position();
        inline void SetPosition(Vector3 pos);
        inline Vector3 Tangent();
        inline double DualLength();
        inline int GlobalIndex();
    };

    class PolyCurveNetwork {
        friend struct CurveVertex;

        public:
        Eigen::MatrixXd positions;
        
        PolyCurveNetwork(std::vector<Vector3> &ps, std::vector<std::array<size_t, 2>> &es);
        PolyCurveNetwork(Eigen::MatrixXd &ps, std::vector<std::array<size_t, 2>> &es);
        void InitStructs(std::vector<std::array<size_t, 2>> &es);
        void FindComponents();

        inline Vector3 Position(CurveVertex* v) {
            return SelectRow(positions, v->id);
        }

        inline void SetPosition(CurveVertex* v, Vector3 newPos) {
            SetPosition(v->id, newPos);
        }

        inline int NumVertices() {
            return nVerts;
        }

        inline int NumComponents() {
            return verticesByComponent.size();
        }

        inline CurveVertex* GetVertex(int i) {
            return vertices[i];
        }

        inline CurveEdge* GetEdge(int i) {
            return edges[i];
        }

        void BoundingCube(Vector3 &center, double &width);
        Vector3 Barycenter();
        double TotalLength();
        PolyCurveNetwork* Coarsen(MultigridOperator &op);

        NullSpaceProjector* constraintProjector;
        template<typename T>
        void AddConstraintProjector(GradientConstraints<T> &constraints) {
            if (constraintProjector) {
                delete constraintProjector;
            }
            constraintProjector = new NullSpaceProjector(constraints);
        }

        private:
        int nVerts;
        std::vector<CurveVertex*> vertices;
        std::vector<CurveEdge*> edges;
        std::vector<std::vector<CurveEdge*>> adjacency;

        std::vector<std::vector<CurveVertex*>> verticesByComponent;
        std::vector<std::vector<CurveEdge*>> edgesByComponent;

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

    inline size_t CurveVertex::numEdges() {
        return curve->adjacency[id].size();
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
        for (size_t i = 0; i < numEdges(); i++) {
            tangent += edge(i)->Tangent();
        }
        tangent = tangent.normalize();
        return tangent;
    }
    
    inline double CurveVertex::DualLength() {
        double length = 0;
        for (size_t i = 0; i < numEdges(); i++) {
            length += edge(i)->Length();
        }
        return length / numEdges();
    }

    inline int CurveVertex::GlobalIndex() {
        return id;
    }

} // namespace LWS

