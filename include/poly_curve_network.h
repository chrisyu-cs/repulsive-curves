#pragma once

#include <vector>
#include "geometrycentral/utilities/vector3.h"
#include "Eigen/Core"
#include "utils.h"

namespace LWS
{
    using namespace geometrycentral;
    class PolyCurveNetwork;
    struct CurveVertex;
    struct CurveEdge;


    struct CurvePoint {
        int pIndex;
        PolyCurveNetwork* curve;
        bool operator ==(const CurvePoint &other);
        bool operator !=(const CurvePoint &other);
        Vector3 Position();
        void SetPosition(Vector3 pos);
        Vector3 Tangent();
        double DualLength();

        int GlobalIndex();

        CurvePoint Next();
        CurvePoint Prev();
    };

    struct CurveVertex {
        CurveEdge* nextEdge;
        CurveEdge* prevEdge;
        int id;
    };

    struct CurveEdge {
        CurveVertex* nextVert;
        CurveVertex* prevVert;
        int id;
    };

    class PolyCurveNetwork {
        public:
        Eigen::MatrixXd positions;
        
        PolyCurveNetwork(std::vector<Vector3> &ps, std::vector<std::array<size_t, 2>> &es);

        inline Vector3 Position(int i) {
            return SelectRow(positions, i);
        }

        inline void SetPosition(int i, Vector3 newPos) {
            SetRow(positions, i, newPos);
        }

        inline int NumVertices() {
            return nVerts;
        }

        inline CurvePoint GetCurvePoint(int i) {
            return CurvePoint{i, this};
        }

        inline CurveVertex* GetVertex(int i) {
            return vertices[i];
        }

        inline CurveEdge* GetEdge(int i) {
            return edges[i];
        }

        void BoundingCube(Vector3 &center, double &width);
        Vector3 VertexTangent(int i);
        double DualLength(int i);
        int PrevVertex(int pt);
        int NextVertex(int pt);

        private:
        int nVerts;
        std::vector<CurveVertex*> vertices;
        std::vector<CurveEdge*> edges;
    };

} // namespace LWS

