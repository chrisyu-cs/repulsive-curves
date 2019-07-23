#pragma once

#include "geometrycentral/utilities/vector3.h"

namespace LWS {
    class PolyCurve;
    using namespace geometrycentral;

    struct PointOnCurve {
        int pIndex;
        PolyCurve* curve;
        bool operator ==(const PointOnCurve &other);
        bool operator !=(const PointOnCurve &other);
        Vector3 Position();
        void SetPosition(Vector3 pos);
        Vector3 Tangent();
        double DualLength();

        int GlobalIndex();

        PointOnCurve Next();
        PointOnCurve Prev();
    };

    class PolyCurve {
        public:
        std::vector<Vector3> positions;
        PolyCurve(int numVerts);
        PolyCurve(std::vector<Vector3> ps);
        int index;
        int offset;

        int NumVertices();
        int NextVertex(int v);
        int PrevVertex(int v);
        PointOnCurve GetCurvePoint(int v);

        Vector3 Position(int v);
        void SetPosition(int v, Vector3 pos);

        Vector3 BisectorNormal(int index);
        Vector3 LengthWeightedNormal(int index);
        Vector3 EdgeNormal(int index);
        Vector3 VertexTangent(int index);
        double TotalLength();
        double DualLength(int i);
        Vector3 Barycenter();
        Vector3 TotalLengthGradient(int i);
    };

    class PolyCurveGroup {
        public:
        std::vector<PolyCurve*> curves;
        void AddCurve(PolyCurve* c);
        int NumVertices();
        void BoundingCube(Vector3 &center, double &width);
        PointOnCurve GetCurvePoint(int v); 
        Vector3 Barycenter();
        double TotalLength();
        int NextIndexInCurve(int v);
        int GlobalIndex(PointOnCurve c);
    };
}