#pragma once

#include "geometrycentral/utilities/vector3.h" 

namespace LWS {

    using namespace geometrycentral;

    enum class BodyType {
        Vertex, Edge
    };

    struct VertexBody {
        Vector3 position;
        Vector3 tangent;
        double mass;
        int vertIndex;
    };

    struct PosTan {
        Vector3 position;
        Vector3 tangent;
        void Print();
    };

    struct VertexBody6D {
        PosTan pt;
        double mass;
        int vertIndex1;
        int vertIndex2;
        
        BodyType type();
    };

    PosTan postan_max(PosTan v1, PosTan v2);
    PosTan postan_min(PosTan v1, PosTan v2);

    Vector3 CartesianToSpherical(Vector3 cartesian);

    Vector3 SphericalToCartesian(Vector3 spherical);
}