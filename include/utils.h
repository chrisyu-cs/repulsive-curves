#pragma once

#include "geometrycentral/utilities/vector3.h"
#include "geometrycentral/surface/vertex_position_geometry.h"

namespace LWS {

    using namespace geometrycentral;

    Vector3 heVector(surface::VertexPositionGeometry* geom, surface::Halfedge he);

    Vector3 vector_max(Vector3 a, Vector3 b);
    Vector3 vector_min(Vector3 a, Vector3 b);
    Vector3 vector_abs(Vector3 a);

    double std_vector_dot(std::vector<double> &x, std::vector<double> &y);
    double std_vector_sum_entries(std::vector<double> &x);

    class Utils {
        public:
        static long currentTimeMilliseconds();
    };
}
