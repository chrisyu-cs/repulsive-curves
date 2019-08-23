#pragma once

#include "geometrycentral/utilities/vector3.h"
#include "geometrycentral/surface/vertex_position_geometry.h"

namespace LWS {

    using namespace geometrycentral;

    Vector3 heVector(surface::VertexPositionGeometry* geom, surface::Halfedge he);

    Vector3 vector_max(Vector3 a, Vector3 b);
    Vector3 vector_min(Vector3 a, Vector3 b);
    Vector3 vector_abs(Vector3 a);

    // Dot product of two vectors of numbers
    double std_vector_dot(std::vector<double> &x, std::vector<double> &y);
    // Sum of all entries in the given vector
    double std_vector_sum_entries(std::vector<double> &x);
    // Adds vectors x and y (per-entry) and stores result in result
    void std_vector_add(std::vector<double> &x, std::vector<double> &y, std::vector<double> &result);

    class Utils {
        public:
        static long currentTimeMilliseconds();
    };
}
