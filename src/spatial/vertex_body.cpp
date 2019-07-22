#include "spatial/vertex_body.h"
#include "utils.h"

namespace LWS {

    Vector3 CartesianToSpherical(Vector3 c) {
        double r = norm(c);
        double theta = acos(c.z / r);
        double phi = atan2(c.y, c.x);
        return Vector3{r, theta, phi};
    }

    Vector3 SphericalToCartesian(Vector3 s) {
        // Spherical coordinates are stored as (r, theta, phi)
        double x = s.x * sin(s.y) * cos(s.z);
        double y = s.x * sin(s.y) * sin(s.z);
        double z = s.x * cos(s.y);
        return Vector3{x, y, z};
    }

    void PosTan::Print() {
        std::cout << position << "; " << tangent << std::endl;
    }

    PosTan postan_max(PosTan v1, PosTan v2) {
        return PosTan{vector_max(v1.position, v2.position), vector_max(v1.tangent, v2.tangent)};
    }

    PosTan postan_min(PosTan v1, PosTan v2) {
        return PosTan{vector_min(v1.position, v2.position), vector_min(v1.tangent, v2.tangent)};
    }
 
    BodyType VertexBody6D::type() {
        if (vertIndex2 < 0) {
            return BodyType::Vertex;
        }
        else {
            return BodyType::Edge;
        }
    }
}
