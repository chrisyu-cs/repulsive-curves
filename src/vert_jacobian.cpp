#include "vert_jacobian.h"

namespace LWS {

    Vector3 VertJacobian::LeftMultiply(Vector3 v) {
        // Interpret v as a row vector multiplied on the left.
        // Then each entry is just the dot of v with the corresponding column.
        double x = dot(v, directional_x);
        double y = dot(v, directional_y);
        double z = dot(v, directional_z);
        return Vector3{x, y, z};
    }


    void VertJacobian::Print() {
        std::cout << directional_x << "\n" << directional_y << "\n" << directional_z << std::endl;
    }
    
    double VertJacobian::Norm() {
        double n2 = norm2(directional_x) + norm2(directional_y) + norm2(directional_z);
        return sqrt(n2);
    }

    Vector3 VertJacobian::RightMultiply(Vector3 v) {
        double x = directional_x.x * v.x + directional_y.x * v.y + directional_z.x * v.z;
        double y = directional_x.y * v.x + directional_y.y * v.y + directional_z.y * v.z;
        double z = directional_x.z * v.x + directional_y.z * v.y + directional_z.z * v.z;
        return Vector3{x, y, z};
    }

    VertJacobian outer_product_to_jacobian(Vector3 v1, Vector3 v2) {
        Vector3 col1 = v1 * v2.x;
        Vector3 col2 = v1 * v2.y;
        Vector3 col3 = v1 * v2.z;
        return VertJacobian{col1, col2, col3};
    }

    VertJacobian operator+(const VertJacobian& a, const VertJacobian& b) {
        Vector3 dir_x = a.directional_x + b.directional_x;
        Vector3 dir_y = a.directional_y + b.directional_y;
        Vector3 dir_z = a.directional_z + b.directional_z;

        return VertJacobian{dir_x, dir_y, dir_z};
    }

    VertJacobian operator-(const VertJacobian& a, const VertJacobian& b) {
        return a + (b * -1);
    }

    VertJacobian operator*(const VertJacobian& a, double c) {
        return VertJacobian{a.directional_x * c, a.directional_y * c, a.directional_z * c};
    }
  
    VertJacobian operator*(double c, const VertJacobian& a) {
        return a * c;
    }
}