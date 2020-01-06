#include "implicit_surface.h"

#include "geometrycentral/utilities/vector2.h"

namespace LWS {

    ImplicitSurface::~ImplicitSurface() {}

    ImplicitSphere::ImplicitSphere(double r, Vector3 c) {
        radius = r;
        center = c;
    }

    double ImplicitSphere::SignedDistance(Vector3 point) {
        return (point - center).norm() - radius;
    }

    Vector3 ImplicitSphere::GradientOfDistance(Vector3 point) {
        return (point - center).normalize();
    }

    double ImplicitSphere::BoundingDiameter() {
        return radius * 2;
    }

    Vector3 ImplicitSphere::BoundingCenter() {
        return center;
    }

    ImplicitTorus::ImplicitTorus(double major, double minor, Vector3 c) {
        majorRadius = major;
        minorRadius = minor;
        center = c;
    }

    double ImplicitTorus::SignedDistance(Vector3 point) {
        // Vector3 disp = point - center;
        // Vector2 planar{disp.x, disp.z};
        // Vector2 q{norm(planar) - majorRadius, disp.y};
        // return norm(q) - minorRadius;

        Vector3 p = point - center;
        double normAll = norm2(p);
        double normXZ = norm(Vector2{p.x, p.z});
        return sqrt(normAll - 2 * majorRadius * normXZ + majorRadius * majorRadius) - minorRadius;
    }

    Vector3 ImplicitTorus::GradientOfDistance(Vector3 point) {
        Vector3 p = point - center;
        double normAll = norm2(p);
        double normXZ = norm(Vector2{p.x, p.z});
        double sqTerm = sqrt(normAll - 2 * majorRadius * normXZ + majorRadius * majorRadius);

        double partialX = 2 * p.x - (2 * p.x * majorRadius) / normXZ;
        double partialY = 2 * p.y;
        double partialZ = 2 * p.z - (2 * p.z * majorRadius) / normXZ;

        return (1.0 / (2 * sqTerm)) * Vector3{partialX, partialY, partialZ};
    }

    double ImplicitTorus::BoundingDiameter() {
        return majorRadius * 2 + minorRadius * 2;
    }

    Vector3 ImplicitTorus::BoundingCenter() {
        return center;
    }
    
    YZeroPlane::YZeroPlane() {}

    double YZeroPlane::SignedDistance(Vector3 point) {
        return point.y;
    }

    Vector3 YZeroPlane::GradientOfDistance(Vector3 point) {
        return Vector3{0, 1, 0};
    }

    double YZeroPlane::BoundingDiameter() {
        return 2;
    }

    Vector3 YZeroPlane::BoundingCenter() {
        return Vector3{0, 0.01, 0};
    }

    ImplicitDoubleTorus::ImplicitDoubleTorus(double r2) {
        radius2 = r2;
    }

    double ImplicitDoubleTorus::SignedDistance(Vector3 p) {
        double x = p.x + 2;
        double y = p.y;
        double z = p.z;
        double inner = (x * (x - 1) * (x - 1) * (x - 2) + y * y);
        return inner * inner + z * z - radius2;
    }

    Vector3 ImplicitDoubleTorus::GradientOfDistance(Vector3 point) {
        double x = point.x + 2;
        double y = point.y;
        double z = point.z;

        double x2 = x * x;
        double x3 = x2 * x;
        double x4 = x3 * x;
        double y2 = y * y;
        double polyn = x4 - 4 * x3 + 5 * x2 - 2 * x + y2;

        double partialX = 2 * polyn * (4 * x3 - 12 * x2 + 10 * x - 2);
        double partialY = 4 * y * polyn;
        double partialZ = 2 * z;

        return Vector3{partialX, partialY, partialZ};
    }

    double ImplicitDoubleTorus::BoundingDiameter() {
        return 4.5;
    }

    Vector3 ImplicitDoubleTorus::BoundingCenter() {
        return Vector3{0, 0, 0};
    }

}

