#include "implicit_surface.h"

namespace LWS {

    ImplicitSurface::~ImplicitSurface() {}

    ImplicitSphere::ImplicitSphere(double r) {
        radius = r;
        center = Vector3{0, 0, 0};
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

}

