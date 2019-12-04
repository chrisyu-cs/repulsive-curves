#include "implicit_surface.h"

namespace LWS {

    ImplicitSurface::~ImplicitSurface() {}

    ImplicitSphere::ImplicitSphere(double r) {
        radius = r;
    }

    double ImplicitSphere::SignedDistance(Vector3 point) {
        return point.norm() - radius;
    }

    Vector3 ImplicitSphere::GradientOfDistance(Vector3 point) {
        return point.normalize();
    }

}

