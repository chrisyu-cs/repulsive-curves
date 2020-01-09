#include "implicit_surface.h"
#include "utils.h"

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

    ImplicitUnion::ImplicitUnion(ImplicitSurface* s1, ImplicitSurface* s2) {
        surfaces.push_back(s1);
        surfaces.push_back(s2);
    }

    double ImplicitUnion::SignedDistance(Vector3 p) {
        double minDist = surfaces[0]->SignedDistance(p);
        for (size_t i = 1; i < surfaces.size(); i++) {
            minDist = fmin(minDist, surfaces[i]->SignedDistance(p));
        }
        return minDist;
    }

    Vector3 ImplicitUnion::GradientOfDistance(Vector3 point) {
        double minDist = surfaces[0]->SignedDistance(point);
        double min_i = 0;

        for (size_t i = 1; i < surfaces.size(); i++) {
            double dist = surfaces[i]->SignedDistance(point);
            if (dist < minDist) {
                minDist = dist;
                min_i = i;
            }
        }

        return surfaces[min_i]->GradientOfDistance(point);
    }

    double ImplicitUnion::BoundingDiameter() {
        Vector3 center = BoundingCenter();
        double maxRadius = 0;
        
        for (size_t i = 0; i < surfaces.size(); i++) {
            Vector3 c_i = surfaces[i]->BoundingCenter();
            double d_i = surfaces[i]->BoundingDiameter();
            Vector3 disp = (center - c_i);
            double l1dist = fmax(fabs(disp.x), fmax(fabs(disp.y), fabs(disp.z)));

            maxRadius = fmax(maxRadius, l1dist + d_i / 2);
        }
        return maxRadius * 2;
    }

    Vector3 ImplicitUnion::BoundingCenter() {
        Vector3 mins = surfaces[0]->BoundingCenter();
        Vector3 maxs = mins;

        for (size_t i = 1; i < surfaces.size(); i++) {
            mins = vector_min(mins, surfaces[i]->BoundingCenter());
            maxs = vector_max(maxs, surfaces[i]->BoundingCenter());
        }

        return (mins + maxs) / 2;
    }

    ImplicitSmoothUnion::ImplicitSmoothUnion(ImplicitSurface* s1, ImplicitSurface* s2, double blendFactor) {
        surface1 = s1;
        surface2 = s2;
        k = blendFactor;
    }

    double ImplicitSmoothUnion::SignedDistance(Vector3 p) {
        double d1 = surface1->SignedDistance(p);
        double d2 = surface2->SignedDistance(p);
        float h = clamp(0.5 + 0.5 * (d2 - d1) / k, 0.0, 1.0);
        return ((1 - h) * d2 + h * d1) - k * h * (1 - h);
    }


    Vector3 ImplicitSmoothUnion::GradientOfDistance(Vector3 point) {
        double d1 = surface1->SignedDistance(point);
        double d2 = surface2->SignedDistance(point);
        float h = clamp(0.5 + 0.5 * (d2 - d1) / k, 0.0, 1.0);

        Vector3 deriv_d1 = surface1->GradientOfDistance(point);
        Vector3 deriv_d2 = surface2->GradientOfDistance(point);
        // We need derivative of a clamp, which is just 0 outside of (0, 1)
        Vector3 deriv_h{0, 0, 0};
        if (h > 0 && h < 1) {
            deriv_h = (1.0 / (2 * k)) * (deriv_d2 - deriv_d1);
        }

        Vector3 deriv_term1 = (-deriv_h * d2 + (1 - h) * deriv_d2);
        Vector3 deriv_term2 = (deriv_h * d1 + h * deriv_d1);
        Vector3 deriv_term3 = -k * (deriv_h - 2 * h * deriv_h);

        return deriv_term1 + deriv_term2 + deriv_term3;
    }

    double ImplicitSmoothUnion::BoundingDiameter() {
        Vector3 center = BoundingCenter();
        double maxRadius = 0;

        ImplicitSurface* surfaces[2] = {surface1, surface2};
        
        for (size_t i = 0; i < 2; i++) {
            Vector3 c_i = surfaces[i]->BoundingCenter();
            double d_i = surfaces[i]->BoundingDiameter();
            Vector3 disp = (center - c_i);
            double l1dist = fmax(fabs(disp.x), fmax(fabs(disp.y), fabs(disp.z)));

            maxRadius = fmax(maxRadius, l1dist + d_i / 2);
        }
        return maxRadius * 2;
    }

    Vector3 ImplicitSmoothUnion::BoundingCenter() {
        Vector3 mins = surface1->BoundingCenter();
        Vector3 maxs = mins;
        mins = vector_min(mins, surface2->BoundingCenter());
        maxs = vector_max(maxs, surface2->BoundingCenter());

        return (mins + maxs) / 2;
    }


}

