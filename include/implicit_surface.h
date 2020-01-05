#pragma once

#include "geometrycentral/utilities/vector3.h"

namespace LWS {

    using namespace geometrycentral;

    class ImplicitSurface {
        public:
        virtual ~ImplicitSurface();
        virtual double SignedDistance(Vector3 point) = 0;
        virtual Vector3 GradientOfDistance(Vector3 point) = 0;
        virtual double BoundingDiameter() = 0;
        virtual Vector3 BoundingCenter() = 0;
    };

    class ImplicitSphere : public ImplicitSurface {
        public:
        ImplicitSphere(double r);
        virtual double SignedDistance(Vector3 point);
        virtual Vector3 GradientOfDistance(Vector3 point);
        virtual double BoundingDiameter();
        virtual Vector3 BoundingCenter();

        private:
        double radius;
        Vector3 center;
    };

    class ImplicitTorus : public ImplicitSurface {
        public:
        ImplicitTorus(double major, double minor);
        virtual double SignedDistance(Vector3 point);
        virtual Vector3 GradientOfDistance(Vector3 point);
        virtual double BoundingDiameter();
        virtual Vector3 BoundingCenter();

        private:
        double majorRadius;
        double minorRadius;
        Vector3 center;
    };

    class YZeroPlane : public ImplicitSurface {
        public:
        YZeroPlane();
        virtual double SignedDistance(Vector3 point);
        virtual Vector3 GradientOfDistance(Vector3 point);
        virtual double BoundingDiameter();
        virtual Vector3 BoundingCenter();
    };
}