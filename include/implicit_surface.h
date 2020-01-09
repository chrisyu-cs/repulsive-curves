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
        ImplicitSphere(double r, Vector3 c);
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
        ImplicitTorus(double major, double minor, Vector3 c);
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

    class ImplicitDoubleTorus : public ImplicitSurface {
        public:
        ImplicitDoubleTorus(double r2);
        virtual double SignedDistance(Vector3 point);
        virtual Vector3 GradientOfDistance(Vector3 point);
        virtual double BoundingDiameter();
        virtual Vector3 BoundingCenter();

        private:
        double radius2;
    };

    class ImplicitUnion : public ImplicitSurface {
        public:
        ImplicitUnion(ImplicitSurface* s1, ImplicitSurface* s2);
        virtual double SignedDistance(Vector3 point);
        virtual Vector3 GradientOfDistance(Vector3 point);
        virtual double BoundingDiameter();
        virtual Vector3 BoundingCenter();

        private:
        std::vector<ImplicitSurface*> surfaces;
    };

    class ImplicitSmoothUnion : public ImplicitSurface {
        public:
        ImplicitSmoothUnion(ImplicitSurface* s1, ImplicitSurface* s2, double blendFactor);
        virtual double SignedDistance(Vector3 point);
        virtual Vector3 GradientOfDistance(Vector3 point);
        virtual double BoundingDiameter();
        virtual Vector3 BoundingCenter();
        
        private:
        ImplicitSurface* surface1;
        ImplicitSurface* surface2;
        double k;
    };
}