#pragma once

#include "geometrycentral/utilities/vector3.h"

namespace LWS {

    using namespace geometrycentral;

    class ImplicitSurface {
        public:
        virtual ~ImplicitSurface();
        virtual double SignedDistance(Vector3 point) = 0;
        virtual Vector3 GradientOfDistance(Vector3 point) = 0;
    };

    class ImplicitSphere : public ImplicitSurface {
        public:
        ImplicitSphere(double r);
        virtual double SignedDistance(Vector3 point);
        virtual Vector3 GradientOfDistance(Vector3 point);

        private:
        double radius;
    };
}