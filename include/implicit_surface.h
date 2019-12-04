#pragma once

#include "poly_curve_network.h"

namespace LWS {

    class ImplicitSurface {
        public:
        virtual ~ImplicitSurface();
        virtual double SignedDistance(Vector3 point);
        virtual Vector3 GradientOfDistance(Vector3 point);
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