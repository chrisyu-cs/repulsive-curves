#pragma once

#include <Eigen/Core>
#include "poly_curve_network.h"

namespace LWS {

    class CurvePotential {
        public:
        CurvePotential();
        virtual ~CurvePotential();
        virtual void AddGradient(PolyCurveNetwork* curves, Eigen::MatrixXd &gradient);
    };

    class TotalLengthPotential : public CurvePotential {
        public:
        TotalLengthPotential(double wt);
        virtual void AddGradient(PolyCurveNetwork* curves, Eigen::MatrixXd &gradient);

        private:
        double weight;
    };

}