#pragma once

#include <Eigen/Core>
#include "poly_curve_network.h"

namespace LWS {

    class CurvePotential {
        public:
        CurvePotential();
        virtual ~CurvePotential();
        virtual double CurrentValue(PolyCurveNetwork* curves);
        virtual void AddGradient(PolyCurveNetwork* curves, Eigen::MatrixXd &gradient);
    };

    class TotalLengthPotential : public CurvePotential {
        public:
        TotalLengthPotential(double wt);
        virtual double CurrentValue(PolyCurveNetwork* curves);
        virtual void AddGradient(PolyCurveNetwork* curves, Eigen::MatrixXd &gradient);

        private:
        double weight;
    };

    class LengthDifferencePotential : public CurvePotential {
        public:
        LengthDifferencePotential(double wt);
        virtual double CurrentValue(PolyCurveNetwork* curves);
        virtual void AddGradient(PolyCurveNetwork* curves, Eigen::MatrixXd &gradient);

        private:
        double LenDiff(PolyCurveNetwork* curves, int i);
        double weight;
    };

    // class AreaPotential : public CurvePotential {
    //     public:
    //     AreaPotential(double wt);
    //     virtual double CurrentValue(PolyCurveNetwork* curves);
    //     virtual void AddGradient(PolyCurveNetwork* curves, Eigen::MatrixXd &gradient);

    //     private:
    //     double weight;  
    // };

}