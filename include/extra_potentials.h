#pragma once

#include <Eigen/Core>
#include "poly_curve_network.h"
#include "vert_jacobian.h"

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

    class PinBendingPotential : public CurvePotential {
        public:
        PinBendingPotential(double wt);
        virtual double CurrentValue(PolyCurveNetwork* curves);
        virtual void AddGradient(PolyCurveNetwork* curves, Eigen::MatrixXd &gradient);

        private:
        double weight;
    };

    class VectorField {
        public:
        VectorField();
        virtual ~VectorField();
        virtual Vector3 Sample(Vector3 x);
        virtual VertJacobian SpatialDerivative(Vector3 x);
    };

    class ConstantVectorField : public VectorField {
        public:
        ConstantVectorField(Vector3 v);
        virtual Vector3 Sample(Vector3 x);
        virtual VertJacobian SpatialDerivative(Vector3 x);
        private:
        Vector3 c;
    };

    class CircularVectorField : public VectorField {
        public:
        CircularVectorField();
        virtual Vector3 Sample(Vector3 x);
        virtual VertJacobian SpatialDerivative(Vector3 x);
        private:
        Vector3 dirDeriv(Vector3 x, Vector3 dir);
    };

    class InterestingVectorField : public VectorField {
        public:
        InterestingVectorField();
        virtual Vector3 Sample(Vector3 x);
        virtual VertJacobian SpatialDerivative(Vector3 x);
    };

    class VectorFieldPotential : public CurvePotential {
        public:
        VectorFieldPotential(double wt, VectorField* vf);
        ~VectorFieldPotential();
        virtual double CurrentValue(PolyCurveNetwork* curves);
        virtual void AddGradient(PolyCurveNetwork* curves, Eigen::MatrixXd &gradient);
        private:
        double weight;
        VectorField* field;
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