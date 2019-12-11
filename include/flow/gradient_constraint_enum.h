#pragma once

#include "gradient_constraints.h"

namespace LWS {
    enum class ConstraintType {
        Barycenter, EdgeLengths, TotalLength, Pins, TangentPins, Surface
    };

    class PolyCurveNetwork;
    int NumRowsForConstraint(ConstraintType type, PolyCurveNetwork* curve);
    std::string NameOfConstraint(ConstraintType type);

    class VariableConstraintSet : public GradientConstraints<VariableConstraintSet> {
        private:
        PolyCurveNetwork* curves;

        public:
        VariableConstraintSet(PolyCurveNetwork* c) {
            curves = c;
        }

        void AddTriplets(std::vector<Eigen::Triplet<double>> &triplets) const;
        int NumConstraintRows() const;
        int NumExpectedCols() const;
        void SetTargetValues(Eigen::VectorXd &targets) const;
        void NegativeConstraintValues(Eigen::VectorXd &b, Eigen::VectorXd &targets) const;

        int startIndexOfConstraint(ConstraintType type);
        int rowsOfConstraint(ConstraintType type);
    };
}