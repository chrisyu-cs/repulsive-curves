#pragma once

#include <Eigen/Sparse>
#include "gradient_constraints.h"
#include "constraint_functions.h"

namespace LWS {

    class PolyCurveNetwork;

    class EdgeLengthConstraint : public GradientConstraints<EdgeLengthConstraint> {
        private:
        PolyCurveNetwork* curves;

        public:
        EdgeLengthConstraint(PolyCurveNetwork* c) {
            curves = c;
        }

        void AddTriplets(std::vector<Eigen::Triplet<double>> &triplets) const;
        int NumConstraintRows() const;
        int NumExpectedCols() const;
        void SetTargetValues(Eigen::VectorXd &targets) const;
        void NegativeConstraintValues(Eigen::VectorXd &b, Eigen::VectorXd &targets) const;
    };

    class BarycenterConstraint3X : public GradientConstraints<BarycenterConstraint3X> {
        private:
        PolyCurveNetwork* curves;

        public:
        BarycenterConstraint3X(PolyCurveNetwork* c) {
            curves = c;
        }

        void AddTriplets(std::vector<Eigen::Triplet<double>> &triplets) const;
        int NumConstraintRows() const;
        int NumExpectedCols() const;
        void SetTargetValues(Eigen::VectorXd &targets) const;
        void NegativeConstraintValues(Eigen::VectorXd &b, Eigen::VectorXd &targets) const;
    };

    class BarycenterConstraint : public GradientConstraints<BarycenterConstraint> {
        private:
        PolyCurveNetwork* curves;

        public:
        BarycenterConstraint(PolyCurveNetwork* c) {
            curves = c;
        }

        void AddTriplets(std::vector<Eigen::Triplet<double>> &triplets) const;
        int NumConstraintRows() const;
        int NumExpectedCols() const;
        void SetTargetValues(Eigen::VectorXd &targets) const;
        void NegativeConstraintValues(Eigen::VectorXd &b, Eigen::VectorXd &targets) const;
    };
}