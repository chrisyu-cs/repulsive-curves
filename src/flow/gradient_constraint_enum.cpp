#include "flow/gradient_constraint_enum.h"

#include "poly_curve_network.h"
#include "flow/constraint_functions.h"

namespace LWS {

    int NumRowsForConstraint(ConstraintType type, PolyCurveNetwork* curve) {
        switch (type) {
            case ConstraintType::Barycenter:
            return 3;
            case ConstraintType::EdgeLengths:
            return curve->NumEdges();
            case ConstraintType::Pins:
            return curve->NumPins();
            default:
            std::cerr << "Called NumRowsForConstraint on an unimplemented constraint type" << std::endl;
            return 0;
        }
    }
        
    void AddTripletsOfConstraint(ConstraintType type, PolyCurveNetwork* curve,
    std::vector<Eigen::Triplet<double>> &triplets, int start) {
        switch (type) {
            case ConstraintType::Barycenter:
            ConstraintFunctions::AddBarycenterTriplets3X(curve, triplets, start);
            break;
            case ConstraintType::EdgeLengths:
            ConstraintFunctions::AddEdgeLengthTriplets(curve, triplets, start);
            break;
            case ConstraintType::Pins:
            ConstraintFunctions::AddPinTriplets(curve, triplets, start);
            break;
            default:
            std::cerr << "Called AddTripletsOfConstraint on an unimplemented constraint type" << std::endl;
            break;
        }
    }

    void SetTargetsOfConstraint(ConstraintType type, PolyCurveNetwork* curve, Eigen::VectorXd &targets, int start) {
        switch (type) {
            case ConstraintType::Barycenter:
            ConstraintFunctions::SetBarycenterTargets(curve, targets, start);
            break;
            case ConstraintType::EdgeLengths:
            ConstraintFunctions::SetEdgeLengthTargets(curve, targets, start);
            break;
            case ConstraintType::Pins:
            ConstraintFunctions::SetPinTargets(curve, targets, start);
            break;
            default:
            std::cerr << "Called SetTargetsOfConstraint on an unimplemented constraint type" << std::endl;
            break;
        }
    }

    void NegativeViolationsOfConstraint(ConstraintType type, PolyCurveNetwork* curve,
    Eigen::VectorXd &b, Eigen::VectorXd &targets, int start) {
        switch (type) {
            case ConstraintType::Barycenter:
            ConstraintFunctions::NegativeBarycenterViolation(curve, b, targets, start);
            break;
            case ConstraintType::EdgeLengths:
            ConstraintFunctions::NegativeEdgeLengthViolation(curve, b, targets, start);
            break;
            case ConstraintType::Pins:
            ConstraintFunctions::NegativePinViolation(curve, b, targets, start);
            break;
            default:
            std::cerr << "Called NegativeViolationsOfConstraint on an unimplemented constraint type" << std::endl;
            break;
        }
    }

    void VariableConstraintSet::AddTriplets(std::vector<Eigen::Triplet<double>> &triplets) const {
        int startIndex = 0;
        for (ConstraintType &type : curves->appliedConstraints) {
            AddTripletsOfConstraint(type, curves, triplets, startIndex);
            startIndex += NumRowsForConstraint(type, curves);
        }
    }

    int VariableConstraintSet::NumConstraintRows() const {
        int rows = 0;
        for (ConstraintType &type : curves->appliedConstraints) {
            rows += NumRowsForConstraint(type, curves);
        }
        return rows;
    }

    int VariableConstraintSet::NumExpectedCols() const {
        return curves->NumVertices() * 3;
    }
    
    void VariableConstraintSet::SetTargetValues(Eigen::VectorXd &targets) const {
        int startIndex = 0;
        for (ConstraintType &type : curves->appliedConstraints) {
            SetTargetsOfConstraint(type, curves, targets, startIndex);
            startIndex += NumRowsForConstraint(type, curves);
        }

    }

    void VariableConstraintSet::NegativeConstraintValues(Eigen::VectorXd &b, Eigen::VectorXd &targets) const {
        int startIndex = 0;
        for (ConstraintType &type : curves->appliedConstraints) {
            NegativeViolationsOfConstraint(type, curves, b, targets, startIndex);
            startIndex += NumRowsForConstraint(type, curves);
        }
    }

}
