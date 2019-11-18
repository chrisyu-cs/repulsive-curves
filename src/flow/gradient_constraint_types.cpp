#include "flow/gradient_constraint_types.h"
#include "poly_curve_network.h"

namespace LWS {

    void LengthsAndPinsConstraint::AddTriplets(std::vector<Eigen::Triplet<double>> &triplets) const {
        int nEdges = curves->NumEdges();

        // Add the edge length rows
        ConstraintFunctions::AddEdgeLengthTriplets(curves, triplets, 0);
        // Add rows for pinned vertices
        ConstraintFunctions::AddPinTriplets(curves, triplets, nEdges);
    }

    int LengthsAndPinsConstraint::NumConstraintRows() const {
        // Each pinned vertex gets 3 rows
        int nPinRows = 3 * curves->NumPins();
        // Plus 3 rows for barycenter, and |E| rows for edge lengths
        return curves->NumEdges() + nPinRows;
    }

    int LengthsAndPinsConstraint::NumExpectedCols() const {
        return curves->NumVertices() * 3;
    }

    void LengthsAndPinsConstraint::SetTargetValues(Eigen::VectorXd &targets) const {
        int nEdges = curves->NumEdges();
        // Set target edge lengths to current lengths
        ConstraintFunctions::SetEdgeLengthTargets(curves, targets, 0);
        // Pin targets to current positions
        ConstraintFunctions::SetPinTargets(curves, targets, nEdges);
    }

    void LengthsAndPinsConstraint::NegativeConstraintValues(Eigen::VectorXd &b, Eigen::VectorXd &targets) const {
        ConstraintFunctions::NegativeEdgeLengthViolation(curves, b, targets, 0);
        ConstraintFunctions::NegativePinViolation(curves, b, targets, curves->NumEdges());
    }

    void PinsConstraint::AddTriplets(std::vector<Eigen::Triplet<double>> &triplets) const {
        int nEdges = curves->NumEdges();

        // Add rows for pinned vertices
        ConstraintFunctions::AddPinTriplets(curves, triplets, 0);
    }

    int PinsConstraint::NumConstraintRows() const {
        // Each pinned vertex gets 3 rows
        return 3 * curves->NumPins();
    }

    int PinsConstraint::NumExpectedCols() const {
        return curves->NumVertices() * 3;
    }

    void PinsConstraint::SetTargetValues(Eigen::VectorXd &targets) const {
        ConstraintFunctions::SetPinTargets(curves, targets, 0);
    }

    void PinsConstraint::NegativeConstraintValues(Eigen::VectorXd &b, Eigen::VectorXd &targets) const {
        ConstraintFunctions::NegativePinViolation(curves, b, targets, 0);
    }

    void EdgeLengthConstraint::AddTriplets(std::vector<Eigen::Triplet<double>> &triplets) const {
        // Add barycenter constraints
        ConstraintFunctions::AddBarycenterTriplets3X(curves, triplets, 0);
        // Add the edge length rows
        ConstraintFunctions::AddEdgeLengthTriplets(curves, triplets, 3);
    }

    int EdgeLengthConstraint::NumConstraintRows() const {
        // 3 rows for barycenter, and |E| rows for edge lengths
        return curves->NumEdges() + 3;
    }

    int EdgeLengthConstraint::NumExpectedCols() const {
        return curves->NumVertices() * 3;
    }

    void EdgeLengthConstraint::SetTargetValues(Eigen::VectorXd &targets) const {
        // Fix barycenter at origin
        ConstraintFunctions::SetBarycenterTargets(curves, targets, 0);
        // Set target edge lengths to current lengths
        ConstraintFunctions::SetEdgeLengthTargets(curves, targets, 3);
    }

    void EdgeLengthConstraint::NegativeConstraintValues(Eigen::VectorXd &b, Eigen::VectorXd &targets) const {
        ConstraintFunctions::NegativeBarycenterViolation(curves, b, targets, 0);
        ConstraintFunctions::NegativeEdgeLengthViolation(curves, b, targets, 3);
    }

    void BarycenterConstraint3X::AddTriplets(std::vector<Eigen::Triplet<double>> &triplets) const {
        // Add barycenter constraints
        ConstraintFunctions::AddBarycenterTriplets3X(curves, triplets, 0);
    }

    int BarycenterConstraint3X::NumConstraintRows() const {
        // 3 rows for barycenter
        return 3;
    }

    int BarycenterConstraint3X::NumExpectedCols() const {
        return curves->NumVertices() * 3;
    }

    void BarycenterConstraint3X::SetTargetValues(Eigen::VectorXd &targets) const {
        // Fix barycenter at origin
        targets(0) = 0;
        targets(1) = 0;
        targets(2) = 0;
    }

    void BarycenterConstraint3X::NegativeConstraintValues(Eigen::VectorXd &b, Eigen::VectorXd &targets) const {
        ConstraintFunctions::NegativeBarycenterViolation(curves, b, targets, 0);
    }

    void BarycenterConstraint::AddTriplets(std::vector<Eigen::Triplet<double>> &triplets) const {
        std::cerr << "1D barycenter constraint currently disabled" << std::endl;
        throw 1;
    }

    int BarycenterConstraint::NumExpectedCols() const {
        return curves->NumVertices();
    }

    int BarycenterConstraint::NumConstraintRows() const {
        // 1 row for barycenter
        return 1;
    }

    void BarycenterConstraint::SetTargetValues(Eigen::VectorXd &targets) const {
        // Fix barycenter at origin
        targets(0) = 0;
    }

    void BarycenterConstraint::NegativeConstraintValues(Eigen::VectorXd &b, Eigen::VectorXd &targets) const {
        std::cerr << "1D barycenter constraint currently disabled" << std::endl;
        throw 1;
    }

    
}