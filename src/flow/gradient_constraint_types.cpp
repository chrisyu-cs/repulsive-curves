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
        int nPinRows = 3 * curves->pinnedVertices.size();
        // Plus 3 rows for barycenter, and |E| rows for edge lengths
        return curves->NumEdges() + nPinRows;
    }

    int LengthsAndPinsConstraint::NumExpectedCols() const {
        return curves->NumVertices() * 3;
    }

    void LengthsAndPinsConstraint::SetTargetValues(Eigen::VectorXd &targets) const {
        // Set target edge lengths to current lengths
        int nEdges = curves->NumEdges();
        for (int i = 0; i < nEdges; i++) {
            targets(i) = curves->GetEdge(i)->Length();
        }

        // Set pinned vertices to current positions
        int nPins = curves->pinnedVertices.size();
        int pinBase = nEdges;
        for (int i = 0; i < nPins; i++) {
            Vector3 p_i = curves->GetPinnedVertex(i)->Position();
            targets(pinBase + 3 * i    ) = p_i.x;
            targets(pinBase + 3 * i + 1) = p_i.y;
            targets(pinBase + 3 * i + 2) = p_i.z;
        }
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
        return 3 * curves->pinnedVertices.size();
    }

    int PinsConstraint::NumExpectedCols() const {
        return curves->NumVertices() * 3;
    }

    void PinsConstraint::SetTargetValues(Eigen::VectorXd &targets) const {
        // Set pinned vertices to current positions
        int nPins = curves->pinnedVertices.size();
        for (int i = 0; i < nPins; i++) {
            Vector3 p_i = curves->GetPinnedVertex(i)->Position();
            targets(3 * i    ) = p_i.x;
            targets(3 * i + 1) = p_i.y;
            targets(3 * i + 2) = p_i.z;
        }
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
        targets(0) = 0;
        targets(1) = 0;
        targets(2) = 0;

        // Set target edge lengths to current lengths
        int nEdges = curves->NumEdges();
        for (int i = 0; i < nEdges; i++) {
            targets(3 + i) = curves->GetEdge(i)->Length();
        }
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