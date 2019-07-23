#pragma once

#include "tpe_energy_sc.h"
#include "sobo_slobo.h"
#include "mesh_helpers.h"
#include "poly_curve.h"
#include "spatial/tpe_bvh.h"
#include "spatial/block_cluster_tree.h"

namespace LWS {

    struct CoordinateLUs {
        Eigen::PartialPivLU<Eigen::MatrixXd> lu_x;
        Eigen::PartialPivLU<Eigen::MatrixXd> lu_y;
        Eigen::PartialPivLU<Eigen::MatrixXd> lu_z;
    };

    class TPEFlowSolverSC {
        public:
        TPEFlowSolverSC(PolyCurveGroup* p);
        double CurrentEnergy(SpatialTree *root = 0);
        inline double CurrentEnergyDirect();
        double CurrentEnergyBH(SpatialTree *root);
        void FillGradientSingle(std::vector<Vector3> &gradients, int i, int j);
        void FillGradientVectorDirect(std::vector<Vector3> &gradients);

        void FillGradientVectorBH(SpatialTree *root, std::vector<Vector3> &gradients);

        void FillConstraintVector(std::vector<Vector3> &gradients);
        bool StepBoundaryNaive(double h);
        bool StepBoundaryLS();
        bool StepSobolevProjLS(bool useBH);

        double ComputeAndProjectGradient(std::vector<Vector3> &gradients);
        double ComputeAndProjectGradient(std::vector<Vector3> &gradients, Eigen::MatrixXd &A, Eigen::PartialPivLU<Eigen::MatrixXd> &lu);
        double ProjectSoboSloboGradient(Eigen::PartialPivLU<Eigen::MatrixXd> &lu, std::vector<Vector3> &gradients);
        
        // Fill a Sobolev-Slobodeckij inner product matrix with a barycenter constraint
        void SoboSloboMatrix(Eigen::MatrixXd &A);
        void FillVertLengthConstraintMatrix(Eigen::MatrixXd &A, int baseIndex);

        void CompareMatrixVectorProduct();

        void ExpandMatrix3x(Eigen::MatrixXd &A, Eigen::MatrixXd &B);

        void SaveCurrentPositions();
        void RestoreOriginalPositions();
        double LineSearchStep(std::vector<Vector3> &gradients, double gradDot = 1, SpatialTree* root = 0);
        double LineSearchStep(std::vector<Vector3> &gradients, double initGuess, int doublingLimit, double gradDot, SpatialTree* root);
        double LSBackproject(std::vector<Vector3> &gradients, double initGuess,
            Eigen::PartialPivLU<Eigen::MatrixXd> &lu, double gradDot, SpatialTree* root);

        bool useEdgeLengthConstraint;
        int matrixNumRows();

        private:
        double ls_step_threshold;
        double backproj_threshold;
        CoordinateLUs coord_lus;
        PolyCurveGroup* curves;
        std::vector<Vector3> originalPositions;
        std::vector<Vector3> l2gradients;
        std::vector<Vector3> vertGradients;
        std::vector<Vector3> vertConstraints;
        std::vector<double> initialLengths;
        double alpha;
        double beta;
        void SetGradientStep(std::vector<Vector3> gradient, double delta);
        double FillLengthConstraintViolations(Eigen::VectorXd &b, int baseIndex);
        bool BackprojectConstraints(Eigen::PartialPivLU<Eigen::MatrixXd> &lu);
    };
}
