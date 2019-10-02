#pragma once

#include "tpe_energy_sc.h"
#include "sobo_slobo.h"
#include "mesh_helpers.h"
#include "poly_curve.h"
#include "spatial/tpe_bvh.h"
#include "product/block_cluster_tree.h"

#include "product/matrix_free.h"

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
        void FillGradientSingle(Eigen::MatrixXd &gradients, int i, int j);
        void FillGradientVectorDirect(Eigen::MatrixXd &gradients);

        void FillGradientVectorBH(SpatialTree *root, Eigen::MatrixXd &gradients);

        void FillConstraintVector(Eigen::MatrixXd &gradients);
        bool StepNaive(double h);
        bool StepLS();
        bool StepSobolevLS(bool useBH);
        bool StepSobolevLSIterative();

        double ComputeAndProjectGradient(Eigen::MatrixXd &gradients);
        double ComputeAndProjectGradient(Eigen::MatrixXd &gradients, Eigen::MatrixXd &A, Eigen::PartialPivLU<Eigen::MatrixXd> &lu);
        double ProjectSoboSloboGradient(Eigen::PartialPivLU<Eigen::MatrixXd> &lu, Eigen::MatrixXd &gradients);
        
        void ExpandMatrix3x(Eigen::MatrixXd &A, Eigen::MatrixXd &B);

        void SaveCurrentPositions();
        void RestoreOriginalPositions();
        double LineSearchStep(Eigen::MatrixXd &gradients, double gradDot = 1, SpatialTree* root = 0);
        double LineSearchStep(Eigen::MatrixXd &gradients, double initGuess, int doublingLimit, double gradDot, SpatialTree* root);
        double LSBackproject(Eigen::MatrixXd &gradients, double initGuess,
            Eigen::PartialPivLU<Eigen::MatrixXd> &lu, double gradDot, SpatialTree* root);

        bool useEdgeLengthConstraint;
        int matrixNumRows();

        private:
        double ls_step_threshold;
        double backproj_threshold;
        CoordinateLUs coord_lus;
        PolyCurveGroup* curves;
        std::vector<Vector3> originalPositions;
        Eigen::MatrixXd l2gradients;
        Eigen::MatrixXd vertGradients;
        Eigen::MatrixXd vertConstraints;
        std::vector<double> initialLengths;
        double alpha;
        double beta;
        void SetGradientStep(Eigen::MatrixXd gradient, double delta);
        double FillLengthConstraintViolations(Eigen::VectorXd &b, int baseIndex);
        bool BackprojectConstraints(Eigen::PartialPivLU<Eigen::MatrixXd> &lu);
    };
}
