#pragma once

#include "tpe_energy_sc.h"
#include "sobo_slobo.h"
#include "spatial/tpe_bvh.h"
#include "product/block_cluster_tree.h"

#include "libgmultigrid/matrix_free.h"
#include "libgmultigrid/multigrid_hierarchy.h"
#include "multigrid/constraint_projector_domain.h"
#include "flow/gradient_constraint_enum.h"

#include "obstacles/obstacle.h"
#include "extra_potentials.h"
#include <iostream>

namespace LWS {

    struct CoordinateLUs {
        Eigen::PartialPivLU<Eigen::MatrixXd> lu_x;
        Eigen::PartialPivLU<Eigen::MatrixXd> lu_y;
        Eigen::PartialPivLU<Eigen::MatrixXd> lu_z;
    };


    class TPEFlowSolverSC {
        public:
        using ConstraintClassType = VariableConstraintSet;
        std::vector<Obstacle*> obstacles;
        std::vector<CurvePotential*> potentials;

        TPEFlowSolverSC(PolyCurveNetwork* p, double a, double b);
        ~TPEFlowSolverSC();
        ConstraintClassType constraint;

        void ReplaceCurve(PolyCurveNetwork* new_p);
        void EnablePerformanceLog(std::string logFile);
        void ClosePerformanceLog();

        void UpdateTargetLengths();
        void SetTotalLengthScaleTarget(double scale);
        void SetEdgeLengthScaleTarget(double scale);
        void MoveLengthTowardsTarget();
        bool TargetLengthReached();

        double CurrentEnergy(SpatialTree *root = 0);
        double TPEnergyDirect();
        double TPEnergyBH(SpatialTree *root);
        void FillGradientSingle(Eigen::MatrixXd &gradients, int i, int j);
        void FillGradientVectorDirect(Eigen::MatrixXd &gradients);
        void FillGradientVectorBH(SpatialTree *root, Eigen::MatrixXd &gradients);
        void AddAllGradients(SpatialTree* root, Eigen::MatrixXd &gradients);
        
        inline void SetExponents(double a, double b) {
            alpha = a;
            beta = b;
        }

        bool StepNaive(double h);
        bool StepLS(bool useBH);
        bool StepLSConstrained(bool useBH, bool useBackproj);
        bool StepSobolevLS(bool useBH, bool useBackproj);
        bool StepSobolevLSIterative(double epsilon, bool useBackproj);

        double ProjectGradient(Eigen::MatrixXd &gradients, Eigen::MatrixXd &A, Eigen::PartialPivLU<Eigen::MatrixXd> &lu);
        double ProjectSoboSloboGradient(Eigen::PartialPivLU<Eigen::MatrixXd> &lu, Eigen::MatrixXd &gradients);
        void GetSecondDerivative(SpatialTree* tree_root, Eigen::MatrixXd &projected1, double epsilon, Eigen::MatrixXd &secondDeriv);

        template<typename Domain, typename Smoother>
        inline void BackprojectMultigrid(MultigridHierarchy<Domain>* solver, PolyCurveNetwork* curveNetwork, Eigen::VectorXd &phi, Eigen::MatrixXd &output, double tol) {
            int nRows = solver->NumRows();
            // Flatten the gradient matrix into a long vector
            Eigen::VectorXd B_pinv_phi(nRows);
            curveNetwork->constraintProjector->ApplyBPinv(phi, B_pinv_phi);
            Product::MatrixReplacement<typename Domain::MultType> multiplier(solver->GetTopLevelMultiplier(), nRows);
            Eigen::VectorXd GB_phi = multiplier * B_pinv_phi;
            // Solve Gv = b by solving PGPv = Pb
            GB_phi = curveNetwork->constraintProjector->ProjectToNullspace(GB_phi);
            Eigen::VectorXd v = solver->template VCycleSolve<Smoother>(GB_phi, tol);
            v = B_pinv_phi - v;
            VectorXdIntoMatrix(v, output);
        }

        template<typename Domain, typename Smoother>
        double ProjectGradientMultigrid(Eigen::MatrixXd &gradients, MultigridHierarchy<Domain>* solver, Eigen::MatrixXd &output, double tol);
        template<typename Domain, typename Smoother>
        double LSBackprojectMultigrid(Eigen::MatrixXd &gradient, double initGuess, MultigridHierarchy<Domain>* solver, BVHNode3D* root, double tol);

        void SaveCurrentPositions();
        void RestoreOriginalPositions();
        double LineSearchStep(Eigen::MatrixXd &gradients, double gradDot = 1, BVHNode3D* root = 0, bool resetStep = false);
        double LineSearchStep(Eigen::MatrixXd &gradients, double initGuess, int doublingLimit, double gradDot, BVHNode3D* root);
        double CircleSearchStep(Eigen::MatrixXd &P_dot, Eigen::MatrixXd &P_ddot, Eigen::MatrixXd &G, BVHNode3D* root);

        double LSBackproject(Eigen::MatrixXd &gradients, double initGuess,
            Eigen::PartialPivLU<Eigen::MatrixXd> &lu, double gradDot, BVHNode3D* root);

        template<typename Domain, typename Smoother>
        double BackprojectConstraintsMultigrid(Eigen::MatrixXd &gradient, MultigridHierarchy<Domain>* solver, double tol);

        inline bool PerformanceLogEnabled() {
            return perfLogEnabled;
        }

        bool soboNormZero;

        private:
        bool perfLogEnabled;
        std::ofstream perfFile;
        bool useEdgeLengthScale;
        bool useTotalLengthScale;
        double lengthScaleStep;
        int iterNum;
        double targetLength;
        double ls_step_threshold;
        double mg_tolerance;
        double backproj_threshold;
        double mg_backproj_threshold;
        double lastStepSize;
        PolyCurveNetwork* curveNetwork;
        Eigen::MatrixXd originalPositionMatrix;
        Eigen::VectorXd constraintTargets;
        Eigen::VectorXd fullDerivVector;
        double alpha;
        double beta;
        void SetGradientStep(Eigen::MatrixXd &gradient, double delta);
        void SetCircleStep(Eigen::MatrixXd &P_dot, Eigen::MatrixXd &K, double sqrt_G, double R, double alpha_delta);
        double BackprojectConstraints(Eigen::PartialPivLU<Eigen::MatrixXd> &lu);
    };

    template<typename Domain, typename Smoother>
    double TPEFlowSolverSC::ProjectGradientMultigrid(Eigen::MatrixXd &gradients, MultigridHierarchy<Domain>* solver, Eigen::MatrixXd &output, double tol) {
        // Flatten the gradient matrix into a long vector
        Eigen::VectorXd gradients3x;
        // Copy only the rows that actually contain gradient vectors
        gradients3x.setZero(curveNetwork->NumVertices() * 3);
        auto block = gradients.block(0, 0, curveNetwork->NumVertices(), 3);
        MatrixIntoVectorX3(block, gradients3x);
        // Project this vector into the constraint null-space:
        // we really want to solve PGPx = Pb
        gradients3x = curveNetwork->constraintProjector->ProjectToNullspace(gradients3x);
        // Solve PGPx = Pb using multigrid
        Eigen::VectorXd sobolevGradients = solver->template VCycleSolve<Smoother>(gradients3x, tol);
        // Compute dot product with unprojected gradient, and copy into results vector
        double soboDot = gradients3x.dot(sobolevGradients); 
        // double dirDot = soboDot / (gradients3x.norm() * sobolevGradients.norm());
        output.setZero();
        VectorXdIntoMatrix(sobolevGradients, output);
        return soboDot;
    }

    template<typename Domain, typename Smoother>
    double TPEFlowSolverSC::BackprojectConstraintsMultigrid(Eigen::MatrixXd &gradient, MultigridHierarchy<Domain>* solver, double tol) {
        int nVerts = curveNetwork->NumVertices();
        Eigen::VectorXd phi(constraint.NumConstraintRows());
        Eigen::MatrixXd correction(nVerts, 3);
        correction.setZero();
        constraint.FillConstraintValues(phi, constraintTargets, 0);
        // Compute and apply the correction
        this->template BackprojectMultigrid<Domain, Smoother>(solver, curveNetwork, phi, correction, tol);
        for (int i = 0; i < nVerts; i++) {
            CurveVertex* p = curveNetwork->GetVertex(i);
            Vector3 cur = p->Position();
            p->SetPosition(cur + SelectRow(correction, i));
        }
        // Add length violations to RHS
        double maxViolation = constraint.FillConstraintValues(phi, constraintTargets, 0);
        std::cout << "  Constraint value = " << maxViolation << std::endl;
        return maxViolation;
    }

    template<typename Domain, typename Smoother>
    double TPEFlowSolverSC::LSBackprojectMultigrid(Eigen::MatrixXd &gradient, double initGuess,
    MultigridHierarchy<Domain>* solver, BVHNode3D* root, double tol) {
        double delta = initGuess;
        int attempts = 0;

        while ((delta > ls_step_threshold || useEdgeLengthScale) && attempts < 10) {
            attempts++;
            SetGradientStep(gradient, delta);
            if (root) {
                // Update the centers of mass to reflect the new positions
                root->recomputeCentersOfMass(curveNetwork);
            }

            for (int c = 0; c < 2; c++) {
                double maxViolation = BackprojectConstraintsMultigrid<Domain, Smoother>(gradient, solver, tol);
                if (maxViolation < mg_backproj_threshold) {
                    std::cout << "Backprojection successful after " << attempts << " attempts" << std::endl;
                    std::cout << "Used " << (c + 1) << " Newton steps on successful attempt" << std::endl;
                    return delta;
                }
            }

            delta /= 2;
        }
        std::cout << "Couldn't make backprojection succeed after " << attempts << " attempts (initial step " << initGuess << ")" << std::endl;
        BackprojectConstraintsMultigrid<Domain, Smoother>(gradient, solver, tol);
        return delta;
    }
}
