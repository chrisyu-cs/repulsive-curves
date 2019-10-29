#pragma once

#include "tpe_energy_sc.h"
#include "sobo_slobo.h"
#include "spatial/tpe_bvh.h"
#include "product/block_cluster_tree.h"

#include "product/matrix_free.h"
#include "multigrid/multigrid_hierarchy.h"

namespace LWS {

    struct CoordinateLUs {
        Eigen::PartialPivLU<Eigen::MatrixXd> lu_x;
        Eigen::PartialPivLU<Eigen::MatrixXd> lu_y;
        Eigen::PartialPivLU<Eigen::MatrixXd> lu_z;
    };

    class TPEFlowSolverSC {
        public:
        TPEFlowSolverSC(PolyCurveNetwork* p, double a, double b);
        double CurrentEnergy(SpatialTree *root = 0);
        double CurrentEnergyDirect();
        double CurrentEnergyBH(SpatialTree *root);
        void FillGradientSingle(Eigen::MatrixXd &gradients, int i, int j);
        void FillGradientVectorDirect(Eigen::MatrixXd &gradients);

        void FillGradientVectorBH(SpatialTree *root, Eigen::MatrixXd &gradients);
        inline void SetExponents(double a, double b) {
            alpha = a;
            beta = b;
        }

        void FillConstraintVector(Eigen::MatrixXd &gradients);
        bool StepNaive(double h);
        bool StepLS();
        bool StepSobolevLS(bool useBH);
        bool StepSobolevLSIterative(double epsilon);

        double ComputeAndProjectGradient(Eigen::MatrixXd &gradients);
        double ComputeAndProjectGradient(Eigen::MatrixXd &gradients, Eigen::MatrixXd &A, Eigen::PartialPivLU<Eigen::MatrixXd> &lu);
        double ProjectSoboSloboGradient(Eigen::PartialPivLU<Eigen::MatrixXd> &lu, Eigen::MatrixXd &gradients);

        template<typename Domain, typename Smoother>
        double ProjectGradientMultigrid(Eigen::MatrixXd &gradients, MultigridHierarchy<Domain>* solver, Eigen::MatrixXd &output);
        template<typename Domain, typename Smoother>
        bool BackprojectMultigridLS(Eigen::MatrixXd &gradient, double initGuess, MultigridHierarchy<Domain>* solver, SpatialTree* root);

        double FillConstraintViolations(Eigen::VectorXd &phi, int base);
        
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
        int iterNum;
        double ls_step_threshold;
        double backproj_threshold;
        CoordinateLUs coord_lus;
        PolyCurveNetwork* curveNetwork;
        std::vector<Vector3> originalPositions;
        std::vector<double> initialLengths;
        double alpha;
        double beta;
        void SetGradientStep(Eigen::MatrixXd gradient, double delta);
        bool BackprojectConstraints(Eigen::PartialPivLU<Eigen::MatrixXd> &lu);
    };

    template<typename Domain, typename Smoother>
    double TPEFlowSolverSC::ProjectGradientMultigrid(Eigen::MatrixXd &gradients, MultigridHierarchy<Domain>* solver, Eigen::MatrixXd &output) {
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
        Eigen::VectorXd sobolevGradients = solver->template VCycleSolve<Smoother>(gradients3x, 1e-2);
        // Compute dot product with unprojected gradient, and copy into results vector
        double dirDot = gradients3x.dot(sobolevGradients) / (gradients3x.norm() * sobolevGradients.norm());
        output.setZero();
        VectorXdIntoMatrix(sobolevGradients, output);
        return dirDot;
    }

    template<typename Domain, typename Smoother>
    bool TPEFlowSolverSC::BackprojectMultigridLS(Eigen::MatrixXd &gradient, double initGuess,
    MultigridHierarchy<Domain>* solver, SpatialTree* root) {
        double delta = initGuess;
        int nVerts = curveNetwork->NumVertices();
        int nEdges = curveNetwork->NumEdges();
        Eigen::MatrixXd correction(nVerts, 3);

        while (delta > ls_step_threshold) {
            SetGradientStep(gradient, delta);
            if (root) {
                // Update the centers of mass to reflect the new positions
                root->recomputeCentersOfMass(curveNetwork);
            }

            Eigen::VectorXd phi(nEdges + 3);

            double maxBefore = FillConstraintViolations(phi, 0);

            // Compute and apply the correction
            solver->template BackprojectMultigrid<Smoother>(curveNetwork, phi, correction);
            for (int i = 0; i < nVerts; i++) {
                CurveVertex* p = curveNetwork->GetVertex(i);
                Vector3 cur = p->Position();
                p->SetPosition(cur + SelectRow(correction, i));
            }

            // Add length violations to RHS
            double maxViolation = FillConstraintViolations(phi, 0);
            std::cout << "  Constraint: " << maxBefore << " -> " << maxViolation << std::endl;

            if (maxViolation < backproj_threshold) {
                return true;
            }
            else {
                delta /= 2;
            }
        }

        std::cout << "Couldn't make backprojection succeed (initial step " << initGuess << ")" << std::endl;
        return false;
    }
}
