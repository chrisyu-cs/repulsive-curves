#include "tpe_flow_sc.h"
#include "utils.h"
#include "product/dense_matrix.h"

namespace LWS {

    TPEFlowSolverSC::TPEFlowSolverSC(PolyCurveNetwork* g, double a, double b) :
    constraint(g),
    originalPositions(g->NumVertices())
    {
        curveNetwork = g;
        alpha = a;
        beta = b;
        ls_step_threshold = 1e-15;
        backproj_threshold = 1e-4;
        iterNum = 0;

        UpdateTargetLengths();
    }

    TPEFlowSolverSC::~TPEFlowSolverSC() {
        for (size_t i = 0; i < obstacles.size(); i++) {
            delete obstacles[i];
        }
        obstacles.clear();
    }

    void TPEFlowSolverSC::UpdateTargetLengths() {
        constraint.UpdateTargetValues(constraintTargets);
    }

    double TPEFlowSolverSC::CurrentEnergy(SpatialTree *root) {
        if (root) return CurrentEnergyBH(root);
        else return CurrentEnergyDirect();
    }

    double TPEFlowSolverSC::CurrentEnergyDirect() {
        return TPESC::tpe_total(curveNetwork, alpha, beta);
    }

    double TPEFlowSolverSC::CurrentEnergyBH(SpatialTree *root) {
        return SpatialTree::TPEnergyBH(curveNetwork, root, alpha, beta);
    }

    void TPEFlowSolverSC::FillGradientVectorDirect(Eigen::MatrixXd &gradients) {
        TPESC::FillGradientVectorDirect(curveNetwork, gradients, alpha, beta);
    }

    void TPEFlowSolverSC::FillGradientVectorBH(SpatialTree *root, Eigen::MatrixXd &gradients) {
        // Use the spatial tree and Barnes-Hut to evaluate the gradient
        SpatialTree::TPEGradientBarnesHut(curveNetwork, root, gradients, alpha, beta);
    }

    void TPEFlowSolverSC::AddObstacleGradients(Eigen::MatrixXd &gradients) {
        for (Obstacle* obs : obstacles) {
            obs->AddGradient(curveNetwork, gradients, alpha, beta);
        }
    }

    bool TPEFlowSolverSC::StepNaive(double h) {
        // Takes a fixed time step h using the L2 gradient
        int nVerts = curveNetwork->NumVertices();
        Eigen::MatrixXd gradients(nVerts, 3);
        gradients.setZero();
        FillGradientVectorDirect(gradients);

        for (int i = 0; i < nVerts; i++) {
            CurveVertex* pt = curveNetwork->GetVertex(i);
            pt->SetPosition(pt->Position() - h * SelectRow(gradients, i));
        }
        return true;
    }

    void TPEFlowSolverSC::SaveCurrentPositions() {
        for (int i = 0; i < curveNetwork->NumVertices(); i++) {
            originalPositions[i] = curveNetwork->GetVertex(i)->Position();
        }
    }

    void TPEFlowSolverSC::RestoreOriginalPositions() {
        for (int i = 0; i < curveNetwork->NumVertices(); i++) {
            curveNetwork->GetVertex(i)->SetPosition(originalPositions[i]);
        }
    }

    void TPEFlowSolverSC::SetGradientStep(Eigen::MatrixXd gradient, double delta) {
        // Write the new vertex positions to the mesh
        // Step every vertex by the gradient times delta
        for (int i = 0; i < curveNetwork->NumVertices(); i++) {
            Vector3 vertGrad = SelectRow(gradient, i);
            curveNetwork->GetVertex(i)->SetPosition(originalPositions[i] - delta * vertGrad);
        }
    }

    double TPEFlowSolverSC::LineSearchStep(Eigen::MatrixXd &gradient, double gradDot, BVHNode3D* root) {
        double gradNorm = gradient.norm();
        //std::cout << "Norm of gradient = " << gradNorm << std::endl;
        double initGuess = (gradNorm > 1) ? 1.0 / gradNorm : 1.0 / sqrt(gradNorm);
        return LineSearchStep(gradient, initGuess, 1, gradDot, root);
    }

    double TPEFlowSolverSC::LineSearchStep(Eigen::MatrixXd &gradient, double initGuess, int doublingLimit,
    double gradDot, BVHNode3D* root) {
        double delta = initGuess;

        // Save initial positions
        SaveCurrentPositions();

        double initialEnergy = CurrentEnergy(root);
        double gradNorm = gradient.norm();
        int numBacktracks = 0, numDoubles = 0;
        double sigma = 0.01f;
        double newEnergy = initialEnergy;

        if (gradNorm < 1e-10) {
            std::cout << "  Gradient is very close to zero" << std::endl;
            return 0;
        }

        // std::cout << "Initial energy " << initialEnergy << std::endl;

        while (delta > ls_step_threshold) {
            SetGradientStep(gradient, delta);
            if (root) {
                // Update the centers of mass to reflect the new positions
                root->recomputeCentersOfMass(curveNetwork);
            }
            newEnergy = CurrentEnergy(root);

            double decrease = initialEnergy - newEnergy;
            double targetDecrease = sigma * delta * gradNorm * gradDot;

            // If the energy hasn't decreased enough to meet the Armijo condition,
            // halve the step size.
            if (decrease < targetDecrease) {
                delta /= 2;
                numBacktracks++;
            }
            else if (decrease >= targetDecrease && numBacktracks == 0 && numDoubles < doublingLimit) {
                delta *= 2;
                numDoubles++;
            }
            // Otherwise, accept the current step.
            else {
                // Empirically we observe that we often need to cut the step size by at least this much
                // before backprojection succeeds, so we just do it now.
                delta *= 0.25;
                SetGradientStep(gradient, delta);
                if (root) {
                    // Update the centers of mass to reflect the new positions
                    root->recomputeCentersOfMass(curveNetwork);
                }
                break;
            }
        }

        if (delta <= ls_step_threshold) {
            std::cout << "Failed to find a non-trivial step after " << numBacktracks << " backtracks" << std::endl;

            // PlotEnergyInDirection(gradient, sigma * gradDot);
            // Restore initial positions if step size goes to 0
            RestoreOriginalPositions();
            return 0;
        }
        else {
            std::cout << "  Energy: " << initialEnergy << " -> " << newEnergy
                << " (step size " << delta << ", " << numBacktracks << " backtracks)" << std::endl;
            return delta;
        }
    }

    double TPEFlowSolverSC::LSBackproject(Eigen::MatrixXd &gradient, double initGuess,
    Eigen::PartialPivLU<Eigen::MatrixXd> &lu, double gradDot, BVHNode3D* root) {
        double delta = initGuess;
        double lastValue = -1;
        int attempts = 1;

        while (delta > ls_step_threshold) {
            SetGradientStep(gradient, delta);
            if (root) {
                // Update the centers of mass to reflect the new positions
                root->recomputeCentersOfMass(curveNetwork);
            }

            double maxValue = BackprojectConstraints(lu);
            
            if (maxValue < backproj_threshold) {
                std::cout << "Backprojection successful after " << attempts << " attempts" << std::endl;
                return delta;
            }
            else {
                // std::cout << "Backproj unsuccessful; halving step size" << std::endl;
                lastValue = maxValue;
                delta /= 2;
                attempts++;
            }
        }
        std::cout << "Couldn't make backprojection succeed after " << attempts << " attempts (initial step " << initGuess << ")" << std::endl;
        BackprojectConstraints(lu);
        return delta;
    }

    bool TPEFlowSolverSC::StepLS() {
        int nVerts = curveNetwork->NumVertices();
        Eigen::MatrixXd gradients(nVerts, 3);
        gradients.setZero();

        FillGradientVectorDirect(gradients);
        bool good_step = LineSearchStep(gradients);

        return good_step;
    }

    bool TPEFlowSolverSC::StepLSConstrained() {
        std::cout << "=== Iteration " << ++iterNum << " ===" << std::endl;
        // Compute gradient
        int nVerts = curveNetwork->NumVertices();
        Eigen::MatrixXd gradients(nVerts, 3);
        gradients.setZero();
        FillGradientVectorDirect(gradients);

        // Set up saddle matrix
        Eigen::MatrixXd A;
        int nRows = constraint.NumConstraintRows() + constraint.NumExpectedCols();
        A.setZero(nRows, nRows);

        // Assemble constraint saddle matrix with identity in the upper-left
        Eigen::MatrixXd mass;
        mass.setIdentity(nVerts * 3, nVerts * 3);
        for (int i = 0; i < nVerts; i++) {
            double m = 1.0 / curveNetwork->GetVertex(i)->DualLength();
            mass(3 * i, 3 * i) = m;
            mass(3 * i + 1, 3 * i + 1) = m;
            mass(3 * i + 2, 3 * i + 2) = m;
        }

        A.block(0, 0, nVerts * 3, nVerts * 3) = mass;
        constraint.FillDenseBlock(A);
        // Factorize it
        Eigen::PartialPivLU<Eigen::MatrixXd> lu;
        lu.compute(A);

        // Project gradient onto constraint differential
        ProjectSoboSloboGradient(lu, gradients);
        double step_size = LineSearchStep(gradients);

        // Backprojection
        LSBackproject(gradients, step_size, lu, 1, 0);

        return (step_size > 0);
    }

    double TPEFlowSolverSC::ProjectSoboSloboGradient(Eigen::PartialPivLU<Eigen::MatrixXd> &lu, Eigen::MatrixXd &gradients) {
        // If using per-edge length constraints, then the matrix has all coordinates merged,
        // so we only need one solve
        int nVerts = curveNetwork->NumVertices();
        Eigen::VectorXd b;
        b.setZero(constraint.NumConstraintRows() + constraint.NumExpectedCols());

        // Fill in RHS with all coordinates
        for (int i = 0; i < nVerts; i++) {
            b(3 * i) = gradients(i, 0);
            b(3 * i + 1) = gradients(i, 1);
            b(3 * i + 2) = gradients(i, 2);
        }
        // Solve for all coordinates
        Eigen::VectorXd ss_grad = lu.solve(b);

        for (int i = 0; i < nVerts; i++) {
            SetRow(gradients, i, Vector3{ss_grad(3 * i), ss_grad(3 * i + 1), ss_grad(3 * i + 2)});
        }

        return 1;
    }

    double TPEFlowSolverSC::BackprojectConstraints(Eigen::PartialPivLU<Eigen::MatrixXd> &lu) {
        int nVerts = curveNetwork->NumVertices();
        Eigen::VectorXd b;
        b.setZero(constraint.NumExpectedCols() + constraint.NumConstraintRows());

        // If using per-edge length constraints, matrix has all coordinates merged,
        // so we only need one solve.

        // Fill RHS with negative constraint values
        double maxViolation = constraint.FillConstraintValues(b, constraintTargets, 3 * nVerts);

        // Solve for correction
        Eigen::VectorXd corr = lu.solve(b);
        // Apply correction
        for (int i = 0; i < nVerts; i++) {
            Vector3 correction{corr(3 * i), corr(3 * i + 1), corr(3 * i + 2)};
            CurveVertex* pt = curveNetwork->GetVertex(i);
            pt->SetPosition(pt->Position() + correction);
        }

        std::cout << "  Constraint value = " << maxViolation << std::endl;

        return maxViolation;
    }

    double TPEFlowSolverSC::ComputeAndProjectGradient(Eigen::MatrixXd &gradients, Eigen::MatrixXd &A, Eigen::PartialPivLU<Eigen::MatrixXd> &lu) {
        size_t nVerts = curveNetwork->NumVertices();
        Eigen::MatrixXd l2gradients = gradients;

        // Assemble the Sobolev gram matrix with constraints
        double ss_start = Utils::currentTimeMilliseconds();
        SobolevCurves::Sobolev3XWithConstraints(curveNetwork, constraint, alpha, beta, A);
        double ss_end = Utils::currentTimeMilliseconds();
        std::cout << "  Assemble saddle matrix: " << (ss_end - ss_start) << " ms" << std::endl;

        // Factorize and solve
        double factor_start = Utils::currentTimeMilliseconds();
        lu.compute(A);
        ProjectSoboSloboGradient(lu, gradients);
        double factor_end = Utils::currentTimeMilliseconds();

        std::cout << "  Solve for descent direction: " << (factor_end - factor_start) << " ms" << std::endl;

        double dot_acc = 0;
        double norm_l2 = l2gradients.norm();
        double norm_w2 = gradients.norm();

        for (int i = 0; i < gradients.rows(); i++) {
            dot_acc += dot(SelectRow(l2gradients, i), SelectRow(gradients, i));
        }

        double dir_dot = dot_acc / (norm_l2 * norm_w2);
        std::cout << "  Directional dot product = " << dir_dot << std::endl;
        std::cout << "  Sobolev gradient norm = " << dot_acc << std::endl;

        return dir_dot;
    }

    bool TPEFlowSolverSC::StepSobolevLS(bool useBH) {
        long start = Utils::currentTimeMilliseconds();

        size_t nVerts = curveNetwork->NumVertices();

        Eigen::MatrixXd vertGradients;
        vertGradients.setZero(nVerts + 1, 3);

        long grad_start = Utils::currentTimeMilliseconds();
        BVHNode3D *tree_root = 0;

        // Barnes-Hut for gradient accumulation
        if (useBH) {
            tree_root = CreateBVHFromCurve(curveNetwork);
            FillGradientVectorBH(tree_root, vertGradients);
        }
        else {
            FillGradientVectorDirect(vertGradients);
        }

        // Add gradient contributions from obstacles
        AddObstacleGradients(vertGradients);
        
        std::cout << "=== Iteration " << ++iterNum << " ===" << std::endl;
        double grad_end = Utils::currentTimeMilliseconds();

        std::cout << "  Assemble gradient " << (useBH ? "(Barnes-Hut)" : "(direct)") << ": " << (grad_end - grad_start) << " ms" << std::endl;

        double length1 = curveNetwork->TotalLength();

        Eigen::MatrixXd A;
        Eigen::PartialPivLU<Eigen::MatrixXd> lu;

        // Compute the Sobolev gradient
        double dot_acc = ComputeAndProjectGradient(vertGradients, A, lu);

        // Take a line search step using this gradient
        double ls_start = Utils::currentTimeMilliseconds();
        double step_size = LineSearchStep(vertGradients, dot_acc, tree_root);
        double ls_end = Utils::currentTimeMilliseconds();
        std::cout << "  Line search: " << (ls_end - ls_start) << " ms" << std::endl;

        if (step_size <= 0) {
            std::cout << "  No line search step taken; stopping." << std::endl;
            return 0;
        }

        // Correct for drift with backprojection
        double bp_start = Utils::currentTimeMilliseconds();
        step_size = LSBackproject(vertGradients, step_size, lu, dot_acc, tree_root);
        double bp_end = Utils::currentTimeMilliseconds();
        std::cout << "  Backprojection: " << (bp_end - bp_start) << " ms" << std::endl;

        if (tree_root) {
            delete tree_root;
        }

        double length2 = curveNetwork->TotalLength();
        std::cout << "Length " << length1 << " -> " << length2 << std::endl;
        long end = Utils::currentTimeMilliseconds();
        std::cout << "Time = " << (end - start) << " ms" << std::endl;

        return step_size > 0;
    }

    bool TPEFlowSolverSC::StepSobolevLSIterative(double epsilon) {
        std::cout << "=== Iteration " << ++iterNum << " ===" << std::endl;
        long all_start = Utils::currentTimeMilliseconds();

        size_t nVerts = curveNetwork->NumVertices();
        BVHNode3D* tree_root = 0;

        Eigen::MatrixXd vertGradients;
        vertGradients.setZero(nVerts, 3);

        // Assemble the L2 gradient
        long bh_start = Utils::currentTimeMilliseconds();
        tree_root = CreateBVHFromCurve(curveNetwork);
        FillGradientVectorBH(tree_root, vertGradients);
        Eigen::MatrixXd l2gradients = vertGradients;
        long bh_end = Utils::currentTimeMilliseconds();
        std::cout << "  Barnes-Hut: " << (bh_end - bh_start) << " ms" << std::endl;

        // Add gradient contributions from obstacles
        AddObstacleGradients(vertGradients);
        std::cout << "Added obstacle gradients" << std::endl;

        // Set up multigrid stuff
        long mg_setup_start = Utils::currentTimeMilliseconds();
        using MultigridDomain = ConstraintProjectorDomain<ConstraintType>;
        using MultigridSolver = MultigridHierarchy<MultigridDomain>;
        MultigridDomain* domain = new MultigridDomain(curveNetwork, alpha, beta, 0.5, epsilon);
        std::cout << "Made domain" << std::endl;
        MultigridSolver* multigrid = new MultigridSolver(domain);
        std::cout << "Made solver" << std::endl;
        long mg_setup_end = Utils::currentTimeMilliseconds();
        std::cout << "  Multigrid setup: " << (mg_setup_end - mg_setup_start) << " ms" << std::endl;

        // Use multigrid to compute the Sobolev gradient
        long mg_start = Utils::currentTimeMilliseconds();
        double dot_acc = ProjectGradientMultigrid<MultigridDomain, MultigridSolver::EigenCG>(vertGradients, multigrid, vertGradients);
        long mg_end = Utils::currentTimeMilliseconds();
        std::cout << "  Multigrid solve: " << (mg_end - mg_start) << " ms" << std::endl;

        // Take a line search step using this gradient
        long ls_start = Utils::currentTimeMilliseconds();
        // double step_size = CircleSearch::CircleSearchStep<MultigridSolver, MultigridSolver::EigenCG>(curveNetwork,
        //     vertGradients, l2gradients, tree_root, multigrid, initialLengths, dot_acc, alpha, beta, 1e-6);
        double step_size = LineSearchStep(vertGradients, dot_acc, tree_root);
        long ls_end = Utils::currentTimeMilliseconds();
        std::cout << "  Line search: " << (ls_end - ls_start) << " ms" << std::endl;

        // Correct for drift with backprojection
        long bp_start = Utils::currentTimeMilliseconds();
        step_size = BackprojectMultigridLS<MultigridDomain, MultigridSolver::EigenCG>(vertGradients, step_size, multigrid, tree_root);
        long bp_end = Utils::currentTimeMilliseconds();
        std::cout << "  Backprojection: " << (bp_end - bp_start) << " ms" << std::endl;

        delete multigrid;
        if (tree_root) delete tree_root;

        long all_end = Utils::currentTimeMilliseconds();
        std::cout << "  Total time: " << (all_end - all_start) << " ms" << std::endl;

        return step_size > 0;
    }
}
