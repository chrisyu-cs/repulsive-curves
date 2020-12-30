#include "tpe_flow_sc.h"
#include "utils.h"
#include "product/dense_matrix.h"

#include "circle_search.h"

namespace LWS {

    TPEFlowSolverSC::TPEFlowSolverSC(PolyCurveNetwork* g, double a, double b) : constraint(g)
    {
        curveNetwork = g;
        alpha = a;
        beta = b;
        ls_step_threshold = 1e-15;
        backproj_threshold = 1e-4;
        iterNum = 0;
        lastStepSize = 0;

        mg_tolerance = 1e-2;

        double averageLength = g->TotalLength();
        averageLength /= g->NumEdges();
        mg_backproj_threshold = fmin(1e-4, averageLength * mg_tolerance * 1e-2);
        std::cout << "Multigrid backprojection threshold set to " << mg_backproj_threshold << std::endl;

        UpdateTargetLengths();
        useEdgeLengthScale = false;
        useTotalLengthScale = false;
        perfLogEnabled = false;
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

    void TPEFlowSolverSC::ReplaceCurve(PolyCurveNetwork* new_p) {
        curveNetwork = new_p;
        constraint = ConstraintClassType(curveNetwork);
        UpdateTargetLengths();
        if (useEdgeLengthScale) {
            double averageLength = curveNetwork->TotalLength() / curveNetwork->NumEdges();
            lengthScaleStep = averageLength / 100;
        }
    }

    void TPEFlowSolverSC::EnablePerformanceLog(std::string logFile) {
        perfFile.open(logFile);
        perfLogEnabled = true;
        std::cout << "Started logging performance to " << logFile << std::endl;
    }

    void TPEFlowSolverSC::ClosePerformanceLog() {
        perfFile.close();
    }

    void TPEFlowSolverSC::SetTotalLengthScaleTarget(double scale) {
        useTotalLengthScale = true;
        int startIndex = constraint.startIndexOfConstraint(ConstraintType::TotalLength);
        if (startIndex < 0) {
            useTotalLengthScale = false;
            std::cout << "No total length constraint; ignoring total length scale" << std::endl;
            return;
        }

        double len = curveNetwork->TotalLength();
        targetLength = scale * len;
        lengthScaleStep = len / 100;
    }

    void TPEFlowSolverSC::SetEdgeLengthScaleTarget(double scale) {
        useEdgeLengthScale = true;

        int startIndex = constraint.startIndexOfConstraint(ConstraintType::EdgeLengths);
        if (startIndex < 0) {
            useEdgeLengthScale = false;
            std::cout << "No edge length constraint; ignoring edge length scale" << std::endl;
            return;
        }
        targetLength = scale * curveNetwork->TotalLength();
        std::cout << "Setting target length to " << targetLength << " (" << scale << " times original)" << std::endl;
        double averageLength = curveNetwork->TotalLength() / curveNetwork->NumEdges();
        lengthScaleStep = averageLength / 100;
    }

    double TPEFlowSolverSC::CurrentEnergy(SpatialTree *root) {
        double energy = 0;

        if (root) energy = TPEnergyBH(root);
        else energy = TPEnergyDirect();

        for (CurvePotential* p : potentials) {
            energy += p->CurrentValue(curveNetwork);
        }

        for (Obstacle* obs : obstacles) {
            energy += obs->ComputeEnergy(curveNetwork);
        }

        return energy;
    }

    double TPEFlowSolverSC::TPEnergyDirect() {
        return TPESC::tpe_total(curveNetwork, alpha, beta);
    }

    double TPEFlowSolverSC::TPEnergyBH(SpatialTree *root)
    {
        return SpatialTree::TPEnergyBH(curveNetwork, root, alpha, beta);
    }

    void TPEFlowSolverSC::FillGradientVectorDirect(Eigen::MatrixXd &gradients) {
        TPESC::FillGradientVectorDirect(curveNetwork, gradients, alpha, beta);
    }

    void TPEFlowSolverSC::FillGradientVectorBH(SpatialTree *root, Eigen::MatrixXd &gradients) {
        // Use the spatial tree and Barnes-Hut to evaluate the gradient
        SpatialTree::TPEGradientBarnesHut(curveNetwork, root, gradients, alpha, beta);
    }

    void TPEFlowSolverSC::AddAllGradients(SpatialTree *tree_root, Eigen::MatrixXd &vertGradients) {
        // Barnes-Hut for gradient accumulation
        if (tree_root) {
            FillGradientVectorBH(tree_root, vertGradients);
        }
        else {
            FillGradientVectorDirect(vertGradients);
        }

        // Add gradient contributions from obstacles
        for (Obstacle* obs : obstacles) {
            obs->AddGradient(curveNetwork, vertGradients);
        }

        // Add gradient contributions from potentials
        for (CurvePotential* p : potentials) {
            p->AddGradient(curveNetwork, vertGradients);
        }
    }
    

    bool TPEFlowSolverSC::StepNaive(double h) {
        // Takes a fixed time step h using the L2 gradient
        int nVerts = curveNetwork->NumVertices();
        Eigen::MatrixXd gradients(nVerts, 3);
        gradients.setZero();
        // FillGradientVectorDirect(gradients);
        AddAllGradients(0, gradients);

        for (int i = 0; i < nVerts; i++) {
            CurveVertex* pt = curveNetwork->GetVertex(i);
            pt->SetPosition(pt->Position() - h * SelectRow(gradients, i));
        }
        return true;
    }

    void TPEFlowSolverSC::SaveCurrentPositions() {
        originalPositionMatrix = curveNetwork->positions;
    }

    void TPEFlowSolverSC::RestoreOriginalPositions() {
        curveNetwork->positions = originalPositionMatrix;
    }

    void TPEFlowSolverSC::SetGradientStep(Eigen::MatrixXd &gradient, double delta) {
        curveNetwork->positions = originalPositionMatrix - delta * gradient;
    }

    double TPEFlowSolverSC::LineSearchStep(Eigen::MatrixXd &gradient, double gradDot, BVHNode3D* root, bool resetStep) {
        double gradNorm = gradient.norm();
        //std::cout << "Norm of gradient = " << gradNorm << std::endl;
        double initGuess = (gradNorm > 1) ? 1.0 / gradNorm : 1.0 / sqrt(gradNorm);
        // Use the step size from the previous iteration, if it exists
        if (!resetStep && lastStepSize > fmax(ls_step_threshold, 1e-5)) {
            initGuess = fmin(lastStepSize * 1.5, initGuess * 4);
        }
        std::cout << "  Starting line search with initial guess " << initGuess << std::endl;
        return LineSearchStep(gradient, initGuess, 0, gradDot, root);
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

    inline double alpha_of_delta(double delta, double alpha_0, double alpha_1, double R) {
        return (1.0 / R) * (alpha_0 * delta + alpha_1 * delta * delta);
    }

    void TPEFlowSolverSC::SetCircleStep(Eigen::MatrixXd &P_dot, Eigen::MatrixXd &K, double sqrt_G, double R, double alpha_delta) {
        curveNetwork->positions = originalPositionMatrix + R * (-P_dot * (sin(alpha_delta) / sqrt_G) + R * K * (1 - cos(alpha_delta)));
    }

    double TPEFlowSolverSC::CircleSearchStep(Eigen::MatrixXd &P_dot, Eigen::MatrixXd &P_ddot, Eigen::MatrixXd &G, BVHNode3D* root) {
        // Save initial positions
        SaveCurrentPositions();
        int nVerts = curveNetwork->NumVertices();
        int nRows = curveNetwork->NumVertices() * 3;

        // TODO: use only Gram matrix or full constraint matrix?
        Eigen::VectorXd P_dot_vec(nRows);
        MatrixIntoVectorX3(P_dot, P_dot_vec);
        Eigen::VectorXd P_ddot_vec(nRows);
        MatrixIntoVectorX3(P_ddot, P_ddot_vec);

        double G_Pd_Pd = P_dot_vec.dot(G.block(0, 0, nRows, nRows) * P_dot_vec);
        double G_Pd_Pdd = P_dot_vec.dot(G.block(0, 0, nRows, nRows) * P_ddot_vec);

        Eigen::VectorXd K_vec = 1.0 / G_Pd_Pd * (-P_ddot_vec + P_dot_vec * (G_Pd_Pdd / G_Pd_Pd));
        double R = 1.0 / sqrt(K_vec.dot(G.block(0, 0, nRows, nRows) * K_vec));
        double alpha_0 = sqrt(G_Pd_Pd);
        double alpha_1 = 0.5 * G_Pd_Pdd / alpha_0;

        Eigen::MatrixXd K(nVerts, 3);
        K.setZero();
        VectorXdIntoMatrix(K_vec, K);

        double delta = -2 * alpha_0 / alpha_1;
        std::cout << "Delta estimate = " << delta << std::endl;
        delta = fmax(0, delta);

        SetCircleStep(P_dot, K, alpha_0, R, alpha_of_delta(delta, alpha_0, alpha_1, R));
        return delta;
    }

    double TPEFlowSolverSC::LSBackproject(Eigen::MatrixXd &gradient, double initGuess,
    Eigen::PartialPivLU<Eigen::MatrixXd> &lu, double gradDot, BVHNode3D* root) {
        double delta = initGuess;
        int attempts = 0;

        while ((delta > ls_step_threshold || useEdgeLengthScale) && attempts < 10) {
            attempts++;
            SetGradientStep(gradient, delta);
            if (root) {
                // Update the centers of mass to reflect the new positions
                root->recomputeCentersOfMass(curveNetwork);
            }

            for (int i = 0; i < 3; i++) {
                double maxValue = BackprojectConstraints(lu);
                if (maxValue < backproj_threshold) {
                    std::cout << "Backprojection successful after " << attempts << " attempts" << std::endl;
                    std::cout << "Used " << (i + 1) << " Newton steps on successful attempt" << std::endl;
                    return delta;
                }
            }
            
            delta /= 2;
        }
        std::cout << "Couldn't make backprojection succeed after " << attempts << " attempts (initial step " << initGuess << ")" << std::endl;
        SetGradientStep(gradient, 0);
        BackprojectConstraints(lu);
        return delta;
    }

    bool TPEFlowSolverSC::StepLS(bool useBH) {
        int nVerts = curveNetwork->NumVertices();
        Eigen::MatrixXd gradients(nVerts, 3);
        gradients.setZero();

        // FillGradientVectorDirect(gradients);
        BVHNode3D *tree_root = 0;
        if (useBH) tree_root = CreateBVHFromCurve(curveNetwork);
        AddAllGradients(tree_root, gradients);
        double gradNorm = gradients.norm();
        double step_size = LineSearchStep(gradients, 1, tree_root);

        delete tree_root;
        soboNormZero = (gradNorm < 1e-4);
        lastStepSize = step_size;
        return (step_size > ls_step_threshold);
    }

    bool TPEFlowSolverSC::StepLSConstrained(bool useBH, bool useBackproj) {
        std::cout << "=== Iteration " << ++iterNum << " ===" << std::endl;
        // Compute gradient
        int nVerts = curveNetwork->NumVertices();
        Eigen::MatrixXd gradients(nVerts, 3);
        gradients.setZero();
        BVHNode3D *tree_root = 0;
        if (useBH) tree_root = CreateBVHFromCurve(curveNetwork);
        AddAllGradients(tree_root, gradients);

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
        double gradNorm = gradients.norm();
        std::cout << "  Norm gradient = " << gradNorm << std::endl;
        double step_size = LineSearchStep(gradients, 1, tree_root);

        // Backprojection
        if (useBackproj) {
            step_size = LSBackproject(gradients, step_size, lu, 1, tree_root);
        }

        if (tree_root) delete tree_root;

        lastStepSize = step_size;
        soboNormZero = (gradNorm < 1e-4);
        return (step_size > ls_step_threshold);
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
        fullDerivVector = ss_grad;

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
        // Compute constraint violation after correction
        maxViolation = constraint.FillConstraintValues(b, constraintTargets, 3 * nVerts);
        std::cout << "  Constraint value = " << maxViolation << std::endl;

        return maxViolation;
    }

    double TPEFlowSolverSC::ProjectGradient(Eigen::MatrixXd &gradients, Eigen::MatrixXd &A, Eigen::PartialPivLU<Eigen::MatrixXd> &lu) {
        size_t nVerts = curveNetwork->NumVertices();
        Eigen::MatrixXd l2gradients = gradients;

        // Assemble the Sobolev gram matrix with constraints
        double ss_start = Utils::currentTimeMilliseconds();
        SobolevCurves::Sobolev3XWithConstraints(curveNetwork, constraint, alpha, beta, A);
        double ss_end = Utils::currentTimeMilliseconds();

        // Factorize and solve
        double factor_start = Utils::currentTimeMilliseconds();
        lu.compute(A);
        ProjectSoboSloboGradient(lu, gradients);
        double factor_end = Utils::currentTimeMilliseconds();

        double soboDot = 0;

        for (int i = 0; i < gradients.rows(); i++) {
            soboDot += dot(SelectRow(l2gradients, i), SelectRow(gradients, i));
        }

        return soboDot;
    }

    void TPEFlowSolverSC::GetSecondDerivative(SpatialTree* tree_root,
    Eigen::MatrixXd &projected1, double epsilon, Eigen::MatrixXd &secondDeriv) {
        Eigen::MatrixXd origPos = curveNetwork->positions;
        int nVerts = curveNetwork->NumVertices();

        // Evaluate second point for circular line search
        double eps = 1e-5;
        // Evaluate new L2 gradients
        curveNetwork->positions -= eps * projected1;
        Eigen::MatrixXd projectedEps;
        projectedEps.setZero(nVerts, 3);
        AddAllGradients(tree_root, projectedEps);
        // Project a second time
        Eigen::MatrixXd A_eps;
        Eigen::PartialPivLU<Eigen::MatrixXd> lu_eps;
        ProjectGradient(projectedEps, A_eps, lu_eps);

        secondDeriv = (projectedEps - projected1) / eps;
        curveNetwork->positions = origPos;
    }

    bool TPEFlowSolverSC::TargetLengthReached() {
        if (useEdgeLengthScale || useTotalLengthScale) {
            double currentLength = curveNetwork->TotalLength();
            return fabs(currentLength - targetLength) <= lengthScaleStep;
        }
        return true;
    }

    void TPEFlowSolverSC::MoveLengthTowardsTarget() {
        if (!useTotalLengthScale && !useEdgeLengthScale) {
            return;
        }

        if (TargetLengthReached()) {
            std::cout << "Target length reached; turning off length scaling" << std::endl;
            useEdgeLengthScale = false;
            useTotalLengthScale = false;
            return;
        }

        if (useEdgeLengthScale) {
            double currentLength = curveNetwork->TotalLength();

            // If we're below the target length, increase edge lengths
            if (currentLength < targetLength) {
                int cStart = constraint.startIndexOfConstraint(ConstraintType::EdgeLengths);
                int nRows = constraint.rowsOfConstraint(ConstraintType::EdgeLengths);

                double diff = targetLength - currentLength;
                double toAdd = lengthScaleStep;
                if (diff < lengthScaleStep * nRows) {
                    toAdd = diff / nRows;
                }

                for (int i = cStart; i < cStart + nRows; i++) {
                    constraintTargets(i) += toAdd;
                }
            }
            else if (currentLength > targetLength) {
                int cStart = constraint.startIndexOfConstraint(ConstraintType::EdgeLengths);
                int nRows = constraint.rowsOfConstraint(ConstraintType::EdgeLengths);

                double diff = currentLength - targetLength;
                double toSubt = lengthScaleStep;
                if (diff < lengthScaleStep * nRows) {
                    toSubt = diff / nRows;
                }

                for (int i = cStart; i < cStart + nRows; i++) {
                    constraintTargets(i) -= toSubt;
                }
            }
        }

        else if (useTotalLengthScale) {
            double currentLength = curveNetwork->TotalLength();

            if (currentLength < targetLength) {
                int cStart = constraint.startIndexOfConstraint(ConstraintType::TotalLength);
                constraintTargets(cStart) += lengthScaleStep;
            }

            else if (currentLength > targetLength) {
                int cStart = constraint.startIndexOfConstraint(ConstraintType::TotalLength);
                constraintTargets(cStart) -= lengthScaleStep;

            }
        }
    }

    bool TPEFlowSolverSC::StepSobolevLS(bool useBH, bool useBackproj) {
        long start = Utils::currentTimeMilliseconds();

        size_t nVerts = curveNetwork->NumVertices();

        Eigen::MatrixXd vertGradients;
        vertGradients.setZero(nVerts, 3);

        // If applicable, move constraint targets
        MoveLengthTowardsTarget();

        // Assemble gradient, either exactly or with Barnes-Hut
        long bh_start = Utils::currentTimeMilliseconds();
        BVHNode3D *tree_root = 0;
        if (useBH) tree_root = CreateBVHFromCurve(curveNetwork);
        AddAllGradients(tree_root, vertGradients);
        Eigen::MatrixXd l2Gradients = vertGradients;

        std::cout << "=== Iteration " << ++iterNum << " ===" << std::endl;
        double bh_end = Utils::currentTimeMilliseconds();

        std::cout << "  Assemble gradient " << (useBH ? "(Barnes-Hut)" : "(direct)") << ": " << (bh_end - bh_start) << " ms" << std::endl;
        std::cout << "  L2 gradient norm = " << l2Gradients.norm() << std::endl;

        double length1 = curveNetwork->TotalLength();

        Eigen::MatrixXd A;
        Eigen::PartialPivLU<Eigen::MatrixXd> lu;

        long project_start = Utils::currentTimeMilliseconds();
        // Compute the Sobolev gradient
        double soboDot = ProjectGradient(vertGradients, A, lu);
        long project_end = Utils::currentTimeMilliseconds();

        std::cout << "  Sobolev gradient norm = " << soboDot << std::endl;
        if (__isnan(soboDot)) {
            std::cout << "Sobolev projection produced NaN; aborting." << std::endl;
            return false;
        }

        double dot_acc = soboDot / (l2Gradients.norm() * vertGradients.norm());
        std::cout << "  Project gradient: " << (project_end - project_start) << " ms" << std::endl;

        // Take a line search step using this gradient
        double ls_start = Utils::currentTimeMilliseconds();
        double step_size = LineSearchStep(vertGradients, dot_acc, tree_root);
        // double step_size = CircleSearchStep(vertGradients, secondDeriv, A, tree_root);
        double ls_end = Utils::currentTimeMilliseconds();
        std::cout << "  Line search: " << (ls_end - ls_start) << " ms" << std::endl;

        if (useEdgeLengthScale && step_size < ls_step_threshold) {
            vertGradients.setZero();
        }

        // Correct for drift with backprojection
        double bp_start = Utils::currentTimeMilliseconds();
        if (useBackproj) {
            step_size = LSBackproject(vertGradients, step_size, lu, dot_acc, tree_root);
        }
        double bp_end = Utils::currentTimeMilliseconds();
        std::cout << "  Backprojection: " << (bp_end - bp_start) << " ms" << std::endl;
        std::cout << "  Final step size = " << step_size << std::endl;

        if (tree_root) {
            delete tree_root;
        }

        double length2 = curveNetwork->TotalLength();
        std::cout << "Length " << length1 << " -> " << length2 << std::endl;
        long end = Utils::currentTimeMilliseconds();
        std::cout << "Time = " << (end - start) << " ms" << std::endl;


        if (perfLogEnabled) {
            double bh_time = bh_end - bh_start;
            double mg_time = project_end - project_start;
            double ls_time = ls_end - ls_start;
            double bp_time = bp_end - bp_start;
            double all_time = end - start;

            perfFile << iterNum << ", " << bh_time << ", " << mg_time << ", " << ls_time << ", " << bp_time << ", " << all_time << std::endl;
        }

        lastStepSize = step_size;
        soboNormZero = (soboDot < 1e-4);

        return step_size > ls_step_threshold;
    }

    bool TPEFlowSolverSC::StepSobolevLSIterative(double epsilon, bool useBackproj) {
        std::cout << "=== Iteration " << ++iterNum << " ===" << std::endl;
        long all_start = Utils::currentTimeMilliseconds();

        size_t nVerts = curveNetwork->NumVertices();
        BVHNode3D* tree_root = 0;

        Eigen::MatrixXd vertGradients;
        vertGradients.setZero(nVerts, 3);

        // If applicable, move constraint targets
        MoveLengthTowardsTarget();

        // Assemble the L2 gradient
        long bh_start = Utils::currentTimeMilliseconds();
        tree_root = CreateBVHFromCurve(curveNetwork);
        AddAllGradients(tree_root, vertGradients);
        Eigen::MatrixXd l2gradients = vertGradients;
        long bh_end = Utils::currentTimeMilliseconds();
        std::cout << "  Barnes-Hut: " << (bh_end - bh_start) << " ms" << std::endl;

        // Set up multigrid stuff
        long mg_setup_start = Utils::currentTimeMilliseconds();
        using MultigridDomain = ConstraintProjectorDomain<ConstraintClassType>;
        using MultigridSolver = MultigridHierarchy<MultigridDomain>;
        double sep = 1.0;
        MultigridDomain* domain = new MultigridDomain(curveNetwork, alpha, beta, sep, epsilon);
        MultigridSolver* multigrid = new MultigridSolver(domain);
        long mg_setup_end = Utils::currentTimeMilliseconds();
        std::cout << "  Multigrid setup: " << (mg_setup_end - mg_setup_start) << " ms" << std::endl;

        // Use multigrid to compute the Sobolev gradient
        long mg_start = Utils::currentTimeMilliseconds();
        double soboDot = ProjectGradientMultigrid<MultigridDomain, MultigridSolver::EigenCG>(vertGradients, multigrid, vertGradients, mg_tolerance);
        double dot_acc = soboDot / (l2gradients.norm() * vertGradients.norm());
        long mg_end = Utils::currentTimeMilliseconds();
        std::cout << "  Multigrid solve: " << (mg_end - mg_start) << " ms" << std::endl;
        std::cout << "  Sobolev gradient norm = " << soboDot << std::endl;

        // Take a line search step using this gradient
        long ls_start = Utils::currentTimeMilliseconds();
        // double step_size = CircleSearch::CircleSearchStep<MultigridSolver, MultigridSolver::EigenCG>(curveNetwork,
        //     vertGradients, l2gradients, tree_root, multigrid, initialLengths, dot_acc, alpha, beta, 1e-6);
        double step_size = LineSearchStep(vertGradients, dot_acc, tree_root);
        long ls_end = Utils::currentTimeMilliseconds();
        std::cout << "  Line search: " << (ls_end - ls_start) << " ms" << std::endl;

        // Correct for drift with backprojection
        long bp_start = Utils::currentTimeMilliseconds();
        if (useBackproj) {
            step_size = LSBackprojectMultigrid<MultigridDomain, MultigridSolver::EigenCG>(vertGradients,
            step_size, multigrid, tree_root, mg_tolerance);
        }
        long bp_end = Utils::currentTimeMilliseconds();
        std::cout << "  Backprojection: " << (bp_end - bp_start) << " ms" << std::endl;
        std::cout << "  Final step size = " << step_size << std::endl;

        delete multigrid;
        if (tree_root) delete tree_root;

        long all_end = Utils::currentTimeMilliseconds();
        std::cout << "  Total time: " << (all_end - all_start) << " ms" << std::endl;

        if (perfLogEnabled) {
            double bh_time = bh_end - bh_start;
            double mg_time = mg_end - mg_setup_start;
            double ls_time = ls_end - ls_start;
            double bp_time = bp_end - bp_start;
            double all_time = all_end - all_start;

            perfFile << iterNum << ", " << bh_time << ", " << mg_time << ", " << ls_time << ", " << bp_time << ", " << all_time << std::endl;
        }

        soboNormZero = (soboDot < 1e-4);

        lastStepSize = step_size;
        return step_size > ls_step_threshold;
    }
}
