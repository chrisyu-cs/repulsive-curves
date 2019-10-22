#include "tpe_flow_sc.h"
#include "utils.h"
#include "product/dense_matrix.h"

#include "circle_search.h"

namespace LWS {

    TPEFlowSolverSC::TPEFlowSolverSC(PolyCurveNetwork* g) :
    originalPositions(g->NumVertices()),
    initialLengths(g->NumEdges())
    {

        curveNetwork = g;
        alpha = 2;
        beta = 4;
        ls_step_threshold = 1e-15;
        backproj_threshold = 1e-3;
        iterNum = 0;

        useEdgeLengthConstraint = true;

        for (int i = 0; i < g->NumEdges(); i++) {
            CurveEdge* p = g->GetEdge(i);
            initialLengths[i] = p->Length();
        }
    }

    double TPEFlowSolverSC::CurrentEnergy(SpatialTree *root) {
        if (root) return CurrentEnergyBH(root);
        else return CurrentEnergyDirect();
    }

    inline double TPEFlowSolverSC::CurrentEnergyDirect() {
        return TPESC::tpe_total(curveNetwork, alpha, beta);
    }

    double TPEFlowSolverSC::CurrentEnergyBH(SpatialTree *root) {
        return SpatialTree::TPEnergyBH(curveNetwork, root, alpha, beta);
    }

    void TPEFlowSolverSC::FillGradientSingle(Eigen::MatrixXd &gradients, int i, int j) {
        if (i == j) return;
        CurveVertex* i_pt = curveNetwork->GetVertex(i);
        CurveVertex* j_pt = curveNetwork->GetVertex(j);

        // Add i and neighbors of i
        std::vector<CurveVertex*> i_pts;
        i_pts.push_back(i_pt);
        for (int e = 0; e < i_pt->numEdges(); e++) {
            i_pts.push_back(i_pt->edge(e)->Opposite(i_pt));
        }

        // Add j and neighbors of j
        std::vector<CurveVertex*> j_pts;
        j_pts.push_back(j_pt);
        for (int e = 0; e < j_pt->numEdges(); e++) {
            j_pts.push_back(j_pt->edge(e)->Opposite(j_pt));
        }

        // Differentiate wrt neighbors of i
        for (CurveVertex* i_n : i_pts) {
            AddToRow(gradients, i_n->GlobalIndex(), TPESC::tpe_grad(i_pt, j_pt, alpha, beta, i_n));
        }
        // Differentiate wrt neighbors of j
        for (CurveVertex* j_n : j_pts) {
            bool noOverlap = true;
            // Only compute this derivative if j_n is not already included in one of the previous pairs
            for (CurveVertex* i_n : i_pts) {
                if (i_n == j_n) noOverlap = false;
            }
            if (noOverlap) {
                AddToRow(gradients, j_n->GlobalIndex(), TPESC::tpe_grad(i_pt, j_pt, alpha, beta, j_n));
            }
        }
    }

    void TPEFlowSolverSC::FillGradientVectorDirect(Eigen::MatrixXd &gradients) {
        int nVerts = curveNetwork->NumVertices();
        // Fill with zeros, so that the constraint entries are 0
        gradients.setZero();
        // Fill vertex entries with accumulated gradients
        for (int i = 0; i < nVerts; i++) {
            for (int j = 0; j < nVerts; j++) {
                if (i == j) continue;
                FillGradientSingle(gradients, i, j);
            }
        }
    }

    void TPEFlowSolverSC::FillGradientVectorBH(SpatialTree *root, Eigen::MatrixXd &gradients) {
        // Use the spatial tree and Barnes-Hut to evaluate the gradient
        SpatialTree::TPEGradientBarnesHut(curveNetwork, root, gradients, alpha, beta);
    }

    void TPEFlowSolverSC::FillConstraintVector(Eigen::MatrixXd &gradients) {
        int nVerts = curveNetwork->NumVertices();

        for (int i = 0; i < gradients.rows(); i++) {
            SetRow(gradients, i, Vector3{0, 0, 0});
        }

        for (int i = 0; i < nVerts; i++) {
            CurveVertex* i_pt = curveNetwork->GetVertex(i);
            Vector3 len_grad = i_pt->TotalLengthGradient();
            SetRow(gradients, i_pt->GlobalIndex(), len_grad);
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

    double TPEFlowSolverSC::LineSearchStep(Eigen::MatrixXd &gradient, double gradDot, SpatialTree* root) {
        double gradNorm = gradient.norm();
        //std::cout << "Norm of gradient = " << gradNorm << std::endl;
        double initGuess = (gradNorm > 1) ? 1.0 / gradNorm : 1.0 / sqrt(gradNorm);
        return LineSearchStep(gradient, initGuess, 1, gradDot, root);
    }

    double TPEFlowSolverSC::LineSearchStep(Eigen::MatrixXd &gradient, double initGuess, int doublingLimit,
    double gradDot, SpatialTree* root) {
        double delta = initGuess;

        // Save initial positions
        SaveCurrentPositions();

        double initialEnergy = CurrentEnergy(root);
        double gradNorm = gradient.norm();
        int numBacktracks = 0, numDoubles = 0;
        double sigma = 0.01f;
        double newEnergy = initialEnergy;

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
            std::cout << "  Energy: " << initialEnergy << " -> " << newEnergy << std::endl;
            return delta;
        }
    }

    double TPEFlowSolverSC::LSBackproject(Eigen::MatrixXd &gradient, double initGuess,
    Eigen::PartialPivLU<Eigen::MatrixXd> &lu, double gradDot, SpatialTree* root) {
        double delta = initGuess;

        while (delta > ls_step_threshold) {
            SetGradientStep(gradient, delta);
            if (root) {
                // Update the centers of mass to reflect the new positions
                root->recomputeCentersOfMass(curveNetwork);
            }

            bool backprojSuccess = BackprojectConstraints(lu);

            if (backprojSuccess) {
                return delta;
            }
            else {
                delta /= 2;
            }
        }
        std::cout << "Couldn't make backprojection succeed (initial step " << initGuess << ")" << std::endl;
        BackprojectConstraints(lu);
        return 0;
    }

    bool TPEFlowSolverSC::StepLS() {
        int nVerts = curveNetwork->NumVertices();
        Eigen::MatrixXd gradients(nVerts, 3);
        gradients.setZero();
        FillGradientVectorDirect(gradients);
        return LineSearchStep(gradients);
    }

    double TPEFlowSolverSC::ProjectSoboSloboGradient(Eigen::PartialPivLU<Eigen::MatrixXd> &lu, Eigen::MatrixXd &gradients) {
        int nVerts = curveNetwork->NumVertices();
        Eigen::VectorXd b;
        b.setZero(matrixNumRows());

        if (useEdgeLengthConstraint) {
            // If using per-edge length constraints, then the matrix has all coordinates merged,
            // so we only need one solve

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

        else {
            // If not using per-edge length constraints, we solve each coordinate separately
            // Solve for x
            b = gradients.col(0);
            Eigen::VectorXd ss_grad_x = lu.solve(b);

            // Solve for y
            b = gradients.col(1);
            Eigen::VectorXd ss_grad_y = lu.solve(b);

            // Solve for z
            b = gradients.col(2);
            Eigen::VectorXd ss_grad_z = lu.solve(b);

            // Copy the projected gradients back into the vector
            for (int i = 0; i < nVerts; i++) {
                SetRow(gradients, i, Vector3{ss_grad_x(i), ss_grad_y(i), ss_grad_z(i)});
            }

            return 1;
        }
    }

    template<typename T>
    void TestMultiply(T &mult, Eigen::VectorXd &xVec, Eigen::VectorXd &result) {
        mult->Multiply(xVec, result);
    }

    void TPEFlowSolverSC::ExpandMatrix3x(Eigen::MatrixXd &A, Eigen::MatrixXd &B) {
        int nRows = A.rows();
        int nCols = A.cols();

        for (int i = 0; i < nRows; i++) {
            for (int j = 0; j < nCols; j++) {
                // Duplicate this value 3 times along the diagonal of a corresponding 3x3 block
                B(3 * i,     3 * j    ) = A(i, j);
                B(3 * i + 1, 3 * j + 1) = A(i, j);
                B(3 * i + 2, 3 * j + 2) = A(i, j);
            }
        }
    }

    bool TPEFlowSolverSC::BackprojectConstraints(Eigen::PartialPivLU<Eigen::MatrixXd> &lu) {
        int nVerts = curveNetwork->NumVertices();
        Eigen::VectorXd b;
        b.setZero(matrixNumRows());

        Vector3 barycenter = curveNetwork->Barycenter();

        if (useEdgeLengthConstraint) {
            // If using per-edge length constraints, matrix has all coordinates merged,
            // so we only need one solve.

            // Place all 3 barycenter coordinates on RHS
            double maxViolation = curveNetwork->FillConstraintViolations(b, initialLengths);
            
            // Solve for correction
            Eigen::VectorXd corr = lu.solve(b);
            // Apply correction
            for (int i = 0; i < nVerts; i++) {
                Vector3 correction{corr(3 * i), corr(3 * i + 1), corr(3 * i + 2)};
                CurveVertex* pt = curveNetwork->GetVertex(i);
                pt->SetPosition(pt->Position() + correction);
            }

            return maxViolation < backproj_threshold;
        }

        else {
            // If not using per-edge length constraints, matrix is per-coordinate and not merged,
            // so we need to solve separately for each coordinate.

            // Solve for x correction
            b(nVerts) = -barycenter.x;
            Eigen::VectorXd corr_x = lu.solve(b);
            // Solve for y correction
            b(nVerts) = -barycenter.y;
            Eigen::VectorXd corr_y = lu.solve(b);
            // Solve for z correction
            b(nVerts) = -barycenter.z;
            Eigen::VectorXd corr_z = lu.solve(b);

            // Apply the correction
            for (int i = 0; i < nVerts; i++) {
                Vector3 correction{corr_x(i), corr_y(i), corr_z(i)};
                CurveVertex* pt = curveNetwork->GetVertex(i);
                pt->SetPosition(pt->Position() + correction);
            }

            return true;
        }
    }

    int TPEFlowSolverSC::matrixNumRows() {
        if (useEdgeLengthConstraint) {
            int nVerts = curveNetwork->NumVertices();
            // 3V rows for all vertex values; 3 rows for barycenter constraint; V rows for edge length constraints
            return 3 * nVerts + 3 + curveNetwork->NumEdges();
        }
        else {
            // V rows for vertex values by coordinate; 1 row for barycenter constraints by coordinate
            return curveNetwork->NumVertices() + 1;
        }
    }

    double TPEFlowSolverSC::ComputeAndProjectGradient(Eigen::MatrixXd &gradients) {
        Eigen::MatrixXd A;
        Eigen::PartialPivLU<Eigen::MatrixXd> lu;
        return ComputeAndProjectGradient(gradients, A, lu);
    }

    double TPEFlowSolverSC::ComputeAndProjectGradient(Eigen::MatrixXd &gradients, Eigen::MatrixXd &A, Eigen::PartialPivLU<Eigen::MatrixXd> &lu) {
        int nRows = matrixNumRows();
        size_t nVerts = curveNetwork->NumVertices();

        Eigen::MatrixXd l2gradients = gradients;
        Eigen::MatrixXd vertConstraints;
        vertConstraints.setZero(curveNetwork->NumEdges() + 1, 3);

        // If we're not using the per-edge length constraints, use total length instead
        if (!useEdgeLengthConstraint) {
            double ss_start = Utils::currentTimeMilliseconds();
            // Fill the Sobo-Slobo matrix, one entry per vertex
            A.setZero(nRows, nRows);
            SobolevCurves::SobolevPlusBarycenter(curveNetwork, alpha, beta, A);
            double ss_end = Utils::currentTimeMilliseconds();

            std::cout << "  Assemble saddle matrix: " << (ss_end - ss_start) << " ms" << std::endl;

            double factor_start = Utils::currentTimeMilliseconds();
            // Factorize it
            lu.compute(A);

            // Solve for the Sobolev gradient
            ProjectSoboSloboGradient(lu, gradients);
            double factor_end = Utils::currentTimeMilliseconds();

            std::cout << "  Solve for descent direction: " << (factor_end - factor_start) << " ms" << std::endl;

            double constr_start = Utils::currentTimeMilliseconds();
            // Compute gradient of total length constraint
            FillConstraintVector(vertConstraints);
            // Solve for Sobolev gradient of total length
            ProjectSoboSloboGradient(lu, vertConstraints);
            // Project out the length gradient direction
            SobolevCurves::SobolevOrthoProjection(gradients, vertConstraints, A);
            double constr_end = Utils::currentTimeMilliseconds();

            std::cout << "  Project total length constraint: " << (constr_end - constr_start) << " ms" << std::endl;
        }

        else {
            double ss_start = Utils::currentTimeMilliseconds();
            Eigen::MatrixXd A_temp;
            // A is the saddle matrix with the correct dimensions (one row for each coordinate of each gradient entry, plus constraints)
            A.setZero(nRows, nRows);
            // A_temp is a smaller saddle matrix (one row for each gradient entry, only gradient + barycenter, with coordinates lumped together)
            A_temp.setZero(nVerts + 1, nVerts + 1);
            SobolevCurves::SobolevPlusBarycenter(curveNetwork, alpha, beta, A_temp);

            // Duplicate the saddle matrix portion into the larger matrix
            ExpandMatrix3x(A_temp, A);
            // Add rows for edge length constraint
            SobolevCurves::AddEdgeLengthConstraints(curveNetwork, A, 3 * nVerts + 3);
            double ss_end = Utils::currentTimeMilliseconds();

            std::cout << "  Assemble saddle matrix: " << (ss_end - ss_start) << " ms" << std::endl;

            double factor_start = Utils::currentTimeMilliseconds();
            // Factorize it
            lu.compute(A);
            // Solve for the Sobolev gradient
            ProjectSoboSloboGradient(lu, gradients);
            double factor_end = Utils::currentTimeMilliseconds();

            std::cout << "  Solve for descent direction: " << (factor_end - factor_start) << " ms" << std::endl;
        }

        double dot_acc = 0;
        double norm_l2 = l2gradients.norm();
        double norm_w2 = gradients.norm();

        for (int i = 0; i < gradients.rows(); i++) {
            dot_acc += dot(SelectRow(l2gradients, i), SelectRow(gradients, i));
        }

        double dir_dot = dot_acc / (sqrt(norm_l2) * sqrt(norm_w2));
        std::cout << "  Directional dot product = " << dir_dot << std::endl;
        std::cout << "  Sobolev gradient norm = " << dot_acc << std::endl;

        return dir_dot;
    }

    bool TPEFlowSolverSC::StepSobolevLS(bool useBH) {
        long start = Utils::currentTimeMilliseconds();

        size_t nVerts = curveNetwork->NumVertices();
        int nRows = matrixNumRows();

        Eigen::MatrixXd vertGradients;
        vertGradients.setZero(nVerts + 1, 3);

        long grad_start = Utils::currentTimeMilliseconds();
        SpatialTree *tree_root = 0;

        // Barnes-Hut for gradient accumulation
        if (useBH) {
            tree_root = CreateBVHFromCurve(curveNetwork);
            FillGradientVectorBH(tree_root, vertGradients);
        }
        else {
            FillGradientVectorDirect(vertGradients);
        }
        
        std::cout << "\n====== Timing ======" << std::endl;
        double grad_end = Utils::currentTimeMilliseconds();

        std::cout << "  Assemble gradient: " << (grad_end - grad_start) << " ms" << std::endl;

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

    double TPEFlowSolverSC::FillConstraintViolations(Eigen::VectorXd &phi) {
        return curveNetwork->FillConstraintViolations(phi, initialLengths);
    }

    bool TPEFlowSolverSC::StepSobolevLSIterative(double epsilon) {
        std::cout << "=== Iteration " << ++iterNum << " ===" << std::endl;
        long all_start = Utils::currentTimeMilliseconds();

        size_t nVerts = curveNetwork->NumVertices();
        int nRows = matrixNumRows();
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

        // Set up multigrid stuff
        long mg_setup_start = Utils::currentTimeMilliseconds();
        using MultigridDomain = EdgeLengthNullProjectorDomain;
        using MultigridSolver = MultigridHierarchy<MultigridDomain>;
        MultigridDomain* domain = new MultigridDomain(curveNetwork, alpha, beta, 0.5, epsilon);
        MultigridSolver* multigrid = new MultigridSolver(domain);
        long mg_setup_end = Utils::currentTimeMilliseconds();
        std::cout << "  Multigrid setup: " << (mg_setup_end - mg_setup_start) << " ms" << std::endl;

        // Use multigrid to compute the Sobolev gradient
        long mg_start = Utils::currentTimeMilliseconds();
        double dot_acc = ProjectGradientMultigrid<MultigridDomain, MultigridSolver::EigenCG>(vertGradients, multigrid, vertGradients);
        long mg_end = Utils::currentTimeMilliseconds();
        std::cout << "  Multigrid solve: " << (mg_end - mg_start) << " ms" << std::endl;

        // Take a line search step using this gradient
        long ls_start = Utils::currentTimeMilliseconds();
        double step_size = CircleSearch::CircleSearchStep<MultigridSolver, MultigridSolver::EigenCG>(curveNetwork,
            vertGradients, l2gradients, tree_root, multigrid, initialLengths, dot_acc, alpha, beta, 1e-6);
        // double step_size = LineSearchStep(vertGradients, dot_acc, tree_root);
        long ls_end = Utils::currentTimeMilliseconds();
        std::cout << "  Line search: " << (ls_end - ls_start) << " ms" << std::endl;

        // Correct for drift with backprojection
        long bp_start = Utils::currentTimeMilliseconds();
        // step_size = BackprojectMultigridLS<MultigridDomain, MultigridSolver::EigenCG>(vertGradients, step_size, multigrid, tree_root);
        long bp_end = Utils::currentTimeMilliseconds();
        std::cout << "  Backprojection: " << (bp_end - bp_start) << " ms" << std::endl;

        delete multigrid;
        if (tree_root) delete tree_root;

        long all_end = Utils::currentTimeMilliseconds();
        std::cout << "  Total time: " << (all_end - all_start) << " ms" << std::endl;

        return step_size > 0;
    }
}
