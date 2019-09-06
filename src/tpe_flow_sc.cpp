#include "tpe_flow_sc.h"
#include "utils.h"
#include "product/dense_matrix.h"

namespace LWS {

    TPEFlowSolverSC::TPEFlowSolverSC(PolyCurveGroup* g) :
    originalPositions(g->NumVertices()),
    l2gradients(g->NumVertices() + 1, 3),
    vertGradients(g->NumVertices() + 1, 3),
    vertConstraints(g->NumVertices() + 1, 3),
    initialLengths(g->NumVertices())
    {
        l2gradients.setZero();
        vertGradients.setZero();
        vertConstraints.setZero();

        curves = g;
        alpha = 3;
        beta = 6;
        ls_step_threshold = 1e-15;
        backproj_threshold = 1e-3;

        useEdgeLengthConstraint = false;

        for (int i = 0; i < g->NumVertices(); i++) {
            PointOnCurve p = g->GetCurvePoint(i);
            Vector3 p1 = p.Position();
            Vector3 p2 = p.Next().Position();
            initialLengths[i] = norm(p1 - p2);
        }
    }

    double TPEFlowSolverSC::CurrentEnergy(SpatialTree *root) {
        if (root) return CurrentEnergyBH(root);
        else return CurrentEnergyDirect();
    }

    inline double TPEFlowSolverSC::CurrentEnergyDirect() {
        return TPESC::tpe_total(curves, alpha, beta);
    }

    double TPEFlowSolverSC::CurrentEnergyBH(SpatialTree *root) {
        int nVerts = curves->NumVertices();
        double fullSum = 0;
        // Loop over all vertices and add up energy contributions
        for (int i = 0; i < nVerts; i++) {
            PointOnCurve i_pt = curves->GetCurvePoint(i);
            double vertSum = 0;
            root->accumulateVertexEnergy(vertSum, i_pt, curves, alpha, beta);
            fullSum += vertSum;
        }
        return fullSum;
    }

    void TPEFlowSolverSC::FillGradientSingle(Eigen::MatrixXd &gradients, int i, int j) {
        if (i == j) return;
        PointOnCurve i_pt = curves->GetCurvePoint(i);
        PointOnCurve j_pt = curves->GetCurvePoint(j);

        PointOnCurve i_prev = i_pt.Prev();
        PointOnCurve i_next = i_pt.Next();
        
        PointOnCurve j_prev = j_pt.Prev();
        PointOnCurve j_next = j_pt.Next();

        AddToRow(gradients, curves->GlobalIndex(i_prev), TPESC::tpe_grad(i_pt, j_pt, alpha, beta, i_prev));
        AddToRow(gradients, curves->GlobalIndex(i_pt), TPESC::tpe_grad(i_pt, j_pt, alpha, beta, i_pt));
        AddToRow(gradients, curves->GlobalIndex(i_next), TPESC::tpe_grad(i_pt, j_pt, alpha, beta, i_next));

        // If j_prev is not already included in one of the above 3
        if (j_prev != i_prev && j_prev != i_pt && j_prev != i_next)
            AddToRow(gradients, curves->GlobalIndex(j_prev), TPESC::tpe_grad(i_pt, j_pt, alpha, beta, j_prev));

        // If j_pt is not already included in one of the above 3
        if (j_pt != i_prev && j_pt != i_pt && j_pt != i_next)
            AddToRow(gradients, curves->GlobalIndex(j_pt), TPESC::tpe_grad(i_pt, j_pt, alpha, beta, j_pt));

        // If j_next is not already included in one of the above 3
        if (j_next != i_prev && j_next != i_pt && j_next != i_next)
            AddToRow(gradients, curves->GlobalIndex(j_next), TPESC::tpe_grad(i_pt, j_pt, alpha, beta, j_next));
    }

    void TPEFlowSolverSC::FillGradientVectorDirect(Eigen::MatrixXd &gradients) {
        int nVerts = curves->NumVertices();
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
        // The single energy term (i, j) affects six vertices:
        // (i_prev, i, i_next, j_prev, j, j_next).
        // We can restructure the computation as follows:
        // for each single 1-ring (i, i_prev, i_next), accumulate the
        // contributions from the gradients of both terms (i, j) and (j, i).
        int nVerts = curves->NumVertices();
        gradients.setZero();

        for (int i = 0; i < nVerts; i++) {
            PointOnCurve i_pt = curves->GetCurvePoint(i);
            root->accumulateTPEGradient(gradients, i_pt, curves, alpha, beta);
        }
    }

    void TPEFlowSolverSC::FillConstraintVector(Eigen::MatrixXd &gradients) {
        int nVerts = curves->NumVertices();

        for (int i = 0; i < gradients.rows(); i++) {
            SetRow(gradients, i, Vector3{0, 0, 0});
        }

        for (int i = 0; i < nVerts; i++) {
            PointOnCurve i_pt = curves->GetCurvePoint(i);
            Vector3 len_grad = i_pt.curve->TotalLengthGradient(i_pt.pIndex);
            SetRow(gradients, curves->GlobalIndex(i_pt), len_grad);
        }
    }

    bool TPEFlowSolverSC::StepNaive(double h) {
        // Takes a fixed time step h using the L2 gradient
        int nVerts = curves->NumVertices();
        Eigen::MatrixXd gradients(nVerts, 3);
        gradients.setZero();
        FillGradientVectorDirect(gradients);

        for (int i = 0; i < nVerts; i++) {
            PointOnCurve pt = curves->GetCurvePoint(i);
            pt.SetPosition(pt.Position() - h * SelectRow(gradients, i));
        }
        return true;
    }

    void TPEFlowSolverSC::SaveCurrentPositions() {
        for (int i = 0; i < curves->NumVertices(); i++) {
            originalPositions[i] = curves->GetCurvePoint(i).Position();
        }
    }

    void TPEFlowSolverSC::RestoreOriginalPositions() {
        for (int i = 0; i < curves->NumVertices(); i++) {
            curves->GetCurvePoint(i).SetPosition(originalPositions[i]);
        }
    }

    void TPEFlowSolverSC::SetGradientStep(Eigen::MatrixXd gradient, double delta) {
        // Write the new vertex positions to the mesh
        // Step every vertex by the gradient times delta
        for (int i = 0; i < curves->NumVertices(); i++) {
            Vector3 vertGrad = SelectRow(gradient, i);
            curves->GetCurvePoint(i).SetPosition(originalPositions[i] - delta * vertGrad);
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

        std::cout << "Initial energy " << initialEnergy << std::endl;

        // // PlotEnergyInDirection(gradient, sigma * gradDot);

        while (delta > ls_step_threshold) {
            SetGradientStep(gradient, delta);
            if (root) {
                // Update the centers of mass to reflect the new positions
                root->recomputeCentersOfMass(curves);
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
                //std::cout << "Took step of size " << delta << " after " << numBacktracks
                //    << " backtracks, " << numDoubles << " doubles" << std::endl;
                //std::cout << "Decreased energy by " << decrease << " (target was " << targetDecrease << ")" << std::endl;
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
            std::cout << "Energy: " << initialEnergy << " -> " << newEnergy << std::endl;
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
                root->recomputeCentersOfMass(curves);
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
        int nVerts = curves->NumVertices();
        Eigen::MatrixXd gradients(nVerts, 3);
        gradients.setZero();
        FillGradientVectorDirect(gradients);
        return LineSearchStep(gradients);
    }

    double TPEFlowSolverSC::ProjectSoboSloboGradient(Eigen::PartialPivLU<Eigen::MatrixXd> &lu, Eigen::MatrixXd &gradients) {
        int nVerts = curves->NumVertices();
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

    void TPEFlowSolverSC::SoboSloboMatrix(Eigen::MatrixXd &A) {
        int nVerts = curves->NumVertices();
        // Fill the top-left block with the gram matrix
        SobolevCurves::FillGlobalMatrix(curves, alpha, beta, A);

        double sumLength = curves->TotalLength();
        double sumW = 0;

        // Fill the bottom row with weights for the constraint
        for (int i = 0; i < nVerts; i++) {
            double areaWeight = curves->GetCurvePoint(i).DualLength() / sumLength;
            sumW += areaWeight;
            // Fill in bottom row and rightmost column
            A(i, nVerts) = areaWeight;
            A(nVerts, i) = areaWeight;
        }
    }

    template<typename T>
    void TestMultiply(T &mult, Eigen::VectorXd &xVec, Eigen::VectorXd &result) {
        mult->Multiply(xVec, result);
    }

    PolyCurve* TPEFlowSolverSC::TestCoarsen(PolyCurve* c) {
        if (c) {
            return c->Coarsen();
        }
        else {
            return curves->curves[0]->Coarsen();
        }
    }

    void TPEFlowSolverSC::CompareMatrixVectorProduct() {
        int nVerts = curves->NumVertices();
        int vecLen = nVerts + 1;

        Eigen::MatrixXd grads;
        grads.setZero(vecLen, 3);
        FillGradientVectorDirect(grads);
        Eigen::VectorXd x = grads.col(0);

        long matrixStart = Utils::currentTimeMilliseconds();
        Eigen::MatrixXd A;
        A.setZero(vecLen, vecLen);
        SoboSloboMatrix(A);
        long matrixEnd = Utils::currentTimeMilliseconds();
        //std::cerr << "Matrix assembly: " << (matrixEnd - matrixStart) << " ms" << std::endl;

        Eigen::VectorXd dense_solution = A.partialPivLu().solve(x);
        long matrixMultEnd = Utils::currentTimeMilliseconds();
        //std::cerr << "Matrix solve: " << (matrixMultEnd - matrixEnd) << " ms" << std::endl;

        long treeStart = Utils::currentTimeMilliseconds();

        BVHNode3D* bvh = CreateEdgeBVHFromCurve(curves);
        BlockClusterTree* tree = new BlockClusterTree(curves, bvh, 0.5, alpha, beta);
        Eigen::VectorXd tree_solution_v;
        tree_solution_v.setZero(vecLen);
        
        long treeEnd = Utils::currentTimeMilliseconds();
        //std::cerr << "Tree assembly: " << (treeEnd - treeStart) << " ms" << std::endl;

        // tree->CompareBlocks();
        tree->SetBlockTreeMode(BlockTreeMode::Barycenter);
        HMatrix matRepl(tree, vecLen);
        //Eigen::GMRES<HMatrix, Eigen::IdentityPreconditioner> gmres;
        //gmres.compute(matRepl);
        A(nVerts, nVerts) = 1;
        Eigen::GMRES<Eigen::MatrixXd, Eigen::IdentityPreconditioner> gmres(A);

        long mult_start = Utils::currentTimeMilliseconds();
        Eigen::VectorXd dense_mult = A * x;
        long mult_time = Utils::currentTimeMilliseconds() - mult_start;

        long tree_mult_start = Utils::currentTimeMilliseconds();
        Eigen::VectorXd tree_mult = matRepl * x;
        long tree_mult_time = Utils::currentTimeMilliseconds() - tree_mult_start;

        std::cerr << "Dense multiplication: " << mult_time << " ms" << std::endl;
        std::cerr << "Tree multiplication: " << tree_mult_time << " ms" << std::endl;

        double mult_relative = 100 * ((dense_mult - tree_mult).norm() / dense_mult.norm());
        std::cerr << "Mult. relative error = " << mult_relative << " percent" << std::endl;

        int maxIters = vecLen * 2;

        std::cerr << "* Approximate maximum time for iterative solve: " << (tree_mult_time * maxIters) << " ms" << std::endl;

        long iterStart = Utils::currentTimeMilliseconds();

        for (int i = 0; i < maxIters / 100; i++) {
            int curMax = (i + 1) * 100;
            gmres.setMaxIterations(curMax);
            tree_solution_v = gmres.solve(x);
            std::cerr << "Done with " << curMax << " iterations" << std::endl;

            double norm_diff = (dense_solution - tree_solution_v).norm();
            double norm_full = dense_solution.norm();
            double relative = 100 * norm_diff / norm_full;

            double residual = (x - A * tree_solution_v).norm() / x.norm();

            std::cout << curMax << ", " << residual << ", " << relative << std::endl;
        }

        std::cerr << "\nGMRES iterations: " << gmres.iterations() << std::endl;
        std::cerr << "Estimated error: " << gmres.error() << std::endl;

        long iterEnd = Utils::currentTimeMilliseconds();
        std::cerr << "Iterative solve: " << (iterEnd - iterStart) << " ms" << std::endl;

        Eigen::VectorXd tree_solution;
        tree_solution.setZero(vecLen);
        for (int i = 0; i < vecLen; i++) {
            tree_solution(i) = tree_solution_v[i];
        }

        Eigen::VectorXd diff = dense_solution - tree_solution;
        double norm_diff = diff.norm();
        double norm_full = dense_solution.norm();
        double relative = 100 * norm_diff / norm_full;

        std::cerr << "Norm dense solve = " << norm_full << "; norm iterative solve = " << tree_solution.norm() << std::endl;
        std::cerr << "Solve relative error = " << relative << " percent" << std::endl;

        delete tree;
        delete bvh;
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

    void TPEFlowSolverSC::FillVertLengthConstraintMatrix(Eigen::MatrixXd &A, int baseIndex) {
        int nVerts = curves->NumVertices();
        for (int i = 0; i < nVerts; i++) {
            PointOnCurve pt1 = curves->GetCurvePoint(i);
            PointOnCurve pt2 = pt1.Next();
            // This is the gradient of edge length wrt pt1; the gradient wrt pt2 is just negative of this.
            Vector3 grad1 = pt1.Position() - pt2.Position();
            grad1 = grad1.normalize();

            int j1 = curves->GlobalIndex(pt1);
            int j2 = curves->GlobalIndex(pt2);
            int curRow = baseIndex + i;

            // Write the three gradient entries for pt1 into the row and column
            A(curRow, 3 * j1    ) = grad1.x;
            A(curRow, 3 * j1 + 1) = grad1.y;
            A(curRow, 3 * j1 + 2) = grad1.z;
            A(3 * j1,     curRow) = grad1.x;
            A(3 * j1 + 1, curRow) = grad1.y;
            A(3 * j1 + 2, curRow) = grad1.z;

            // Similarly write the three gradient entries for pt2 into the same row and column
            A(curRow, 3 * j2    ) = -grad1.x;
            A(curRow, 3 * j2 + 1) = -grad1.y;
            A(curRow, 3 * j2 + 2) = -grad1.z;
            A(3 * j2,     curRow) = -grad1.x;
            A(3 * j2 + 1, curRow) = -grad1.y;
            A(3 * j2 + 2, curRow) = -grad1.z;
        }
    }

    double TPEFlowSolverSC::FillLengthConstraintViolations(Eigen::VectorXd &b, int baseIndex) {
        if (!useEdgeLengthConstraint) return 0;
        int nVerts = curves->NumVertices();

        double maxViolation = 0;

        for (int i = 0; i < nVerts; i++) {
            // Compute current segment length
            PointOnCurve p = curves->GetCurvePoint(i);
            Vector3 p1 = p.Position();
            Vector3 p2 = p.Next().Position();
            double curDist = norm(p1 - p2);

            double negError = initialLengths[i] - curDist;
            b(baseIndex + i) = negError;

            maxViolation = fmax(maxViolation, fabs(negError));
        }
        return maxViolation;
    }

    bool TPEFlowSolverSC::BackprojectConstraints(Eigen::PartialPivLU<Eigen::MatrixXd> &lu) {
        int nVerts = curves->NumVertices();
        Eigen::VectorXd b;
        b.setZero(matrixNumRows());

        Vector3 barycenter = curves->Barycenter();

        if (useEdgeLengthConstraint) {
            // If using per-edge length constraints, matrix has all coordinates merged,
            // so we only need one solve.

            // Place all 3 barycenter coordinates on RHS
            b(3 * nVerts) = -barycenter.x;
            b(3 * nVerts + 1) = -barycenter.y;
            b(3 * nVerts + 2) = -barycenter.z;

            double maxViolation = 0;

            // Add length violations to RHS
            maxViolation = FillLengthConstraintViolations(b, 3 * nVerts + 3);
            // Solve for correction
            Eigen::VectorXd corr = lu.solve(b);
            // Apply correction
            for (int i = 0; i < nVerts; i++) {
                Vector3 correction{corr(3 * i), corr(3 * i + 1), corr(3 * i + 2)};
                PointOnCurve pt = curves->GetCurvePoint(i);
                pt.SetPosition(pt.Position() + correction);
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
                PointOnCurve pt = curves->GetCurvePoint(i);
                pt.SetPosition(pt.Position() + correction);
            }

            return true;
        }
    }

    int TPEFlowSolverSC::matrixNumRows() {
        if (useEdgeLengthConstraint) {
            int nVerts = curves->NumVertices();
            // 3V rows for all vertex values; 3 rows for barycenter constraint; V rows for edge length constraints
            return 3 * nVerts + 3 + nVerts;
        }
        else {
            // V rows for vertex values by coordinate; 1 row for barycenter constraints by coordinate
            return curves->NumVertices() + 1;
        }
    }

    double TPEFlowSolverSC::ComputeAndProjectGradient(Eigen::MatrixXd &gradients) {
        Eigen::MatrixXd A;
        Eigen::PartialPivLU<Eigen::MatrixXd> lu;
        return ComputeAndProjectGradient(gradients, A, lu);
    }

    double TPEFlowSolverSC::ComputeAndProjectGradient(Eigen::MatrixXd &gradients, Eigen::MatrixXd &A, Eigen::PartialPivLU<Eigen::MatrixXd> &lu) {
        int nRows = matrixNumRows();
        size_t nVerts = curves->NumVertices();

        for (int i = 0; i < gradients.rows(); i++) {
            l2gradients.row(i) = gradients.row(i);
        }

        // If we're not using the per-edge length constraints, use total length instead
        if (!useEdgeLengthConstraint) {
            double ss_start = Utils::currentTimeMilliseconds();
            // Fill the Sobo-Slobo matrix, one entry per vertex
            A.setZero(nRows, nRows);
            SoboSloboMatrix(A);
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
            SoboSloboMatrix(A_temp);

            // Duplicate the saddle matrix portion into the larger matrix
            ExpandMatrix3x(A_temp, A);
            // Add rows for edge length constraint
            FillVertLengthConstraintMatrix(A, 3 * nVerts + 3);
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

        size_t nVerts = curves->NumVertices();

        int nRows = matrixNumRows();

        long grad_start = Utils::currentTimeMilliseconds();
        SpatialTree *tree_root = 0;

        // Barnes-Hut for gradient accumulation
        if (useBH) {
            tree_root = CreateBVHFromCurve(curves);
            FillGradientVectorBH(tree_root, vertGradients);
        }
        else {
            FillGradientVectorDirect(vertGradients);
        }
        
        std::cout << "\n====== Timing ======" << std::endl;
        double grad_end = Utils::currentTimeMilliseconds();

        std::cout << "  Assemble gradient: " << (grad_end - grad_start) << " ms" << std::endl;

        double length1 = curves->TotalLength();

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

        double length2 = curves->TotalLength();
        std::cout << "Length " << length1 << " -> " << length2 << std::endl;
        long end = Utils::currentTimeMilliseconds();
        std::cout << "Time = " << (end - start) << " ms" << std::endl;

        return step_size > 0;
    }

    double TPEFlowSolverSC::ProjectGradientIterative(Eigen::MatrixXd &gradients, BlockClusterTree* &blockTree) {

        for (int i = 0; i < gradients.rows(); i++) {
            l2gradients(i) = gradients(i);
        }

        // If we're not using the per-edge length constraints, use total length instead
        if (!useEdgeLengthConstraint) {
            // First solve Ax = b, where b is the gradients
            blockTree->SetBlockTreeMode(BlockTreeMode::Barycenter);
            std::vector<Vector3> sobolevGradients(gradients.rows());
            // TODO: Solve blockTree * sobolevGradients = gradients
            // Now compute gradient of total length constraint
            // FillConstraintVector(vertConstraints);
            // Solve Ax = l, where l is length gradient
            std::vector<Vector3> sobolevLengthGrads(vertConstraints.rows());
            // TODO: Solve blockTree * sobolevLengthGrads = vertConstraints
            blockTree->SetBlockTreeMode(BlockTreeMode::MatrixOnly);
            // TODO: Project out component of Sobolev length gradient from Sobolev energy gradient

        }
        else {
            // TODO: Per-edge case
            return 0;
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

    bool TPEFlowSolverSC::StepSobolevLSIterative() {
        size_t nVerts = curves->NumVertices();

        int nRows = matrixNumRows();

        SpatialTree *tree_root = 0;

        // Assemble the L2 gradient
        tree_root = CreateBVHFromCurve(curves);
        FillGradientVectorBH(tree_root, vertGradients);

        // Make a block cluster tree
        BVHNode3D *edgeTree = CreateEdgeBVHFromCurve(curves);
        BlockClusterTree* blockTree = new BlockClusterTree(curves, edgeTree, 0.25, alpha, beta);

        // Use tree to compute the Sobolev gradient
        double dot_acc = ProjectGradientIterative(vertGradients, blockTree);

        // Take a line search step using this gradient
        double ls_start = Utils::currentTimeMilliseconds();
        // double step_size = LineSearchStep(vertGradients, dot_acc, tree_root);
        double ls_end = Utils::currentTimeMilliseconds();
        std::cout << "  Line search: " << (ls_end - ls_start) << " ms" << std::endl;

        // Correct for drift with backprojection

        delete edgeTree;
        delete blockTree;
        return false;
    }
}
