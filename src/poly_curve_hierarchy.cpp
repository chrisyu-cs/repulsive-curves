#include "poly_curve_hierarchy.h"
#include "spatial/tpe_bvh.h"
#include "product/block_cluster_tree.h"
#include "product/matrix_free.h"

namespace LWS {

    PolyCurveGroupHierarchy::PolyCurveGroupHierarchy(PolyCurveGroup* topLevel, size_t numLevels) {
        levels.push_back(topLevel);

        while (levels.size() < numLevels) {
            AddNextLevel();
        }
    }

    PolyCurveGroupHierarchy::~PolyCurveGroupHierarchy() {
        for (size_t i = 1; i < levels.size(); i++) {
            delete levels[i];
        }
    }

    void PolyCurveGroupHierarchy::AddNextLevel() {
        PolyCurveGroup* lastLevel = levels[levels.size() - 1];
        MultigridOperator prolongOp, sparsifyOp;
        PolyCurveGroup* nextLevel = lastLevel->Coarsen(prolongOp, sparsifyOp);
        prolongOp.lowerSize = nextLevel->NumVertices();
        prolongOp.upperSize = lastLevel->NumVertices();

        std::cout << "Added multigrid level with " << nextLevel->NumVertices() << std::endl;

        prolongationOps.push_back(prolongOp);
        sparsifyOps.push_back(sparsifyOp);
        levels.push_back(nextLevel);
    }

    Eigen::VectorXd PolyCurveGroupHierarchy::VCycleSolve(Eigen::VectorXd b, double sepCoeff, double alpha, double beta) {
        int numLevels = levels.size();
        std::vector<BVHNode3D*> vertexBVHs(numLevels);
        std::vector<BVHNode3D*> edgeBVHs(numLevels);
        std::vector<BlockClusterTree*> trees(numLevels);
        std::vector<HMatrix> hMatrices(numLevels);

        for (int i = 0; i < numLevels; i++) {
            vertexBVHs[i] = CreateBVHFromCurve(levels[i]);
            edgeBVHs[i] = CreateEdgeBVHFromCurve(levels[i]);
            trees[i] = new BlockClusterTree(levels[i], edgeBVHs[i], sepCoeff, alpha, beta);
            // For now, just solve the problem with only the top-left block
            trees[i]->SetBlockTreeMode(BlockTreeMode::MatrixOnly);
            hMatrices[i] = HMatrix(trees[i], levels[i]->NumVertices());
        }

        std::vector<Eigen::VectorXd> residuals(numLevels);
        std::vector<Eigen::VectorXd> solutions(numLevels);

        residuals[0] = b;
        bool done = false;
        int numIters = 0;

        int coarsestRows = levels[numLevels - 1]->NumVertices();
        Eigen::MatrixXd coarse_A(coarsestRows, coarsestRows);
        coarse_A.setZero();

        while (!done) {
            // Propagate residuals downward
            for (int i = 0; i < numLevels - 1; i++) {
                Eigen::GMRES<HMatrix, Eigen::IdentityPreconditioner> gmres;
                gmres.compute(hMatrices[i]);
                gmres.setMaxIterations(30);
                // Run the smoother on this level to get some intermediate solution
                if (numIters == 0) {
                    solutions[i] = gmres.solve(residuals[i]);
                }
                else {
                    solutions[i] = gmres.solveWithGuess(residuals[i], solutions[i]);
                }
                // On the next level, the right-hand side becomes the restricted    
                // residual from this level.
                residuals[i + 1] = prolongationOps[i].mapDownward(residuals[i] - hMatrices[i] * solutions[i], MultigridMode::Barycenter);
            }

            // Assemble full operator on sparsest curve and solve normally
            solutions[numLevels - 1].setZero(residuals[numLevels - 1].rows());
            SobolevCurves::FillGlobalMatrix(levels[numLevels - 1], alpha, beta, coarse_A);
            for (int i = 0; i < coarsestRows; i++) {
                coarse_A(i, i) += 0.1 * levels[numLevels - 1]->GetCurvePoint(i).DualLength();
            }
            solutions[numLevels - 1] = coarse_A.partialPivLu().solve(residuals[numLevels - 1]);

            // Propagate solutions upward
            for (int i = numLevels - 2; i >= 0; i--) {
                Eigen::GMRES<HMatrix, Eigen::IdentityPreconditioner> gmres;
                gmres.compute(hMatrices[i]);
                gmres.setMaxIterations(30);
                // Compute the new initial guess -- old solution from this level, plus
                // new solution from prolongation operator
                Eigen::VectorXd guess = solutions[i] + prolongationOps[i].mapUpward(solutions[i + 1], MultigridMode::Barycenter);
                solutions[i] = gmres.solveWithGuess(residuals[i], guess);
            }

            double overallResidual = (b - hMatrices[0] * solutions[0]).lpNorm<Eigen::Infinity>();

            done = (overallResidual < 1e-4 || numIters > 100);
            numIters++;
            std::cout << "[Iteration " << numIters << "] residual = " << overallResidual << std::endl;
        }

        std::cout << "Did " << numIters << " iterations" << std::endl; 
        return solutions[0];
    }
}
