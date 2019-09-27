#pragma once

#include "multigrid_domain.h"
#include "geometrycentral/utilities/vector3.h"

namespace LWS {

    template<typename T>
    class MultigridHierarchy {
        using Mult = typename T::MultType;

        public:
        std::vector<MultigridDomain<T, Mult>*> levels;
        std::vector<MultigridOperator> prolongationOps;
        std::vector<MultigridOperator> restrictionOps;

        MultigridHierarchy(MultigridDomain<T, Mult>* topLevel, size_t numLevels) {
            levels.push_back(topLevel);
            std::cout << "Constructing hierarchy with " << numLevels << " levels..." << std::endl;

            while (levels.size() < numLevels) {
                AddNextLevel();
            }
        }

        ~MultigridHierarchy() {
            for (size_t i = 0; i < levels.size(); i++) {
                delete levels[i];
            }
        }

        void AddNextLevel() {
            MultigridDomain<T, Mult>* lastLevel = levels[levels.size() - 1];
            MultigridOperator prolongOp, restrictOp;
            MultigridDomain<T, Mult>* nextLevel = lastLevel->Coarsen(prolongOp, restrictOp);
            prolongOp.lowerSize = nextLevel->NumVertices();
            prolongOp.upperSize = lastLevel->NumVertices();
            restrictOp.lowerSize = prolongOp.lowerSize;
            restrictOp.upperSize = prolongOp.upperSize;

            prolongOp.lowerP = nextLevel->GetConstraintProjector();
            prolongOp.upperP = lastLevel->GetConstraintProjector();
            restrictOp.lowerP = nextLevel->GetConstraintProjector();
            restrictOp.upperP = lastLevel->GetConstraintProjector();

            std::cout << "Added multigrid level with " << nextLevel->NumVertices() << std::endl;

            prolongationOps.push_back(prolongOp);
            restrictionOps.push_back(restrictOp);
            levels.push_back(nextLevel);
        }

        Eigen::VectorXd VCycleInitGuess(Eigen::VectorXd b, MultigridMode mode, std::vector<Product::MatrixReplacement<Mult>> &hMatrices) {
            int numLevels = levels.size();
            Eigen::VectorXd coarseB = b;

            std::vector<Eigen::VectorXd> coarseBs(numLevels);
            coarseBs[0] = coarseB;

            // Propagate RHS downward
            for (int i = 0; i < numLevels - 1; i++) {
                coarseB = prolongationOps[i].restrictWithTranspose(coarseB, mode);
                coarseBs[i + 1] = coarseB;
            }

            Eigen::VectorXd sol = levels[numLevels - 1]->DirectSolve(coarseB);

            // Propagate solution upward
            for (int i = numLevels - 2; i >= 0; i--) {
                sol = prolongationOps[i].prolong(sol, mode);
                Eigen::ConjugateGradient<Product::MatrixReplacement<Mult>, Eigen::Lower|Eigen::Upper, Eigen::IdentityPreconditioner> smoother;
                smoother.compute(hMatrices[i]);
                smoother.setMaxIterations(24);
                sol = smoother.solveWithGuess(coarseBs[i], sol);
            }

            return sol;
    }

        typedef Eigen::GMRES<Product::MatrixReplacement<Mult>, Eigen::IdentityPreconditioner> EigenGMRES;
        typedef Eigen::ConjugateGradient<Product::MatrixReplacement<Mult>, Eigen::Lower|Eigen::Upper, Eigen::IdentityPreconditioner> EigenCG;

        // Solve Gx = b, where G is the Sobolev Gram matrix of the top-level curve.
        template<typename Smoother>
        Eigen::VectorXd VCycleSolve(Eigen::VectorXd b) {
            MultigridMode mode = levels[0]->GetMode();

            int numLevels = levels.size();
            std::vector<Eigen::VectorXd> residuals(numLevels);
            std::vector<Eigen::VectorXd> solutions(numLevels);
            std::vector<Product::MatrixReplacement<Mult>> hMatrices(numLevels);

            for (int i = 0; i < numLevels; i++) {
                int levelNRows = levels[i]->NumRows();
                solutions[i].setZero(levelNRows);
                hMatrices[i] = Product::MatrixReplacement<Mult>(levels[i]->GetMultiplier(), levelNRows);
            }

            residuals[0] = b;
            bool done = false;
            int numIters = 0;

            Eigen::VectorXd initGuess = VCycleInitGuess(b, mode, hMatrices);
            Eigen::VectorXd initResidual = b - hMatrices[0] * initGuess;
            double initRelative = initResidual.lpNorm<Eigen::Infinity>() / b.lpNorm<Eigen::Infinity>();
            std::cout << "Initial guess relative residual = " << initRelative << std::endl;

            while (!done) {
                // Propagate residuals downward
                for (int i = 0; i < numLevels - 1; i++) {
                    Smoother smoother;
                    smoother.compute(hMatrices[i]);

                    smoother.setMaxIterations(12);
                    // Run the smoother on this level to get some intermediate solution
                    if (numIters == 0 && i == 0) {
                        // On the first level of the first V-cycle, use the initial guess
                        solutions[i] = smoother.solveWithGuess(residuals[i], initGuess);
                        Eigen::VectorXd presmoothedResid = b - hMatrices[0] * solutions[0];
                        double presmoothedRel = presmoothedResid.lpNorm<Eigen::Infinity>() / b.lpNorm<Eigen::Infinity>();
                        std::cout << "First pre-smoothed relative residual = " << presmoothedRel << std::endl;
                    }
                    else {
                        // On anything below the first level, zero out the initial guess
                        if (i != 0) solutions[i].setZero();
                        // Do pre-smoothing iterations
                        solutions[i] = smoother.solveWithGuess(residuals[i], solutions[i]);
                    }
                    // On the next level, the right-hand side becomes the restricted    
                    // residual from this level.
                    Eigen::VectorXd resid_i = residuals[i] - hMatrices[i] * solutions[i];
                    residuals[i + 1] = prolongationOps[i].restrictWithTranspose(resid_i, mode);
                }

                // Assemble full operator on sparsest curve and solve normally
                solutions[numLevels - 1].setZero(residuals[numLevels - 1].rows());
                solutions[numLevels - 1] = levels[numLevels - 1]->DirectSolve(residuals[numLevels - 1]);

                Eigen::VectorXd resid_bottom = residuals[numLevels - 1] - hMatrices[numLevels - 1] * solutions[numLevels - 1];

                // Propagate solutions upward
                for (int i = numLevels - 2; i >= 0; i--) {
                    Smoother smoother;
                    smoother.compute(hMatrices[i]);
                    smoother.setMaxIterations(12);
                    // Compute the new initial guess -- old solution from this level, plus
                    // new solution from prolongation operator
                    Eigen::VectorXd guess = solutions[i] + prolongationOps[i].prolong(solutions[i + 1], mode);
                    solutions[i] = smoother.solveWithGuess(residuals[i], guess);
                }

                Eigen::VectorXd overallResidual = b - hMatrices[0] * solutions[0];
                double residNorm = overallResidual.lpNorm<Eigen::Infinity>() / b.lpNorm<Eigen::Infinity>();

                numIters++;
                done = (residNorm < 1e-5 || numIters >= 100);
                std::cerr << "[Iteration " << numIters << "] residual = " << residNorm << "     \r" << std::flush;
            }

            std::cout << "\nDid " << numIters << " iterations" << std::endl; 
            return solutions[0];
        }

    };
    
}
