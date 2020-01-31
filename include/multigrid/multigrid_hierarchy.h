#pragma once

#include "multigrid_domain.h"

namespace LWS {

    // The template argument should be a class that extends MultigridDomain.
    template<typename Domain>
    class MultigridHierarchy {
        public:
        using Mult = typename Domain::MultType;
        std::vector<MultigridDomain<Mult>*> levels;
        std::vector<MultigridOperator*> prolongationOps;

        MultigridHierarchy(MultigridDomain<Mult>* topLevel, size_t numLevels) {
            AddLevels(topLevel, numLevels);
        }

        MultigridHierarchy(MultigridDomain<Mult>* topLevel) {
            int nVerts = topLevel->NumVertices();
            int logNumVerts = log2(nVerts) - 4;
            logNumVerts = std::max(1, logNumVerts);
            AddLevels(topLevel, logNumVerts);
        }

        ~MultigridHierarchy() {
            for (size_t i = 0; i < levels.size(); i++) {
                delete levels[i];
            }

            for (size_t i = 0; i < prolongationOps.size(); i++) {
                delete prolongationOps[i];
            }
        }

        void AddLevels(MultigridDomain<Mult>* topLevel, size_t numLevels) {
            levels.push_back(topLevel);
            while (levels.size() < numLevels) {
                AddNextLevel();
            }
        }

        void AddNextLevel() {
            MultigridDomain<Mult>* lastLevel = levels[levels.size() - 1];
            MultigridOperator* prolongOp = lastLevel->MakeNewOperator();
            MultigridDomain<Mult>* nextLevel = lastLevel->Coarsen(prolongOp);

            prolongationOps.push_back(prolongOp);
            levels.push_back(nextLevel);
        }

        inline int NumRows() {
            return levels[0]->NumRows();
        }

        inline Mult* GetTopLevelMultiplier() const {
            return levels[0]->GetMultiplier();
        }

        template<typename Smoother>
        Eigen::VectorXd VCycleInitGuess(Eigen::VectorXd b, std::vector<Product::MatrixReplacement<Mult>> &hMatrices) {
            int numLevels = levels.size();
            Eigen::VectorXd coarseB = b;

            std::vector<Eigen::VectorXd> coarseBs(numLevels);
            coarseBs[0] = coarseB;

            // Propagate RHS downward
            for (int i = 0; i < numLevels - 1; i++) {
                coarseB = prolongationOps[i]->restrictWithTranspose(coarseB);
                coarseBs[i + 1] = coarseB;
            }

            Eigen::VectorXd sol = levels[numLevels - 1]->DirectSolve(coarseB);

            // Propagate solution upward
            for (int i = numLevels - 2; i >= 0; i--) {
                sol = prolongationOps[i]->prolong(sol);
                Smoother smoother;
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
        Eigen::VectorXd VCycleSolve(Eigen::VectorXd b, double tolerance) {
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

            Eigen::VectorXd initGuess = VCycleInitGuess<Smoother>(b, hMatrices);
            Eigen::VectorXd initResidual = b - hMatrices[0] * initGuess;
            double initRelative = initResidual.lpNorm<Eigen::Infinity>() / b.lpNorm<Eigen::Infinity>();
            // std::cout << "Initial residual = " << initRelative << std::endl;
            if (initRelative < tolerance) {
                // std::cout << "  * Initial residual low enough; no V-cycles needed" << std::endl;
                return initGuess;
            }

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
                        // std::cout << "First pre-smoothed relative residual = " << presmoothedRel << std::endl;
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
                    residuals[i + 1] = prolongationOps[i]->restrictWithTranspose(resid_i);
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
                    Eigen::VectorXd guess = solutions[i] + prolongationOps[i]->prolong(solutions[i + 1]);
                    solutions[i] = smoother.solveWithGuess(residuals[i], guess);
                }

                Eigen::VectorXd overallResidual = b - hMatrices[0] * solutions[0];
                double residNorm = overallResidual.lpNorm<Eigen::Infinity>() / b.lpNorm<Eigen::Infinity>();

                numIters++;
                done = (residNorm < tolerance || numIters >= 20);

                std::cerr << "  * [Iteration " << numIters << "] residual = " << residNorm << "     \r" << std::flush;
            }

            std::cout << "\n  * Finished in " << numIters << " iterations" << std::endl; 
            return solutions[0];
        }

    };
    
}
