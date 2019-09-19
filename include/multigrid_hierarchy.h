#pragma once

#include "multigrid_domain.h"
#include "geometrycentral/utilities/vector3.h"

namespace LWS {

    template<typename T, typename Mult>
    class MultigridHierarchy {
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

            std::cout << "Added multigrid level with " << nextLevel->NumVertices() << std::endl;

            prolongationOps.push_back(prolongOp);
            restrictionOps.push_back(restrictOp);
            levels.push_back(nextLevel);
        }

        Eigen::VectorXd VCycleInitGuess(Eigen::VectorXd b, MultigridMode mode) {
            int numLevels = levels.size();
            Eigen::VectorXd coarseB = b;

            // Propagate RHS downward
            for (int i = 0; i < numLevels - 1; i++) {
                coarseB = restrictionOps[i].mapDownward(coarseB, mode);
            }

            Eigen::MatrixXd coarse_A = levels[numLevels - 1]->GetFullMatrix();
            Eigen::VectorXd sol = coarse_A.partialPivLu().solve(coarseB);

            // Propagate solution upward
            for (int i = numLevels - 2; i >= 0; i--) {
                sol = prolongationOps[i].mapUpward(sol, mode);
            }

            return sol;
    }

        // Solve Gx = b, where G is the Sobolev Gram matrix of the top-level curve.
        Eigen::VectorXd VCycleSolve(Eigen::VectorXd b, MultigridMode mode) {
            int numLevels = levels.size();
            std::vector<Eigen::VectorXd> residuals(numLevels);
            std::vector<Eigen::VectorXd> solutions(numLevels);
            std::vector<WrappedMatrix> hMatrices(numLevels);

            for (int i = 0; i < numLevels; i++) {
                int levelNVerts = levels[i]->NumVertices();
                if (mode == MultigridMode::Barycenter) levelNVerts++;
                hMatrices[i] = WrappedMatrix(levels[i]->GetMultiplier(), levelNVerts);
            }

            residuals[0] = b;
            bool done = false;
            int numIters = 0;

            int coarsestRows = levels[numLevels - 1]->NumVertices();
            Eigen::MatrixXd coarse_A;
            coarse_A.setZero();

            int iter = 0;

            Eigen::VectorXd initGuess = VCycleInitGuess(b, mode);
            // std::cout << "Initial guess: " << initGuess.transpose() << std::endl;

            while (!done) {
                // Propagate residuals downward
                for (int i = 0; i < numLevels - 1; i++) {
                    //Eigen::GMRES<HMatrix, Eigen::IdentityPreconditioner> gmres;
                    Eigen::BiCGSTAB<Product::MatrixReplacement<Mult>, Eigen::IdentityPreconditioner> gmres;
                    // Eigen::ConjugateGradient<Product::MatrixReplacement<Mult>, Eigen::Lower|Eigen::Upper, Eigen::IdentityPreconditioner> gmres;
                    gmres.compute(hMatrices[i]);

                    gmres.setMaxIterations(24);
                    // Run the smoother on this level to get some intermediate solution
                    if (numIters == 0) {
                        if (i == 0) {
                            solutions[i] = gmres.solveWithGuess(residuals[i], initGuess);
                        }
                        else {
                            solutions[i] = gmres.solve(residuals[i]);
                        }
                        // solutions[i] = gmres.solve(residuals[i]);
                    }
                    else {
                        // Eigen::VectorXd resid_i = residuals[i] - hMatrices[i] * solutions[i];
                        // std::cout << iter++ << ", " << i << ", " << resid_i.lpNorm<Eigen::Infinity>() << std::endl;
                        solutions[i] = gmres.solveWithGuess(residuals[i], solutions[i]);
                    }
                    // On the next level, the right-hand side becomes the restricted    
                    // residual from this level.
                    Eigen::VectorXd resid_i = residuals[i] - hMatrices[i] * solutions[i];
                    // std::cout << iter++ << ", " << i << ", " << resid_i.lpNorm<Eigen::Infinity>() << std::endl;
                    residuals[i + 1] = prolongationOps[i].mapDownward(resid_i, mode);
                }

                // Assemble full operator on sparsest curve and solve normally
                solutions[numLevels - 1].setZero(residuals[numLevels - 1].rows());
                coarse_A = levels[numLevels - 1]->GetFullMatrix();
                solutions[numLevels - 1] = coarse_A.partialPivLu().solve(residuals[numLevels - 1]);

                Eigen::VectorXd resid_bottom = residuals[numLevels - 1] - hMatrices[numLevels - 1] * solutions[numLevels - 1];
                // std::cout << iter++ << ", " << (numLevels - 1) << ", " << resid_bottom.lpNorm<Eigen::Infinity>() << std::endl;

                // Propagate solutions upward
                for (int i = numLevels - 2; i >= 0; i--) {
                    //Eigen::GMRES<HMatrix, Eigen::IdentityPreconditioner> gmres;
                    Eigen::BiCGSTAB<Product::MatrixReplacement<Mult>, Eigen::IdentityPreconditioner> gmres;
                    // Eigen::ConjugateGradient<Product::MatrixReplacement<Mult>, Eigen::Lower|Eigen::Upper, Eigen::IdentityPreconditioner> gmres;
                    gmres.compute(hMatrices[i]);
                    gmres.setMaxIterations(24);
                    // Compute the new initial guess -- old solution from this level, plus
                    // new solution from prolongation operator
                    Eigen::VectorXd guess = solutions[i] + prolongationOps[i].mapUpward(solutions[i + 1], mode);

                    Eigen::VectorXd resid_before = residuals[i] - hMatrices[i] * guess;
                    // std::cout << iter++ << ", " << i << ", " << resid_before.lpNorm<Eigen::Infinity>() << std::endl;

                    solutions[i] = gmres.solveWithGuess(residuals[i], guess);

                    Eigen::VectorXd resid_i = residuals[i] - hMatrices[i] * solutions[i];
                    // std::cout << iter++ << ", " << i << ", " << resid_i.lpNorm<Eigen::Infinity>() << std::endl;
                }

                Eigen::VectorXd overallResidual = b - hMatrices[0] * solutions[0];
                double residNorm = overallResidual.lpNorm<Eigen::Infinity>();

                numIters++;
                done = (residNorm < 1e-5 || numIters >= 100);
                std::cerr << "[Iteration " << numIters << "] residual = " << residNorm << "     \r" << std::flush;
            }

            std::cout << "\nDid " << numIters << " iterations" << std::endl; 
            return solutions[0];
        }

    };
    
}
