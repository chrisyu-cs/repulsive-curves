#pragma once

#include "multigrid_domain.h"
#include "geometrycentral/utilities/vector3.h"

namespace LWS {

    template<typename T, typename Mult>
    class MultigridHierarchy {
        public:
        std::vector<MultigridDomain<T, Mult>*> levels;
        std::vector<MultigridOperator> prolongationOps;
        std::vector<MultigridOperator> sparsifyOps;

        MultigridHierarchy(MultigridDomain<T, Mult>* topLevel, size_t numLevels) {
            levels.push_back(topLevel);
            std::cout << "Constructing hierarchy..." << std::endl;

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
            MultigridOperator prolongOp, sparsifyOp;
            MultigridDomain<T, Mult>* nextLevel = lastLevel->Coarsen(prolongOp, sparsifyOp);
            prolongOp.lowerSize = nextLevel->NumVertices();
            prolongOp.upperSize = lastLevel->NumVertices();

            std::cout << "Added multigrid level with " << nextLevel->NumVertices() << std::endl;

            prolongationOps.push_back(prolongOp);
            sparsifyOps.push_back(sparsifyOp);
            levels.push_back(nextLevel);
        }

        // Solve Gx = b, where G is the Sobolev Gram matrix of the top-level curve.
        Eigen::VectorXd VCycleSolve(Eigen::VectorXd b, double sepCoeff, double alpha, double beta, MultigridMode mode) {
            int numLevels = levels.size();
            std::vector<Eigen::VectorXd> residuals(numLevels);
            std::vector<Eigen::VectorXd> solutions(numLevels);

            residuals[0] = b;
            bool done = false;
            int numIters = 0;

            int coarsestRows = levels[numLevels - 1]->NumVertices();
            Eigen::MatrixXd coarse_A;
            coarse_A.setZero();

            while (!done) {
                // Propagate residuals downward
                for (int i = 0; i < numLevels - 1; i++) {
                    //Eigen::GMRES<HMatrix, Eigen::IdentityPreconditioner> gmres;
                    Eigen::GMRES<Product::MatrixReplacement<Mult>, Eigen::IdentityPreconditioner> gmres;
                    gmres.compute(levels[i]->GetMultiplier());

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
                    residuals[i + 1] = prolongationOps[i].mapDownward(residuals[i] - levels[i]->GetMultiplier() * solutions[i], mode);
                }

                // Assemble full operator on sparsest curve and solve normally
                solutions[numLevels - 1].setZero(residuals[numLevels - 1].rows());
                coarse_A = levels[numLevels - 1]->GetFullMatrix();
                solutions[numLevels - 1] = coarse_A.partialPivLu().solve(residuals[numLevels - 1]);

                // Propagate solutions upward
                for (int i = numLevels - 2; i >= 0; i--) {
                    //Eigen::GMRES<HMatrix, Eigen::IdentityPreconditioner> gmres;
                    Eigen::GMRES<Product::MatrixReplacement<Mult>, Eigen::IdentityPreconditioner> gmres;
                    gmres.compute(levels[i]->GetMultiplier());
                    gmres.setMaxIterations(30);
                    // Compute the new initial guess -- old solution from this level, plus
                    // new solution from prolongation operator
                    Eigen::VectorXd guess = solutions[i] + prolongationOps[i].mapUpward(solutions[i + 1], mode);
                    solutions[i] = gmres.solveWithGuess(residuals[i], guess);
                }

                Eigen::VectorXd overallResidual = b - levels[0]->GetMultiplier() * solutions[0];
                double residNorm = overallResidual.lpNorm<Eigen::Infinity>();

                numIters++;
                done = (residNorm < 1e-5 || numIters >= 100);
                std::cout << "[Iteration " << numIters << "] residual = " << residNorm << std::endl;
            }

            std::cout << "Did " << numIters << " iterations" << std::endl; 
            return solutions[0];
        }

    };
    
}
