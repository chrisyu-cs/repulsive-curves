#include "preconditioned/iterative_solve.h"

#include <unsupported/Eigen/IterativeSolvers>

namespace LWS
{
    namespace preconditioned
    {
        void SolveIterative(PolyCurveNetwork *curves, const Eigen::VectorXd &v, Eigen::VectorXd &dst, double sep)
        {
            SparseHs *hs = new SparseHs(curves, 3, 6, sep);
            BCTMatrixReplacement fracL;
            fracL.addHs(hs);

            // Eigen::ConjugateGradient<BCTMatrixReplacement, Eigen::Lower | Eigen::Upper, CustomPreconditioner> cg;
            Eigen::GMRES<BCTMatrixReplacement, CustomPreconditioner> cg;
            cg.compute(fracL);

            Eigen::VectorXd temp;
            temp.setZero(v.rows());
            cg.setTolerance(1e-6);
            temp = cg.solveWithGuess(v, temp);
            std::cout << "  * GMRES converged in " << cg.iterations() << " iterations, final residual = " << cg.error() << std::endl;

            dst = temp;

            delete hs;
        }
    } // namespace preconditioned
} // namespace LWS
