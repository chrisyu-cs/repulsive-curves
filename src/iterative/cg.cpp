#include "iterative/cg.h"
#include "utils.h"

namespace LWS {

    void ConjugateGradient::CGSolve(BlockClusterTree &A, std::vector<double> &x, std::vector<double> &b) {
        size_t bsize = b.size();
        std::vector<double> r_k(bsize);
        std::vector<double> p_k(bsize);
        
        // Use p_k as temporary storage for A * x
        A.MultiplyVector(x, p_k);

        // Compute initial residual r_0
        for (size_t i = 0; i < bsize; i++) {
            r_k[i] = b[i] - p_k[i];
        }
        // Set p_0 = r_0
        for (size_t i = 0; i < bsize; i++) {
            p_k[i] = r_k[i];
        }

        int k = 0;
        bool done = false;
        while (!done) {
            CGIteration(A, x, r_k, p_k);
            k++;
        }

        // Once finished, final result is alread stored in x
        std::cout << "CG finished in " << k << " iterations" << std::endl;
    }

    bool ConjugateGradient::CGIteration(BlockClusterTree &A, std::vector<double> &x_k,
    std::vector<double> &r_k, std::vector<double> &p_k) {
        // Save size of vector
        size_t size = x_k.size();
        std::vector<double> A_times_pk(size);

        // Compute alpha_k
        double rkT_rk = std_vector_dot(r_k, r_k);
        A.MultiplyVector(p_k, A_times_pk);
        double pkT_A_pk = std_vector_dot(p_k, A_times_pk);

        double alpha_k = rkT_rk / pkT_A_pk;

        // Update x and r
        for (size_t i = 0; i < size; i++) {
            x_k[i] = x_k[i] + alpha_k * p_k[i];
            r_k[i] = r_k[i] - alpha_k * A_times_pk[i];
        }

        // If residual is sufficiently small, we can stop
        double rk1T_rk1 = std_vector_dot(r_k, r_k);
        if (rk1T_rk1 < 1e-6) {
            return true;
        }

        // Otherwise compute beta and update p
        double beta_k = rk1T_rk1 / rkT_rk;
        for (size_t i = 0; i < size; i++) {
            p_k[i] = r_k[i] + beta_k * p_k[i];
        }

        return false;
    }

}
