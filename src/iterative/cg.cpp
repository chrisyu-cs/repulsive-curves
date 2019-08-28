#include "iterative/cg.h"
#include "utils.h"

namespace LWS {

    void ConjugateGradient::CGSolve(BlockClusterTree *tree, std::vector<double> &x, std::vector<double> &b) {
        size_t bsize = b.size();
        std::vector<double> r_k(bsize);
        std::vector<double> p_k(bsize);
        
        // Use p_k as temporary storage for A * x
        tree->Multiply(x, p_k);

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
            // std::cout << k << ", ";
            done = CGIteration(tree, x, r_k, p_k);
            k++;
        }

        // Once finished, final result is alread stored in x
        std::cout << "CG finished in " << k << " iterations" << std::endl;
    }

    void ConjugateGradient::CGSolveComponents(BlockClusterTree *tree, std::vector<Vector3> &x, std::vector<Vector3> &b) {
        std::vector<double> x_comp(x.size());
        std::vector<double> b_comp(b.size());

        // Solve system for x coordinate
        for (size_t i = 0; i < x.size(); i++) {
            x_comp[i] = x[i].x;
            b_comp[i] = b[i].x;
        }

        CGSolve(tree, x_comp, b_comp);

        for (size_t i = 0; i < b.size(); i++) {
            x[i].x = x_comp[i];
        }

        // Solve system for y coordinate
        for (size_t i = 0; i < x.size(); i++) {
            x_comp[i] = x[i].y;
            b_comp[i] = b[i].y;
        }

        CGSolve(tree, x_comp, b_comp);

        for (size_t i = 0; i < b.size(); i++) {
            x[i].y = x_comp[i];
        }

        // Solve system for z coordinate
        for (size_t i = 0; i < x.size(); i++) {
            x_comp[i] = x[i].z;
            b_comp[i] = b[i].z;
        }

        CGSolve(tree, x_comp, b_comp);

        for (size_t i = 0; i < b.size(); i++) {
            x[i].z = x_comp[i];
        }
    }

    bool ConjugateGradient::CGIteration(BlockClusterTree *tree, std::vector<double> &x_k,
    std::vector<double> &r_k, std::vector<double> &p_k) {
        // Save size of vector
        size_t size = x_k.size();
        std::vector<double> A_times_pk(size);

        // Compute alpha_k
        double rkT_rk = std_vector_dot(r_k, r_k);
        tree->Multiply(p_k, A_times_pk);
        double pkT_A_pk = std_vector_dot(p_k, A_times_pk);

        double alpha_k = rkT_rk / pkT_A_pk;

        // Update x and r
        for (size_t i = 0; i < size; i++) {
            x_k[i] = x_k[i] + alpha_k * p_k[i];
            r_k[i] = r_k[i] - alpha_k * A_times_pk[i];
        }

        // If residual is sufficiently small, we can stop
        double rk1T_rk1 = std_vector_dot(r_k, r_k);
        std::cout << "Residual = " << rk1T_rk1 << "\r" << std::flush;
        
        // std::cout << rk1T_rk1 << std::endl;
        
        if (rk1T_rk1 < 1e-3) {
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
