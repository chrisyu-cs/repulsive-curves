#pragma once

#include "../product/block_cluster_tree.h"
#include "utils.h"

namespace LWS {

    class BiCGStab {
        private:
        static std::vector<double> A_pk;
        static std::vector<double> s_k;
        static std::vector<double> A_sk;

        public:
        // Solve the system Ax = b, using whatever is in x as an initial guess.
        // Stores the solution in x once finished.
        template<typename T>
        static void Solve(VectorMultiplier<T> *A, std::vector<double> &x, std::vector<double> &b)
        {
            A_pk.resize(x.size());
            s_k.resize(x.size());
            A_sk.resize(x.size());

            std::vector<double> r_k(x.size());
            std::vector<double> r_0_prime(x.size());
            std::vector<double> p_k(x.size());

            // Temporarily store A*x_0 in r_k
            A->Multiply(x, r_k);

            // Now let r_0 = b - A*x_0
            for (size_t i = 0; i < r_k.size(); i++) {
                r_k[i] = b[i] - r_k[i];
            }

            for (size_t i = 0; i < r_k.size(); i++) {
                // Choose r_0' such that r_0 * r_0' != 0;
                // an easy way to do this is just r_0' = r_0
                r_0_prime[i] = r_k[i];
                // Set p_0 = r_0
                p_k[i] = r_k[i];
            }

            // Now do the iterations
            int k = 0;
            bool done = false;
            while (!done) {
                done = Iteration(A, x, r_k, p_k, r_0_prime, 1e-6);
                k++;
            }

            // Once finished, final result is alread stored in x
            std::cout << "\nBiCGStab finished in " << k << " iterations" << std::endl;
        }
        
        template<typename T>
        static bool Iteration(VectorMultiplier<T> *A, std::vector<double> &x_k,
            std::vector<double> &r_k, std::vector<double> &p_k, std::vector<double> &r_0_prime, double epsilon)
        {
            // Precompute A * p_k
            A->Multiply(p_k, A_pk);
            // Precompute r_k * r_0 as well
            double rk_r0 = std_vector_dot(r_k, r_0_prime);

            // alpha_k = (r_j * r_0') / ((A * p_k) * r_0')
            double alpha_k = rk_r0 / std_vector_dot(A_pk, r_0_prime);
            
            // s_k = r_k - alpha_k A p_k
            for (size_t i = 0; i < s_k.size(); i++) {
                s_k[i] = r_k[i] - alpha_k * A_pk[i];
            }

            // If ||s_k|| is small, stop
            double norm_sk = sqrt(std_vector_dot(s_k, s_k));
            if (norm_sk < epsilon) {
                // x_{k+1} = x_k + alpha_k * p_k 
                for (size_t i = 0; i < x_k.size(); i++) {
                    x_k[i] = x_k[i] + alpha_k * p_k[i];
                }
                return true;
            }

            // Precompute A * s_j
            A->Multiply(s_k, A_sk);

            // omega_k = ((A * s_j) * s_j) / ((A * s_j) * (A * s_j))
            double omega_k = std_vector_dot(A_sk, s_k) / std_vector_dot(A_sk, A_sk);

            // Update x_k in place, since this is the last time we need to refer to it
            // x_{k+1} = x_k + alpha_k * p_k + omega_k * s_k
            for (size_t i = 0; i < x_k.size(); i++) {
                x_k[i] = x_k[i] + alpha_k * p_k[i] + omega_k * s_k[i];
            }

            // Update r_k in place as well
            // r_{k+1} = s_k - omega_k * A * s_k
            for (size_t i = 0; i < r_k.size(); i++) {
                r_k[i] = s_k[i] - omega_k * A_sk[i];
            }

            // If ||r_{k+1}|| < epsilon, break
            double norm_rk1 = sqrt(std_vector_dot(r_k, r_k));
            std::cout << "Residual = " << norm_rk1 << ", rk * r0 = " << rk_r0 << "          \r" << std::flush;
            if (norm_rk1 < epsilon) {
                return true;
            }

            // beta_k = (alpha_k / omega_k) * (r_{k+1} * r_0') / (r_k * r_0')
            double rk1_r0 = std_vector_dot(r_k, r_0_prime);
            double beta_k = (alpha_k / omega_k) * (rk1_r0) / (rk_r0);

            // Update p_k in place
            // p_{k+1} = r_{k+1} + beta_k * (p_k - omega_j * A * p_k)
            for (size_t i = 0; i < p_k.size(); i++) {
                p_k[i] = r_k[i] + beta_k * (p_k[i] - omega_k * A_pk[i]);
            }

            // Restart procedure: if |r_{k+1} * r_k| is too small, then reset
            // the r_0 values
            if (fabs(rk1_r0) < 1e-6) {
                for (size_t i = 0; i < p_k.size(); i++) {
                    // r_0' = r_{k+1}
                    r_0_prime[i] = r_k[i];
                    // p_{k+1} = r_{k+1}
                    p_k[i] = r_k[i];
                }
            }
            return false;
        }
    };

}
