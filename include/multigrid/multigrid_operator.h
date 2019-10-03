#pragma once

#include <vector>
#include <Eigen/Core>
#include <Eigen/Sparse>
#include "nullspace_projector.h"

namespace LWS {
    struct IndexedMatrix {
        Eigen::SparseMatrix<double> M;
        int fineOffset;
        int coarseOffset;
    };

    enum class ProlongationMode {
        MatrixOnly,
        MatrixAndProjector,
        Barycenter,
        Matrix3AndProjector,
        Matrix3AndBarycenter,
        Matrix3AndEdgeConstraints
    };
    
    template<typename V>
    V ApplyPinv(Eigen::SparseMatrix<double> &J, V &x) {
        Eigen::SparseMatrix<double> JTJ = J.transpose() * J;
        Eigen::SparseLU<Eigen::SparseMatrix<double>, Eigen::COLAMDOrdering<int>> JTJ_solver;
        JTJ_solver.analyzePattern(JTJ);
        JTJ_solver.factorize(JTJ);

        return JTJ_solver.solve(J.transpose() * x);
    }

    template<> inline
    Eigen::MatrixXd ApplyPinv(Eigen::SparseMatrix<double> &J, Eigen::MatrixXd &xs) {
        Eigen::MatrixXd out(J.cols(), xs.cols());
        for (int i = 0; i < xs.cols(); i++) {
            Eigen::VectorXd c_i = xs.col(i);
            out.col(i) = ApplyPinv(J, c_i);
        }
        return out;
    }

    class MultigridOperator {
        public:
        MultigridOperator();
        int lowerSize;
        int upperSize;
        NullSpaceProjector* lowerP;
        NullSpaceProjector* upperP;
        
        std::vector<IndexedMatrix> matrices;
        std::vector<IndexedMatrix> edgeMatrices;
        
        Eigen::VectorXd prolong(Eigen::VectorXd v, ProlongationMode mode);
        Eigen::VectorXd restrictWithTranspose(Eigen::VectorXd v, ProlongationMode mode);
        Eigen::VectorXd restrictWithPinv(Eigen::VectorXd v, ProlongationMode mode);

        private:
        void checkAndSetOutputSize(Eigen::VectorXd &in, Eigen::VectorXd &out, int upperSize, int lowerSize, ProlongationMode mode);

        template<typename V, typename Dest>
        void prolongVerts1X(V &in, Dest &out) {
            for (size_t i = 0; i < matrices.size(); i++) {
                int outputStart = matrices[i].fineOffset;
                int inputStart = matrices[i].coarseOffset;
                int outputRows = matrices[i].M.rows();
                int inputRows = matrices[i].M.cols();

                Eigen::Map<Eigen::VectorXd> in_slice(in.data() + inputStart, inputRows);
                Eigen::Map<Eigen::VectorXd> out_slice(out.data() + outputStart, outputRows);
                out_slice = matrices[i].M * in_slice;
            }
        }

        template<typename V, typename Dest>
        void prolongVerts3X(V &in, Dest &out) {
            for (size_t i = 0; i < matrices.size(); i++) {
                int outputStart = matrices[i].fineOffset * 3;
                int inputStart = matrices[i].coarseOffset * 3;
                int outputRows = matrices[i].M.rows();
                int inputRows = matrices[i].M.cols();

                for (int c = 0; c < 3; c++) {
                    Eigen::Map<Eigen::VectorXd, 0, Eigen::InnerStride<3>> in_slice(in.data() + inputStart + c, inputRows);
                    Eigen::Map<Eigen::VectorXd, 0, Eigen::InnerStride<3>> out_slice(out.data() + outputStart + c, outputRows);
                    out_slice = matrices[i].M * in_slice;
                }
            }
        }

        template<typename V, typename Dest>
        void prolongBarycenter3X(V &in, Dest &out) {
            // Assume the first 3V rows are occupied by vertex data
            int lowerBase = lowerSize * 3;
            int upperBase = upperSize * 3;

            // Then just copy the next three entries
            out(upperBase) = in(lowerBase);
            out(upperBase + 1) = in(lowerBase + 1);
            out(upperBase + 2) = in(lowerBase + 2);
        }

        template<typename V, typename Dest>
        void prolongEdgeConstraints(V &in, Dest &out) {
            // Assume the first 3V rows are occupied by vertex data
            int lowerBase = lowerSize * 3 + 3;
            int upperBase = upperSize * 3 + 3;

            for (size_t i = 0; i < edgeMatrices.size(); i++) {
                int outputStart = upperBase + edgeMatrices[i].fineOffset;
                int inputStart = lowerBase + edgeMatrices[i].coarseOffset;
                int outputRows = edgeMatrices[i].M.rows();
                int inputRows = edgeMatrices[i].M.cols();

                Eigen::Map<Eigen::VectorXd> in_slice(in.data() + inputStart, inputRows);
                Eigen::Map<Eigen::VectorXd> out_slice(out.data() + outputStart, outputRows);
                out_slice = edgeMatrices[i].M * in_slice;
            }
        }

        template<typename V, typename Dest>
        void restrictVertsTranspose1X(V &in, Dest &out) {
            for (size_t i = 0; i < matrices.size(); i++) {
                int outputStart = matrices[i].coarseOffset;
                int inputStart = matrices[i].fineOffset;
                int outputRows = matrices[i].M.cols();
                int inputRows = matrices[i].M.rows();

                Eigen::Map<Eigen::VectorXd> in_slice(in.data() + inputStart, inputRows);
                Eigen::Map<Eigen::VectorXd> out_slice(out.data() + outputStart, outputRows);
                out_slice = matrices[i].M.transpose() * in_slice;
            }
        }

        template<typename V, typename Dest>
        void restrictVertsTranspose3X(V &in, Dest &out) {
            for (size_t i = 0; i < matrices.size(); i++) {
                int outputStart = matrices[i].coarseOffset * 3;
                int inputStart = matrices[i].fineOffset * 3;
                int outputRows = matrices[i].M.cols();
                int inputRows = matrices[i].M.rows();

                for (int c = 0; c < 3; c++) {
                    Eigen::Map<Eigen::VectorXd, 0, Eigen::InnerStride<3>> in_slice(in.data() + inputStart + c, inputRows);
                    Eigen::Map<Eigen::VectorXd, 0, Eigen::InnerStride<3>> out_slice(out.data() + outputStart + c, outputRows);
                    out_slice = matrices[i].M.transpose() * in_slice;
                }
            }
        }

        template<typename V, typename Dest>
        void restrictBarycenter3X(V &in, Dest &out) {
            // Assume the first 3V rows are occupied by vertex data
            int lowerBase = lowerSize * 3;
            int upperBase = upperSize * 3;

            // Then just copy the next three entries
            out(lowerBase) = in(upperBase);
            out(lowerBase + 1) = in(upperBase + 1);
            out(lowerBase + 2) = in(upperBase + 2);
        }

        template<typename V, typename Dest>
        void restrictEdgeConstraintsTranspose(V &in, Dest &out) {
            // Assume the first 3V rows are occupied by vertex data
            int lowerBase = lowerSize * 3 + 3;
            int upperBase = upperSize * 3 + 3;

            for (size_t i = 0; i < edgeMatrices.size(); i++) {
                int outputStart = lowerBase + edgeMatrices[i].coarseOffset;
                int inputStart = upperBase + edgeMatrices[i].fineOffset;
                int outputRows = edgeMatrices[i].M.cols();
                int inputRows = edgeMatrices[i].M.rows();

                Eigen::Map<Eigen::VectorXd> in_slice(in.data() + inputStart, inputRows);
                Eigen::Map<Eigen::VectorXd> out_slice(out.data() + outputStart, outputRows);
                out_slice = edgeMatrices[i].M.transpose() * in_slice;
            }
        }

        template<typename V, typename Dest>
        void restrictVertsPinv1X(V &in, Dest &out) {
            for (size_t i = 0; i < matrices.size(); i++) {
                int outputStart = matrices[i].coarseOffset;
                int inputStart = matrices[i].fineOffset;
                int outputRows = matrices[i].M.cols();
                int inputRows = matrices[i].M.rows();

                Eigen::VectorXd in_slice = Eigen::Map<Eigen::VectorXd>(in.data() + inputStart, inputRows);
                Eigen::Map<Eigen::VectorXd> out_slice(out.data() + outputStart, outputRows);
                out_slice = ApplyPinv(matrices[i].M, in_slice);
            }
        }

        template<typename V, typename Dest>
        void restrictVertsPinv3X(V &in, Dest &out) {
            for (size_t i = 0; i < matrices.size(); i++) {
                int outputStart = matrices[i].coarseOffset * 3;
                int inputStart = matrices[i].fineOffset * 3;
                int outputRows = matrices[i].M.cols();
                int inputRows = matrices[i].M.rows();

                for (int c = 0; c < 3; c++) {
                    Eigen::VectorXd in_slice = Eigen::Map<Eigen::VectorXd, 0, Eigen::InnerStride<3>>(in.data() + inputStart + c, inputRows);
                    Eigen::Map<Eigen::VectorXd, 0, Eigen::InnerStride<3>> out_slice(out.data() + outputStart + c, outputRows);
                    out_slice = ApplyPinv(matrices[i].M, in_slice);
                }
            }
        }

        template<typename V, typename Dest>
        void restrictEdgeConstraintsPinv(V &in, Dest &out) {
            // Assume the first 3V rows are occupied by vertex data
            int lowerBase = lowerSize * 3 + 3;
            int upperBase = upperSize * 3 + 3;

            for (size_t i = 0; i < edgeMatrices.size(); i++) {
                int outputStart = lowerBase + edgeMatrices[i].coarseOffset;
                int inputStart = upperBase + edgeMatrices[i].fineOffset;
                int outputRows = edgeMatrices[i].M.cols();
                int inputRows = edgeMatrices[i].M.rows();

                Eigen::VectorXd in_slice = Eigen::Map<Eigen::VectorXd>(in.data() + inputStart, inputRows);
                Eigen::Map<Eigen::VectorXd> out_slice(out.data() + outputStart, outputRows);
                out_slice = ApplyPinv(edgeMatrices[i].M, in_slice);
            }
        }
    };
}