#pragma once

#include "libgmultigrid/multigrid_operator.h"

namespace LWS {

    class MatrixProjectorOperator : public MultigridOperator {
        public:

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

        void checkAndSetSize(int expectedInput, int expectedOutput, Eigen::VectorXd &v, Eigen::VectorXd &out) {
            out.setZero(expectedOutput);
            if (v.rows() != expectedInput) {
                std::cerr << "Input doesn't match expected size: " << v.rows() << " vs " << expectedInput << std::endl;
                exit(1);
            }
        }

        Eigen::VectorXd prolong(Eigen::VectorXd v) {
            Eigen::VectorXd out;
            int expectedInput = 3 * lowerSize;
            int expectedOutput = 3 * upperSize;
            checkAndSetSize(expectedInput, expectedOutput, v, out);

            v = lowerP->ProjectToNullspace(v);
            prolongVerts3X(v, out);
            out = upperP->ProjectToNullspace(out);
            return out;
        }

        Eigen::VectorXd restrictWithTranspose(Eigen::VectorXd v) {
            Eigen::VectorXd out;
            int expectedInput = 3 * upperSize;
            int expectedOutput = 3 * lowerSize;
            checkAndSetSize(expectedInput, expectedOutput, v, out);

            v = upperP->ProjectToNullspace(v);
            restrictVertsTranspose3X(v, out);
            out = lowerP->ProjectToNullspace(out);
            return out;
        }

        Eigen::VectorXd restrictWithPinv(Eigen::VectorXd v) {
            Eigen::VectorXd out;
            int expectedInput = 3 * upperSize;
            int expectedOutput = 3 * lowerSize;
            checkAndSetSize(expectedInput, expectedOutput, v, out);

            v = upperP->ProjectToNullspace(v);
            restrictVertsPinv3X(v, out);
            out = lowerP->ProjectToNullspace(out);
            return out;
        }
    };
}