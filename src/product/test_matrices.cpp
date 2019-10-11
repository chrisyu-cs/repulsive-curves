#include "product/test_matrices.h"

#include <iostream>

namespace LWS {
    Eigen::MatrixXd TestMatrices::LaplacianDirichlet1D(int rows) {
        Eigen::MatrixXd M;
        M.setZero(rows, rows);

        double mass = 1.0 / rows;

        for (int i = 0; i < rows; i++) {
            if (i == 0) {
                M(i, i) = 2;
                M(i, i + 1) = -1;
            }
            else if (i == rows - 1) {
                M(i, i - 1) = -1;
                M(i, i) = 2;
            }
            else {
                M(i, i) = 2;
                M(i, i - 1) = -1;
                M(i, i + 1) = -1;
            }
        }

        return M * mass;
    }

    Eigen::MatrixXd TestMatrices::LaplacianNeumann1D(int rows) {
        Eigen::MatrixXd M;
        M.setZero(rows, rows);

        double mass = 1.0 / rows;

        for (int i = 0; i < rows; i++) {
            if (i == 0) {
                M(i, i) = 1;
                M(i, i + 1) = -1;
            }
            else if (i == rows - 1) {
                M(i, i) = 1;
                M(i, i - 1) = -1;
            }
            else {
                M(i, i - 1) = -1;
                M(i, i) = 2;
                M(i, i + 1) = -1;
            }
        }

        M *= mass;

        for (int i = 0; i < rows; i++) {
            M(i, i) += 1e-4 * mass;
        }

        return M;
    }

    Eigen::MatrixXd TestMatrices::LaplacianSaddle1D(int rows) {
        Eigen::MatrixXd M;
        M.setZero(rows + 1, rows + 1);

        double mass = 1.0 / rows;

        for (int i = 0; i < rows; i++) {
            if (i == 0) {
                M(i, i) = -1;
                M(i, i + 1) = 1;
            }
            else if (i == rows - 1) {
                M(i, i) = -1;
                M(i, i - 1) = 1;
            }
            else {
                M(i, i - 1) = 1;
                M(i, i) = -2;
                M(i, i + 1) = 1;
            }
        }

        for (int i = 0; i < rows; i++) {
            M(rows, i) = mass;
            M(i, rows) = mass;
        }

        return M;
    }

    Eigen::MatrixXd TestMatrices::CurveLaplacian(int rows) {
        Eigen::MatrixXd M;
        M.setZero(rows, rows);

        for (int i = 0; i < rows; i++) {
            M(i, (i - 1 + rows) % rows) = -1;
            M(i, i) = 2 + 1e-3;
            M(i, (i + 1) % rows) = -1;
        }

        return M;
    }

    // Eigen::MatrixXd TestMatrices::CurveMetricLaplacian(PolyCurveGroup* curves, double epsilon) {
    //     int rows = curves->NumVertices();

    //     std::vector<Eigen::Triplet<double>> triplets;
    //     Eigen::VectorXd masses(rows);

    //     for (int i = 0; i < rows; i++) {
    //         double mass = curves->GetCurvePoint(i).DualLength();
    //         triplets.push_back(Eigen::Triplet<double>(i, i, -1.0 / mass));
    //         triplets.push_back(Eigen::Triplet<double>(i, (i+1) % rows, 1.0 / mass));
    //         masses(i) = mass;
    //     }

    //     Eigen::SparseMatrix<double> G;
    //     G.resize(rows, rows);
    //     G.setFromTriplets(triplets.begin(), triplets.end());

    //     Eigen::SparseMatrix<double> laplacian = G.transpose() * masses.asDiagonal() * G;
    //     Eigen::MatrixXd L = laplacian.toDense();

    //     for (int i = 0; i < rows; i++) {
    //         L(i, i) += epsilon * masses(i);
    //     }

    //     return L;
    // }

    
}