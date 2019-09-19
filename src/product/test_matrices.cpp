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

        return M;
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
                M(i, i) = 2.001;
                M(i, i + 1) = -1;
            }
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

    
}