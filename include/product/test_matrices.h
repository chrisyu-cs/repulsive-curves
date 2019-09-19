#pragma once

#include <Eigen/Core>

namespace LWS {
    class TestMatrices {
        public:
        static Eigen::MatrixXd LaplacianDirichlet1D(int rows);
        static Eigen::MatrixXd LaplacianNeumann1D(int rows);
        static Eigen::MatrixXd LaplacianSaddle1D(int rows);
        static Eigen::MatrixXd CurveLaplacian(int rows);
    };
}
