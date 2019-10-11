#pragma once

#include <Eigen/Core>
#include <Eigen/Sparse>

namespace LWS {
    class TestMatrices {
        public:
        static Eigen::MatrixXd LaplacianDirichlet1D(int rows);
        static Eigen::MatrixXd LaplacianNeumann1D(int rows);
        static Eigen::MatrixXd LaplacianSaddle1D(int rows);
        static Eigen::MatrixXd CurveLaplacian(int rows);
        
        // static Eigen::MatrixXd CurveMetricLaplacian(PolyCurveGroup* curves, double epsilon);
    };
}
