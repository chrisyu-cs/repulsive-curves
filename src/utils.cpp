#include "utils.h"
#include <chrono>

namespace LWS {

    Vector3 vector_max(Vector3 a, Vector3 b) {
        return Vector3{fmax(a.x, b.x), fmax(a.y, b.y), fmax(a.z, b.z)};
    }

    Vector3 vector_min(Vector3 a, Vector3 b) {
        return Vector3{fmin(a.x, b.x), fmin(a.y, b.y), fmin(a.z, b.z)};
    }

    Vector3 vector_abs(Vector3 a) {
        return Vector3{fabs(a.x), fabs(a.y), fabs(a.z)};
    }

    Vector3 heVector(surface::VertexPositionGeometry* geom, surface::Halfedge he) {
        return geom->vertexPositions[he.twin().vertex()] - geom->vertexPositions[he.vertex()];
    }

    long Utils::currentTimeMilliseconds() {
        using namespace std::chrono;
        return duration_cast<milliseconds>(system_clock::now().time_since_epoch()).count();
    }

    double std_vector_dot(std::vector<double> &x, std::vector<double> &y) {
        int s = x.size();
        double sum = 0;
        for (int i = 0; i < s; i++) {
            sum += x[i] * y[i];
        }
        return sum;
    }

    double std_vector_sum_entries(std::vector<double> &x) {
        int s = x.size();
        double sum = 0;
        for (int i = 0; i < s; i++) {
            sum += x[i];
        }
        return sum;
    }

    void std_vector_add(std::vector<double> &x, std::vector<double> &y, std::vector<double> &result) {
        int s = x.size();
        for (int i = 0; i < s; i++) {
            result[i] = x[i] + y[i];
        }
    }

    Eigen::VectorXd VectorToVectorXd(std::vector<double> &x) {
        Eigen::VectorXd v;
        v.setZero(x.size());
        for (size_t i = 0; i < x.size(); i++) {
            v(i) = x[i];
        }
        return v;
    }

    Eigen::MatrixXd Vector3ToMatrixXd(std::vector<Vector3> &x) {
        Eigen::MatrixXd v;
        v.setZero(x.size(), 3);
        for (size_t i = 0; i < x.size(); i++) {
            v(i, 0) = x[i].x;
            v(i, 1) = x[i].y;
            v(i, 2) = x[i].z;
        }
        return v;
    }

}
