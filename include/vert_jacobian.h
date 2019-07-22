#pragma once

#include "boundary_derivatives.h"
#include "boundary_loop.h"

namespace LWS {
    
    struct VertJacobian {
        Vector3 directional_x;
        Vector3 directional_y;
        Vector3 directional_z;

        // Multiply by a column vector on the right-hand side
        Vector3 RightMultiply(Vector3 v);
        // Multiply by a row vector on the left-hand side
        Vector3 LeftMultiply(Vector3 v);
        void Print();
        double Norm();

        // Add two matrices
        friend VertJacobian operator+(const VertJacobian& a, const VertJacobian& b);
        // Subtract two matrices
        friend VertJacobian operator-(const VertJacobian& a, const VertJacobian& b);
        // Multiply by a scalar
        friend VertJacobian operator*(const VertJacobian& a, double c);
        friend VertJacobian operator*(double c, const VertJacobian& a);
    };

    VertJacobian operator+(const VertJacobian& a, const VertJacobian& b);
    VertJacobian operator-(const VertJacobian& a, const VertJacobian& b);
    VertJacobian operator*(const VertJacobian& a, double c);
    VertJacobian operator*(double c, const VertJacobian& a);
    VertJacobian outer_product_to_jacobian(Vector3 v1, Vector3 v2);

}
