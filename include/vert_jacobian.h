#pragma once

#include "geometrycentral/utilities/vector3.h"
#include <Eigen/Core>

namespace LWS {
using namespace geometrycentral;

struct VertJacobian {
    Vector3 directional_x;
    Vector3 directional_y;
    Vector3 directional_z;

    // Multiply by a column vector on the right-hand side
    Vector3 RightMultiply(const Vector3& v) const;
    // Multiply by a row vector on the left-hand side
    Vector3 LeftMultiply(const Vector3& v) const;
    void Print() const;
    double Norm() const;
    void FillMatrix(Eigen::Matrix3d& matrix) const;

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
VertJacobian outer_product_to_jacobian(const Vector3& v1, const Vector3& v2);
}
