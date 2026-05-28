#pragma once

#include "poly_curve_network.h"

namespace LWS {
class Obstacle {
public:
    virtual ~Obstacle() = 0;
    virtual void AddGradient(const PolyCurveNetwork* curves, Eigen::MatrixXd& gradient) const = 0;
    virtual double ComputeEnergy(const PolyCurveNetwork* curves) const = 0;
};
}