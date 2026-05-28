#pragma once

#include "obstacles/obstacle.h"

namespace LWS {
class SphereObstacle : public Obstacle {
public:
    Vector3 center;
    double radius;
    double p;
    SphereObstacle(const Vector3& c, double r, double p_exp);
    virtual ~SphereObstacle() override;
    virtual void AddGradient(const PolyCurveNetwork* curves, Eigen::MatrixXd& gradient) const override;
    virtual double ComputeEnergy(const PolyCurveNetwork* curves) const override;

    inline Vector3 ClosestPoint(Vector3 input) const {
        const Vector3 dir = (input - center).normalize();
        return center + radius * dir;
    }
};
}

