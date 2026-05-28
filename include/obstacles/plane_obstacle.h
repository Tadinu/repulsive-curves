#pragma once

#include "obstacles/obstacle.h"

namespace LWS {
class PlaneObstacle : public Obstacle {
public:
    Vector3 center;
    Vector3 normal;
    double p;
    double weight;
    PlaneObstacle(const Vector3& c, const Vector3& n, double p_exp, double wt);
    virtual ~PlaneObstacle() override;
    virtual void AddGradient(const PolyCurveNetwork* curves, Eigen::MatrixXd& gradient) const override;
    virtual double ComputeEnergy(const PolyCurveNetwork* curves) const override;

    inline Vector3 ClosestPoint(const Vector3& input) const {
        const Vector3 fromInput = center - input;
        const Vector3 normalProj = normal * dot(normal, fromInput);
        return input + normalProj;
    }
};
}
