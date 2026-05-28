#include "obstacles/sphere_obstacle.h"
#include "tpe_energy_sc.h"

namespace LWS {
SphereObstacle::SphereObstacle(const Vector3& c, double r, double p_exp) {
    center = c;
    radius = r;
    p = p_exp;
}

SphereObstacle::~SphereObstacle() {
}

double SphereObstacle::ComputeEnergy(const PolyCurveNetwork* curves) const {
    const int nVerts = curves->NumVertices();
    double sumE = 0;

    for (int i = 0; i < nVerts; i++) {
        const CurveVertex* v_i = curves->GetVertex(i);
        const Vector3 p_i = v_i->Position();
        // If we're very close to the center of the sphere, gradient is 0
        if ((p_i - center).norm() < 1e-6) continue;
        // Find the closest point on the plane
        const Vector3 nearest = ClosestPoint(p_i);
        // Simulate an energy contribution of 1 / r^(b - a)
        //Vector3 toPoint = nearest - p_i;
        const double distance = (nearest - p_i).norm();
        sumE += 1.0 / pow(distance, p);
    }

    return sumE;
}

void SphereObstacle::AddGradient(const PolyCurveNetwork* curves, Eigen::MatrixXd& gradient) const {
    const int nVerts = curves->NumVertices();

    for (int i = 0; i < nVerts; i++) {
        const CurveVertex* v_i = curves->GetVertex(i);
        const Vector3 p_i = v_i->Position();
        // If we're very close to the center of the sphere, gradient is 0
        if ((p_i - center).norm() < 1e-6) continue;
        // Find the closest point on the plane
        const Vector3 nearest = ClosestPoint(p_i);
        // Simulate an energy contribution of 1 / r^(b - a)
        Vector3 toPoint = nearest - p_i;
        const double distance = (nearest - p_i).norm();
        toPoint /= distance;
        const Vector3 grad_i = toPoint * p / pow(distance, p + 1);

        // Add an antipodal component
        // double oppDistance = 2 * radius - distance;
        // grad_i += -toPoint * 1.0 / pow(oppDistance, p);

        AddToRow(gradient, v_i->GlobalIndex(), grad_i);
    }
}
}
