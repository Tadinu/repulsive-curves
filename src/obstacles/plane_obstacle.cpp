#include "obstacles/plane_obstacle.h"
#include "tpe_energy_sc.h"

namespace LWS {
PlaneObstacle::PlaneObstacle(const Vector3& c, const Vector3& n, double p_exp, double wt) {
    center = c;
    normal = n.normalize();
    p = p_exp;
    weight = wt;
}

PlaneObstacle::~PlaneObstacle() {
}

double PlaneObstacle::ComputeEnergy(const PolyCurveNetwork* curves) const {
    const int nVerts = curves->NumVertices();
    double sumE = 0;

    for (int i = 0; i < nVerts; i++) {
        const CurveVertex* v_i = curves->GetVertex(i);
        const Vector3 p_i = v_i->Position();
        // Find the closest point on the plane
        const Vector3 nearest = ClosestPoint(p_i);
        // Simulate an energy contribution of 1 / r^p
        //const Vector3 toPoint = nearest - p_i;
        const double distance = (nearest - p_i).norm();
        const double energy = 1.0 / pow(distance, p);
        sumE += energy;
    }

    return weight * sumE;
}

void PlaneObstacle::AddGradient(const PolyCurveNetwork* curves, Eigen::MatrixXd& gradient) const {
    const int nVerts = curves->NumVertices();

    for (int i = 0; i < nVerts; i++) {
        const CurveVertex* v_i = curves->GetVertex(i);
        const Vector3 p_i = v_i->Position();
        // Find the closest point on the plane
        const Vector3 nearest = ClosestPoint(p_i);
        // Simulate an energy contribution of 1 / r^p
        Vector3 toPoint = nearest - p_i;
        const double distance = (nearest - p_i).norm();
        toPoint /= distance;
        const Vector3 grad_i = toPoint * p / pow(distance, p + 1);

        AddToRow(gradient, v_i->GlobalIndex(), weight * grad_i);
    }
}
}