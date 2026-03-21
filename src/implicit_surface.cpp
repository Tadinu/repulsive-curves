#include "implicit_surface.h"
#include "utils.h"

#include "geometrycentral/utilities/vector2.h"

namespace LWS {
ImplicitSurface::~ImplicitSurface() {
}

ImplicitSphere::ImplicitSphere(double r, const Vector3& c) {
    radius = r;
    center = c;
}

double ImplicitSphere::SignedDistance(const Vector3& point) const {
    return (point - center).norm() - radius;
}

Vector3 ImplicitSphere::GradientOfDistance(const Vector3& point) const {
    return (point - center).normalize();
}

double ImplicitSphere::BoundingDiameter() const {
    return 2 * radius;
}

Vector3 ImplicitSphere::BoundingCenter() const {
    return center;
}

ImplicitTorus::ImplicitTorus(double major, double minor, const Vector3& c) {
    majorRadius = major;
    minorRadius = minor;
    center = c;
}

double ImplicitTorus::SignedDistance(const Vector3& point) const {
    // Vector3 disp = point - center;
    // Vector2 planar{disp.x, disp.z};
    // Vector2 q{norm(planar) - majorRadius, disp.y};
    // return norm(q) - minorRadius;

    const Vector3 p = point - center;
    const double normAll = norm2(p);
    const double normXZ = norm(Vector2{p.x, p.z});
    return sqrt(normAll - 2 * majorRadius * normXZ + majorRadius * majorRadius) - minorRadius;
}

Vector3 ImplicitTorus::GradientOfDistance(const Vector3& point) const {
    const Vector3 p = point - center;
    const double normAll = norm2(p);
    const double normXZ = norm(Vector2{p.x, p.z});
    const double sqTerm = sqrt(normAll - 2 * majorRadius * normXZ + majorRadius * majorRadius);

    const double partialX = 2 * p.x - (2 * p.x * majorRadius) / normXZ;
    const double partialY = 2 * p.y;
    const double partialZ = 2 * p.z - (2 * p.z * majorRadius) / normXZ;

    return (1.0 / (2 * sqTerm)) * Vector3{partialX, partialY, partialZ};
}

double ImplicitTorus::BoundingDiameter() const {
    return majorRadius * 2 + minorRadius * 2;
}

Vector3 ImplicitTorus::BoundingCenter() const {
    return center;
}

YZeroPlane::YZeroPlane() {
}

double YZeroPlane::SignedDistance(const Vector3& point) const {
    return point.y;
}

Vector3 YZeroPlane::GradientOfDistance(const Vector3& point) const {
    return Vector3{0, 1, 0};
}

double YZeroPlane::BoundingDiameter() const {
    return 2;
}

Vector3 YZeroPlane::BoundingCenter() const {
    return Vector3{0, 0.01, 0};
}

ImplicitDoubleTorus::ImplicitDoubleTorus(double r2) {
    radius2 = r2;
}

double ImplicitDoubleTorus::SignedDistance(const Vector3& p) const {
    const double x = p.x + 2;
    const double y = p.y;
    const double z = p.z;
    const double inner = (x * (x - 1) * (x - 1) * (x - 2) + y * y);
    return inner * inner + z * z - radius2;
}

Vector3 ImplicitDoubleTorus::GradientOfDistance(const Vector3& point) const {
    const double x = point.x + 2;
    const double y = point.y;
    const double z = point.z;

    const double x2 = x * x;
    const double x3 = x2 * x;
    const double x4 = x3 * x;
    const double y2 = y * y;
    const double polyn = x4 - 4 * x3 + 5 * x2 - 2 * x + y2;

    const double partialX = 2 * polyn * (4 * x3 - 12 * x2 + 10 * x - 2);
    const double partialY = 4 * y * polyn;
    const double partialZ = 2 * z;

    return Vector3{partialX, partialY, partialZ};
}

double ImplicitDoubleTorus::BoundingDiameter() const {
    return 4.5;
}

Vector3 ImplicitDoubleTorus::BoundingCenter() const {
    return Vector3::zero();
}

ImplicitUnion::ImplicitUnion(ImplicitSurface* s1, ImplicitSurface* s2) {
    surfaces.push_back(s1);
    surfaces.push_back(s2);
}

double ImplicitUnion::SignedDistance(const Vector3& p) const {
    double minDist = surfaces[0]->SignedDistance(p);
    for (size_t i = 1; i < surfaces.size(); i++) {
        minDist = fmin(minDist, surfaces[i]->SignedDistance(p));
    }
    return minDist;
}

Vector3 ImplicitUnion::GradientOfDistance(const Vector3& point) const {
    double minDist = surfaces[0]->SignedDistance(point);
    double min_i = 0;

    for (size_t i = 1; i < surfaces.size(); i++) {
        const double dist = surfaces[i]->SignedDistance(point);
        if (dist < minDist) {
            minDist = dist;
            min_i = i;
        }
    }

    return surfaces[min_i]->GradientOfDistance(point);
}

double ImplicitUnion::BoundingDiameter() const {
    const Vector3 center = BoundingCenter();
    double maxRadius = 0;

    for (const auto* surface : surfaces) {
        const Vector3 c_i = surface->BoundingCenter();
        const double d_i = surface->BoundingDiameter();
        const Vector3 disp = (center - c_i);
        const double l1dist = fmax(fabs(disp.x), fmax(fabs(disp.y), fabs(disp.z)));

        maxRadius = fmax(maxRadius, l1dist + d_i / 2);
    }
    return maxRadius * 2;
}

Vector3 ImplicitUnion::BoundingCenter() const {
    Vector3 mins = surfaces[0]->BoundingCenter();
    Vector3 maxs = mins;

    for (const auto* surface : surfaces) {
        mins = vector_min(mins, surface->BoundingCenter());
        maxs = vector_max(maxs, surface->BoundingCenter());
    }

    return (mins + maxs) / 2;
}

ImplicitSmoothUnion::ImplicitSmoothUnion(ImplicitSurface* s1, ImplicitSurface* s2, double blendFactor) {
    surface1 = s1;
    surface2 = s2;
    k = blendFactor;
}

double ImplicitSmoothUnion::SignedDistance(const Vector3& p) const {
    const double d1 = surface1->SignedDistance(p);
    const double d2 = surface2->SignedDistance(p);
    const float h = clamp(0.5 + 0.5 * (d2 - d1) / k, 0.0, 1.0);
    return ((1 - h) * d2 + h * d1) - k * h * (1 - h);
}


Vector3 ImplicitSmoothUnion::GradientOfDistance(const Vector3& point) const {
    const double d1 = surface1->SignedDistance(point);
    const double d2 = surface2->SignedDistance(point);
    const float h = clamp(0.5 + 0.5 * (d2 - d1) / k, 0.0, 1.0);

    const Vector3 deriv_d1 = surface1->GradientOfDistance(point);
    const Vector3 deriv_d2 = surface2->GradientOfDistance(point);
    // We need derivative of a clamp, which is just 0 outside of (0, 1)
    Vector3 deriv_h = Vector3::zero();
    if (h > 0 && h < 1) {
        deriv_h = (1.0 / (2 * k)) * (deriv_d2 - deriv_d1);
    }

    const Vector3 deriv_term1 = (-deriv_h * d2 + (1 - h) * deriv_d2);
    const Vector3 deriv_term2 = (deriv_h * d1 + h * deriv_d1);
    const Vector3 deriv_term3 = -k * (deriv_h - 2 * h * deriv_h);

    return deriv_term1 + deriv_term2 + deriv_term3;
}

double ImplicitSmoothUnion::BoundingDiameter() const {
    const Vector3 center = BoundingCenter();
    double maxRadius = 0;

    const ImplicitSurface* surfaces[2] = {surface1, surface2};

    for (size_t i = 0; i < 2; i++) {
        const Vector3 c_i = surfaces[i]->BoundingCenter();
        const double d_i = surfaces[i]->BoundingDiameter();
        const Vector3 disp = (center - c_i);
        const double l1dist = fmax(fabs(disp.x), fmax(fabs(disp.y), fabs(disp.z)));

        maxRadius = fmax(maxRadius, l1dist + d_i / 2);
    }
    return maxRadius * 2;
}

Vector3 ImplicitSmoothUnion::BoundingCenter() const {
    Vector3 mins = surface1->BoundingCenter();
    Vector3 maxs = mins;
    mins = vector_min(mins, surface2->BoundingCenter());
    maxs = vector_max(maxs, surface2->BoundingCenter());

    return (mins + maxs) / 2;
}
}

