#pragma once

#include "utils.h"
#include "geometrycentral/utilities/vector3.h"

namespace LWS {
using namespace geometrycentral;

class ImplicitSurface {
public:
    virtual ~ImplicitSurface();
    virtual double SignedDistance(const Vector3& point) const = 0;
    virtual Vector3 GradientOfDistance(const Vector3& point) const = 0;
    virtual double BoundingDiameter() const = 0;
    virtual Vector3 BoundingCenter() const = 0;
};

class ImplicitSphere : public ImplicitSurface {
public:
    ImplicitSphere(double r, const Vector3& c);
    virtual double SignedDistance(const Vector3& point) const override;
    virtual Vector3 GradientOfDistance(const Vector3& point) const override;
    virtual double BoundingDiameter() const override;
    virtual Vector3 BoundingCenter() const override;

private:
    double radius = 0.;
    Vector3 center = Vector3::zero();
};

class ImplicitTorus : public ImplicitSurface {
public:
    ImplicitTorus(double major, double minor, const Vector3& c);
    virtual double SignedDistance(const Vector3& point) const override;
    virtual Vector3 GradientOfDistance(const Vector3& point) const override;
    virtual double BoundingDiameter() const override;
    virtual Vector3 BoundingCenter() const override;

private:
    double majorRadius = 0.;
    double minorRadius = 0.;
    Vector3 center = Vector3::zero();
};

class YZeroPlane : public ImplicitSurface {
public:
    YZeroPlane();
    virtual double SignedDistance(const Vector3& point) const override;
    virtual Vector3 GradientOfDistance(const Vector3& point) const override;
    virtual double BoundingDiameter() const override;
    virtual Vector3 BoundingCenter() const override;
};

class ImplicitDoubleTorus : public ImplicitSurface {
public:
    ImplicitDoubleTorus(double r2);
    virtual double SignedDistance(const Vector3& point) const override;
    virtual Vector3 GradientOfDistance(const Vector3& point) const override;
    virtual double BoundingDiameter() const override;
    virtual Vector3 BoundingCenter() const override;

private:
    double radius2 = 0.;
};

class ImplicitUnion : public ImplicitSurface {
public:
    ImplicitUnion(ImplicitSurface* s1, ImplicitSurface* s2);
    virtual double SignedDistance(const Vector3& point) const override;
    virtual Vector3 GradientOfDistance(const Vector3& point) const override;
    virtual double BoundingDiameter() const override;
    virtual Vector3 BoundingCenter() const override;

private:
    std::vector<ImplicitSurface*> surfaces;
};

class ImplicitSmoothUnion : public ImplicitSurface {
public:
    ImplicitSmoothUnion(ImplicitSurface* s1, ImplicitSurface* s2, double blendFactor);
    virtual double SignedDistance(const Vector3& point) const override;
    virtual Vector3 GradientOfDistance(const Vector3& point) const override;
    virtual double BoundingDiameter() const override;
    virtual Vector3 BoundingCenter() const override;

private:
    ImplicitSurface* surface1 = nullptr;
    ImplicitSurface* surface2 = nullptr;
    double k = 0.;
};
}