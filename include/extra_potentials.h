#pragma once

#include <Eigen/Core>
#include "poly_curve_network.h"
#include "vert_jacobian.h"

namespace LWS {
class CurvePotential {
public:
    CurvePotential();
    virtual ~CurvePotential();
    virtual double CurrentValue(const PolyCurveNetwork* curves) const;
    virtual void AddGradient(const PolyCurveNetwork* curves, Eigen::MatrixXd& gradient) const;
};

class TotalLengthPotential : public CurvePotential {
public:
    explicit TotalLengthPotential(double wt);
    virtual double CurrentValue(const PolyCurveNetwork* curves) const override;
    virtual void AddGradient(const PolyCurveNetwork* curves, Eigen::MatrixXd& gradient) const override;

private:
    double weight;
};

class LengthDifferencePotential : public CurvePotential {
public:
    explicit LengthDifferencePotential(double wt);
    virtual double CurrentValue(const PolyCurveNetwork* curves) const override;
    virtual void AddGradient(const PolyCurveNetwork* curves, Eigen::MatrixXd& gradient) const override;

private:
    double LenDiff(const PolyCurveNetwork* curves, int i) const;
    double weight;
};

class PinBendingPotential : public CurvePotential {
public:
    explicit PinBendingPotential(double wt);
    virtual double CurrentValue(const PolyCurveNetwork* curves) const override;
    virtual void AddGradient(const PolyCurveNetwork* curves, Eigen::MatrixXd& gradient) const override;

private:
    double weight;
};

class VectorField {
public:
    VectorField();
    virtual ~VectorField();
    virtual Vector3 Sample(const Vector3& x) const;
    virtual VertJacobian SpatialDerivative(const Vector3& x) const;
};

class ConstantVectorField : public VectorField {
public:
    explicit ConstantVectorField(const Vector3& v);
    virtual Vector3 Sample(const Vector3& x) const override;
    virtual VertJacobian SpatialDerivative(const Vector3& x) const override;

private:
    Vector3 c;
};

class CircularVectorField : public VectorField {
public:
    CircularVectorField();
    virtual Vector3 Sample(const Vector3& x) const override;
    virtual VertJacobian SpatialDerivative(const Vector3& x) const override;

private:
    Vector3 dirDeriv(const Vector3& x, const Vector3& dir) const;
};

class InterestingVectorField : public VectorField {
public:
    InterestingVectorField();
    virtual Vector3 Sample(const Vector3& x) const override;
    virtual VertJacobian SpatialDerivative(const Vector3& x) const override;
};

class VectorFieldPotential : public CurvePotential {
public:
    VectorFieldPotential(double wt, VectorField* vf);
    ~VectorFieldPotential() override;
    virtual double CurrentValue(const PolyCurveNetwork* curves) const override;
    virtual void AddGradient(const PolyCurveNetwork* curves, Eigen::MatrixXd& gradient) const override;

private:
    double weight = 0.;
    VectorField* field = nullptr;
};

// class AreaPotential : public CurvePotential {
//     public:
//     AreaPotential(double wt);
//     virtual double CurrentValue(const PolyCurveNetwork* curves)const override;
//     virtual void AddGradient(const PolyCurveNetwork* curves, Eigen::MatrixXd &gradient)const override;

//     private:
//     double weight= 0.;
// };
}