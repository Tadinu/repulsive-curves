#pragma once

#include "geometrycentral/utilities/vector3.h"

namespace LWS {
using namespace geometrycentral;

enum class BodyType {
    Vertex, Edge
};

struct VertexBody {
    Vector3 position;
    Vector3 tangent;
    double mass;
    int vertIndex;
};

struct PosTan {
    Vector3 position;
    Vector3 tangent;
    void Print() const;
};

struct VertexBody6D {
    PosTan pt;
    double mass;
    int elementIndex;
    BodyType type;
};

PosTan postan_max(const PosTan& v1, const PosTan& v2);
PosTan postan_min(const PosTan& v1, const PosTan& v2);

Vector3 CartesianToSpherical(const Vector3& cartesian);

Vector3 SphericalToCartesian(const Vector3& spherical);
}