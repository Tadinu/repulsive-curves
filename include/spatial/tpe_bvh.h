#pragma once

#include "spatial_tree.h"
#include "../tpe_energy_sc.h"
#include "geometrycentral/utilities/vector2.h"
#include "geometrycentral/surface/halfedge_mesh.h"
#include "geometrycentral/surface/vertex_position_geometry.h"
#include "utils.h"

#include <fstream>
#include <Eigen/Core>

namespace LWS {
struct BHPlotData {
    double theta = 0.;
    double error = 0.;
    double gradientNorm = 0.;
    double minWidth = 0.;
    double maxWidth = 0.;
};

class BVHNode3D : public SpatialTree {
public:
    static int globalID;
    int thisNodeID = 0;
    int numNodes = 0;

    inline void recursivelyAssignIDs() {
        thisNodeID = globalID++;
        for (BVHNode3D* child : children) {
            child->recursivelyAssignIDs();
        }
    }

    inline void assignIDs() {
        globalID = 1;
        recursivelyAssignIDs();
    }

    inline void printIDs(std::ofstream& stream, int parentID = 0) const {
        if (parentID > 0) {
            stream << thisNodeID << ", " << parentID << std::endl;
        }
        for (BVHNode3D* child : children) {
            child->printIDs(stream, thisNodeID);
        }
    }

    // Build a BVH of the given points
    BVHNode3D(const std::vector<VertexBody6D>& points, int axis, BVHNode3D* root, bool splitTangents);
    virtual ~BVHNode3D() override;

    double totalMass = 0.;
    Vector3 centerOfMass;
    Vector3 averageTangent;
    std::vector<int> clusterIndices;
    BVHNode3D* bvhRoot = nullptr;
    Eigen::VectorXd fullMasses;

    // Fields for use by matrix-vector products; not used
    // by any BVH functions.
    double V_I = 0.;
    double B_I = 0.;
    double aIJ_VJ = 0.;

    inline void zeroMVFields() {
        V_I = 0;
        B_I = 0;
        aIJ_VJ = 0;
    }

    inline void recursivelyZeroMVFields() {
        zeroMVFields();
        if (!isLeaf) {
            for (BVHNode3D* child : children) {
                child->recursivelyZeroMVFields();
            }
        }
    }

    void findCurveSegments(const std::vector<VertexBody6D>& points, PolyCurveNetwork* curves);

    // Copy the new weights from the curves
    void refreshWeightsVector(const PolyCurveNetwork* curves, BodyType bType);

    // Recursively recompute all centers of mass in this tree
    template <typename T>
    void recomputeCentersOfMass(const T* curves);

    // Compute the total energy contribution from a single vertex
    virtual void accumulateVertexEnergy(double& result, const CurveVertex* i_pt,
                                        const PolyCurveNetwork* curves,
                                        double alpha, double beta) override;
    virtual void accumulateTPEGradient(Eigen::MatrixXd& gradients, const CurveVertex* i_pt,
                                       const PolyCurveNetwork* curves, double alpha, double beta) override;
    int NumElements() const;

    virtual double bodyEnergyEvaluation(const CurveVertex* i_pt, double alpha, double beta) const;
    virtual Vector3 bodyForceEvaluation(const CurveVertex* i_pt, double alpha, double beta) const;

    Vector3 exactGradient(const CurveVertex* basePoint, const PolyCurveNetwork* curves, double alpha,
                          double beta) const;

    PosTan minBound() const;
    PosTan maxBound() const;
    Vector3 BoxCenter() const;
    std::vector<BVHNode3D*> children;

    // for visualization, assign each node a range of indices
    // determined by the (post)ordering of its leaves
    int indexStart = 0;
    int indexEnd = 0;

    int indexNodes(int index) {
        indexStart = indexEnd = index;
        if (isLeaf) {
            indexEnd = index + 1;
            return index + 1;
        } else {
            for (BVHNode3D* child : children) {
                index = child->indexNodes(index);
            }
            indexEnd = index;
            return index;
        }
    }

    void accumulateChildren(std::vector<VertexBody6D>& result);
    Vector2 viewspaceBounds(const Vector3& point) const;

    inline void fillClusterMassVector(Eigen::VectorXd& w) const {
        w.setZero(clusterIndices.size());
        for (size_t i = 0; i < clusterIndices.size(); i++) {
            w(i) = bvhRoot->fullMasses(clusterIndices[i]);
        }
    }

    inline bool IsLeaf() const {
        return isLeaf;
    }

    inline bool IsEmpty() const {
        return isEmpty;
    }

    inline int VertexIndex() const {
        return body.elementIndex;
    }

    inline int countBadNodes() const {
        if (isEmpty) return 0;
        else if (isLeaf) {
            if (testTangent()) return 0;
            else return 1;
        } else {
            int sum = 0;
            if (!testTangent()) sum++;
            for (BVHNode3D* child : children) {
                sum += child->countBadNodes();
            }
            return sum;
        }
    }

    inline double nodeRatio(double d) const {
        // Compute diagonal distance from corner to corner
        // double diag = norm(maxCoords.position - minCoords.position);
        const Vector3 diag = maxCoords.position - minCoords.position;
        double maxCoord = fmax(diag.x, fmax(diag.y, diag.z));
        // double spatialR = diag.norm() / 2;
        return diag.norm() / d;

        // Vector3 tanDiag = maxCoords.tangent - minCoords.tangent;
        // double maxTanCoord = fmax(tanDiag.x, fmax(tanDiag.y, tanDiag.z));
        // double tangentR = tanDiag.norm() / 2;
        // return fmax(spatialR / d, tangentR);
    }

    bool shouldUseCell(const Vector3& vertPos) const;

    VertexBody6D body;

    //private:
    int numElements = 0;
    static double AxisSplittingPlane(const std::vector<VertexBody6D>& points, int axis);

    inline bool testTangent() const {
        const Vector3 tanDiag = maxCoords.tangent - minCoords.tangent;
        const double r = tanDiag.norm() / 2;
        return r < thresholdTheta;
    }

    template <typename T>
    void setLeafData(const T* curves);

    int splitAxis = 0;
    double splitPoint = 0.;
    bool isEmpty = false;
    bool isLeaf = false;
    PosTan minCoords;
    PosTan maxCoords;
    double thresholdTheta = 0.;
};

template <typename T>
void BVHNode3D::setLeafData(const T* curves) {
    std::cerr << "Type not supported" << std::endl;
    throw 1;
}

template <>
inline void BVHNode3D::setLeafData(const PolyCurveNetwork* curves) {
    if (body.type == BodyType::Vertex) {
        const CurveVertex* p = curves->GetVertex(body.elementIndex);
        body.mass = p->DualLength();
        body.pt.position = p->Position();
        body.pt.tangent = p->Tangent();
    } else if (body.type == BodyType::Edge) {
        const CurveEdge* p1 = curves->GetEdge(body.elementIndex);

        // Mass of an edge is its length
        body.mass = p1->Length();
        // Use midpoint as center of mass
        body.pt.position = p1->Midpoint();
        // Tangent direction is normalized edge vector
        body.pt.tangent = p1->Tangent();
    }

    totalMass = body.mass;
    centerOfMass = body.pt.position;
    averageTangent = body.pt.tangent;

    minCoords = PosTan{body.pt.position, body.pt.tangent};
    maxCoords = minCoords;
}

template <>
inline void BVHNode3D::setLeafData(
    const std::pair<std::shared_ptr<geometrycentral::surface::HalfedgeMesh>,
                    std::shared_ptr<geometrycentral::surface::VertexPositionGeometry>>* pair) {
    const auto& mesh = pair->first;
    const auto& geom = pair->second;

    using namespace geometrycentral;
    using namespace surface;

    if (body.type == BodyType::Vertex) {
        const Vertex v = mesh->vertex(body.elementIndex);
        body.mass = geom->vertexDualAreas[v];
        body.pt.position = geom->vertexPositions[v];
        body.pt.tangent = geom->vertexNormals[v];
    } else {
        std::cerr << "Element types besides vertex are not supported for meshes" << std::endl;
        throw 1;
    }

    totalMass = body.mass;
    centerOfMass = body.pt.position;
    averageTangent = body.pt.tangent;

    minCoords = PosTan{body.pt.position, body.pt.tangent};
    maxCoords = minCoords;
}

template <typename T>
inline void BVHNode3D::recomputeCentersOfMass(const T* curves) {
    if (isEmpty) {
        totalMass = 0;
        numElements = 0;
    }
    // For a leaf, just set centers and bounds from the one body
    else if (isLeaf) {
        setLeafData(curves);
        numElements = 1;
    } else {
        // Recursively compute bounds for all children
        for (auto* child : children) {
            child->recomputeCentersOfMass(curves);
        }

        minCoords = children[0]->minCoords;
        maxCoords = children[0]->maxCoords;

        totalMass = 0;
        centerOfMass = Vector3::zero();
        averageTangent = Vector3::zero();

        // Accumulate max/min over all nonempty children
        numElements = 0;
        for (const auto* child : children) {
            if (!child->isEmpty) {
                minCoords = postan_min(child->minCoords, minCoords);
                maxCoords = postan_max(child->maxCoords, maxCoords);

                totalMass += child->totalMass;
                centerOfMass += child->centerOfMass * child->totalMass;
                averageTangent += child->averageTangent * child->totalMass;
            }
            numElements += child->numElements;
        }

        centerOfMass /= totalMass;
        averageTangent /= totalMass;

        averageTangent = averageTangent.normalize();
    }
}

BVHNode3D* CreateBVHFromCurve(const PolyCurveNetwork* curves);
BVHNode3D* CreateEdgeBVHFromCurve(const PolyCurveNetwork* curves);
BVHNode3D* CreateBVHFromMesh(const std::shared_ptr<geometrycentral::surface::HalfedgeMesh>& mesh,
                             const std::shared_ptr<geometrycentral::surface::VertexPositionGeometry>& geom);
}
