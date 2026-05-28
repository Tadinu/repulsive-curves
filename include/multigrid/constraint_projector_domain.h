#pragma once

#include "libgmultigrid/multigrid_domain.h"
#include "spatial/tpe_bvh.h"
#include "libgmultigrid/matrix_free.h"
#include "product/dense_matrix.h"
#include "constraint_projector_operator.h"

namespace LWS {
template <typename Constraint>
class ConstraintProjectorDomain : public MultigridDomain<BlockClusterTree, MatrixProjectorOperator> {
public:
    PolyCurveNetwork* curves = nullptr;
    BVHNode3D* bvh = nullptr;
    BlockClusterTree* tree = nullptr;
    double alpha = 0., beta = 0.;
    double sepCoeff = 0.;
    int nVerts = 0;
    double epsilon = 0.;
    Constraint constraint;
    bool isTopLevel = false;

    ConstraintProjectorDomain<Constraint>(PolyCurveNetwork* c, double a, double b, double sep,
                                          double diagEps = 0)
        : constraint(c) {
        curves = c;
        alpha = a;
        beta = b;
        sepCoeff = sep;
        nVerts = curves->NumVertices();
        epsilon = diagEps;

        bvh = CreateEdgeBVHFromCurve(curves);
        tree = new BlockClusterTree(curves, bvh, sepCoeff, alpha, beta, epsilon);
        tree->SetBlockTreeMode(BlockTreeMode::Matrix3AndProjector);

        curves->AddConstraintProjector(constraint);
        // std::cout << "Made level with " << nVerts << std::endl;
        isTopLevel = true;
    }

    ~ConstraintProjectorDomain<Constraint>() override {
        if (!isTopLevel) {
            delete curves;
        }
        delete tree;
        delete bvh;
    }

    MultigridDomain<BlockClusterTree, MatrixProjectorOperator>* Coarsen(MatrixProjectorOperator* prolongOp)
    const override {
        PolyCurveNetwork* coarsened = curves->Coarsen(prolongOp);
        ConstraintProjectorDomain<Constraint>* coarseDomain = new ConstraintProjectorDomain<Constraint>(
            coarsened, alpha, beta, sepCoeff, epsilon);
        prolongOp->lowerP = coarseDomain->GetConstraintProjector();
        prolongOp->upperP = GetConstraintProjector();
        coarseDomain->isTopLevel = false;
        return coarseDomain;
    }

    BlockClusterTree* GetMultiplier() const override {
        return tree;
    }

    Eigen::MatrixXd GetFullMatrix() const override {
        Eigen::MatrixXd A;
        SobolevCurves::Sobolev3XWithConstraints<Constraint>(curves, alpha, beta, A);
        return A;
    }

    Eigen::VectorXd DirectSolve(Eigen::VectorXd& b) const override {
        const int fullRows = constraint.SaddleNumRows();
        const int constraintRows = constraint.NumConstraintRows();

        Eigen::VectorXd b_aug(fullRows);
        b_aug.block(0, 0, b.rows(), 1) = b;
        b_aug.block(b.rows(), 0, constraintRows, 1).setZero();

        const Eigen::MatrixXd A = GetFullMatrix();
        b_aug = A.partialPivLu().solve(b_aug);
        return b_aug.block(0, 0, b.rows(), 1);
    }

    int NumVertices() const override {
        return nVerts;
    }

    int NumRows() const override {
        return 3 * nVerts;
    }

    MatrixProjectorOperator* MakeNewOperator() const override {
        return new MatrixProjectorOperator();
    }

    NullSpaceProjector* GetConstraintProjector() const {
        return curves->constraintProjector;
    }

    void UpdateCurvePositions() {
        bvh->recomputeCentersOfMass(curves);
        tree->refreshEdgeWeights();
    }
};
}