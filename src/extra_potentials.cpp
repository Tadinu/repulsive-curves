#include "extra_potentials.h"

#include "tpe_energy_sc.h"

namespace LWS {
// Base class doesn't do anything
CurvePotential::CurvePotential() {
}

CurvePotential::~CurvePotential() {
}

double CurvePotential::CurrentValue(const PolyCurveNetwork* curves) const {
    return 0;
}

void CurvePotential::AddGradient(const PolyCurveNetwork* curves, Eigen::MatrixXd& gradient) const {
}

TotalLengthPotential::TotalLengthPotential(double wt) {
    weight = wt;
}

double TotalLengthPotential::CurrentValue(const PolyCurveNetwork* curves) const {
    double len = curves->TotalLength();
    return weight * 0.5 * len * len;
}

void TotalLengthPotential::AddGradient(const PolyCurveNetwork* curves, Eigen::MatrixXd& gradient) const {
    int nVerts = curves->NumVertices();
    double len = curves->TotalLength();

    Eigen::MatrixXd lenGradient;
    lenGradient.setZero(gradient.rows(), gradient.cols());

    for (int i = 0; i < nVerts; i++) {
        // Every vertex length needs to be differentiated by its neighbors
        CurveVertex* v_i = curves->GetVertex(i);
        int nNeighbors = v_i->numEdges();
        for (int j = 0; j < nNeighbors; j++) {
            CurveVertex* v_neighbor = curves->GetVertex(j);
            Vector3 lenGrad = TPESC::length_wrt_vert(v_i, v_neighbor);
            AddToRow(lenGradient, v_neighbor->GlobalIndex(), weight * lenGrad);
        }
        // Length also needs to be differentiated by itself
        Vector3 selfGrad = TPESC::length_wrt_vert(v_i, v_i);
        AddToRow(lenGradient, v_i->GlobalIndex(), weight * selfGrad);
    }

    gradient += (weight * len) * lenGradient;
}

LengthDifferencePotential::LengthDifferencePotential(double wt) {
    weight = wt;
}

double LengthDifferencePotential::CurrentValue(const PolyCurveNetwork* curves) const {
    int nVerts = curves->NumVertices();
    double diffs2 = 0;
    for (int i = 0; i < nVerts; i++) {
        diffs2 += LenDiff(curves, i);
    }
    return 0.5 * weight * diffs2;
}

double LengthDifferencePotential::LenDiff(const PolyCurveNetwork* curves, int i) const {
    CurveVertex* v_i = curves->GetVertex(i);
    if (v_i->numEdges() != 2) {
        return 0;
    }
    CurveEdge* e1 = v_i->edge(0);
    CurveEdge* e2 = v_i->edge(1);
    double lenDiff = e2->Length() - e1->Length();

    return lenDiff * lenDiff;
}

void LengthDifferencePotential::AddGradient(const PolyCurveNetwork* curves, Eigen::MatrixXd& gradient) const {
    int nVerts = curves->NumVertices();
    double len = curves->TotalLength();

    for (int i = 0; i < nVerts; i++) {
        // Each length difference centered at i affects vertices (i-1, i, i+1).
        const double lDiff = LenDiff(curves, i);
        if (lDiff == 0) continue;

        const CurveVertex* v_i = curves->GetVertex(i);
        const CurveEdge* e1 = v_i->edge(0);
        const CurveEdge* e2 = v_i->edge(1);
        const CurveVertex* v_prev = e1->Opposite(v_i);
        const CurveVertex* v_next = e2->Opposite(v_i);

        // Differentiate (length(e2) - length(e1)) by the 3 vertices
        const Vector3 deriv2_next = TPESC::edge_length_wrt_vert(e2, v_next);
        const Vector3 deriv2_mid = TPESC::edge_length_wrt_vert(e2, v_i);
        const Vector3 deriv1_mid = -TPESC::edge_length_wrt_vert(e1, v_i);
        const Vector3 deriv1_prev = -TPESC::edge_length_wrt_vert(e1, v_prev);

        const double lenWeight = lDiff * weight;

        AddToRow(gradient, v_prev->GlobalIndex(), lenWeight * deriv1_prev);
        AddToRow(gradient, v_i->GlobalIndex(), lenWeight * (deriv1_mid + deriv2_mid));
        AddToRow(gradient, v_next->GlobalIndex(), lenWeight * deriv2_next);
    }
}

// AreaPotential::AreaPotential(double wt) {
//     weight = wt;
// }

PinBendingPotential::PinBendingPotential(double wt) {
    weight = wt;
}

double PinBendingPotential::CurrentValue(const PolyCurveNetwork* curves) const {
    int nVerts = curves->NumVertices();
    double energy = 0;

    for (int i = 0; i < nVerts; i++) {
        CurveVertex* v_i = curves->GetVertex(i);
        if (v_i->numEdges() != 2) continue;
        CurveEdge* e1 = v_i->edge(0);
        CurveEdge* e2 = v_i->edge(1);

        Vector3 vec1 = e1->Vector();
        Vector3 vec2 = e2->Vector();
        double pair = 1 - dot(vec1.normalize(), vec2.normalize());
        energy += pair * pair;
    }

    return 0.5 * weight * energy;
}

void PinBendingPotential::AddGradient(const PolyCurveNetwork* curves, Eigen::MatrixXd& gradient) const {
    int nVerts = curves->NumVertices();

    for (int i = 0; i < nVerts; i++) {
        const CurveVertex* v_i = curves->GetVertex(i);
        if (v_i->numEdges() != 2) continue;
        const CurveEdge* e1 = v_i->edge(0);
        const CurveEdge* e2 = v_i->edge(1);

        const Vector3 vec1 = e1->Vector();
        const Vector3 vec2 = e2->Vector();
        const double pair = 1 - dot(vec1.normalize(), vec2.normalize());

        const CurveVertex* v_prev = e1->Opposite(v_i);
        const CurveVertex* v_next = e2->Opposite(v_i);

        // We want to differentiate dot(vec1_hat, vec2_hat)
        const VertJacobian vec1_deriv_prev = TPESC::edge_tangent_wrt_vert(e1, v_prev);
        const VertJacobian vec1_deriv_i = TPESC::edge_tangent_wrt_vert(e1, v_i);
        const VertJacobian vec2_deriv_i = TPESC::edge_tangent_wrt_vert(e2, v_i);
        const VertJacobian vec2_deriv_next = TPESC::edge_tangent_wrt_vert(e2, v_next);

        const Vector3 deriv_prev = vec1_deriv_prev.LeftMultiply(vec2);
        const Vector3 deriv_i = vec1_deriv_i.LeftMultiply(vec2) + vec2_deriv_i.LeftMultiply(vec1);
        const Vector3 deriv_next = vec2_deriv_next.LeftMultiply(vec1);

        const double pairWeight = pair * weight;
        AddToRow(gradient, v_prev->GlobalIndex(), pairWeight * deriv_prev);
        AddToRow(gradient, v_i->GlobalIndex(), pairWeight * deriv_i);
        AddToRow(gradient, v_next->GlobalIndex(), pairWeight * deriv_next);
    }
}

VectorField::VectorField() {
}

VectorField::~VectorField() {
}

Vector3 VectorField::Sample(const Vector3& x) const {
    return Vector3::zero();
}

VertJacobian VectorField::SpatialDerivative(const Vector3& x) const {
    return VertJacobian{Vector3{0, 0, 0}, Vector3{0, 0, 0}, Vector3{0, 0, 0}};
}

ConstantVectorField::ConstantVectorField(const Vector3& v) {
    c = v.normalize();
}

Vector3 ConstantVectorField::Sample(const Vector3& x) const {
    return c;
}

VertJacobian ConstantVectorField::SpatialDerivative(const Vector3& x) const {
    return VertJacobian{Vector3{0, 0, 0}, Vector3{0, 0, 0}, Vector3{0, 0, 0}};
}

CircularVectorField::CircularVectorField() {
}

Vector3 CircularVectorField::Sample(const Vector3& x) const {
    if (x.x == 0 && x.z == 0) return Vector3{0, 0, 0};

    Vector3 outward{x.x, 0, x.z};
    outward.normalize();
    Vector3 up{0, 1, 0};
    Vector3 circle = cross(up, outward).normalize();
    return circle;
}

VertJacobian CircularVectorField::SpatialDerivative(const Vector3& x) const {
    Vector3 dir_x = dirDeriv(x, Vector3{1, 0, 0});
    Vector3 dir_y = dirDeriv(x, Vector3{0, 1, 0});
    Vector3 dir_z = dirDeriv(x, Vector3{0, 0, 1});
    return VertJacobian{dir_x, dir_y, dir_z};
}

Vector3 CircularVectorField::dirDeriv(const Vector3& x, const Vector3& dir) const {
    Vector3 T = Sample(x);
    double radius = x.norm();
    double rate = dot(T, dir);
    return (rate / radius) * -x.normalize();
}

// Evaluates the gradient of the function f (defined above) at x
Vector3 gradF(const Vector3& x) {
    return x.normalize(); // unit normal to the sphere
}


Vector3 gradGaussian(double a, double b, const Vector3& c, const Vector3& x) {
    const Vector3 r = (x - c) / b;
    return -2. * (a / (b * b)) * exp(-r.norm2()) * (x - c);
}

// Evaluates a tangent vector field at x
Vector3 coolField(const Vector3& x) {
    // constants used to define the vector field (via a sum of Gaussians)
    const int nTerms = 2;
    double as[nTerms] = {0.33570404546960386, 0.6339815359591325}; // amplitudes
    double bs[nTerms] = {0.45645123825263556, 0.39270985340263576}; // widths
    double cxs[nTerms] = {0.6127842849933675, -0.5468629255159502}; // centers x
    double cys[nTerms] = {-0.7513894319788624, -0.8331733676192284}; // centers y
    double czs[nTerms] = {0.24476384858808273, 0.08223794857710300}; // centers z

    // get a vector by summing the gradient of Gaussians
    Vector3 X{0., 0., 0.};
    for (int i = 0; i < nTerms; i++) {
        Vector3 csi{cxs[i], cys[i], czs[i]};
        X += gradGaussian(as[i], bs[i], csi, x);
    }

    // project onto tangent space at x
    Vector3 N = gradF(x).normalize(); // unit surface normal
    X -= dot(X, N) * N;

    // rotate by 90 degrees around normal
    X = cross(X, N);

    // normalize
    X = X.normalize();

    return X;
}

InterestingVectorField::InterestingVectorField() {
}

Vector3 InterestingVectorField::Sample(const Vector3& x) const {
    return coolField(x);
}

VertJacobian InterestingVectorField::SpatialDerivative(const Vector3& x) const {
    return VertJacobian{Vector3{0, 0, 0}, Vector3{0, 0, 0}, Vector3{0, 0, 0}};
}

VectorFieldPotential::VectorFieldPotential(double wt, VectorField* vf) {
    weight = wt;
    field = vf;
}

VectorFieldPotential::~VectorFieldPotential() {
    delete field;
}

double VectorFieldPotential::CurrentValue(const PolyCurveNetwork* curves) const {
    int nEdges = curves->NumEdges();
    double sum = 0;

    for (int i = 0; i < nEdges; i++) {
        CurveEdge* e_i = curves->GetEdge(i);
        Vector3 center = e_i->Midpoint();
        Vector3 tangent = e_i->Tangent();
        Vector3 vf = field->Sample(center);
        Vector3 vxt = cross(vf, tangent);
        sum += dot(vxt, vxt);
    }

    return sum * weight;
}

void VectorFieldPotential::AddGradient(const PolyCurveNetwork* curves, Eigen::MatrixXd& gradient) const {
    int nEdges = curves->NumEdges();

    for (int i = 0; i < nEdges; i++) {
        CurveEdge* e_i = curves->GetEdge(i);
        // Midpoint of edge
        Vector3 center = e_i->Midpoint();
        // Tangent at edge
        Vector3 tangent = e_i->Tangent();
        // Sample the vector field at edge midpoint
        Vector3 vf = field->Sample(center);
        Vector3 vxt = cross(vf, tangent);

        // Cross product matrix from tangent
        Eigen::Matrix3d Y_skw, T_skw, deriv_Y;
        VertJacobian deriv_Y_J = field->SpatialDerivative(center);
        CrossMatrix(vf, Y_skw);
        deriv_Y_J.FillMatrix(deriv_Y);
        CrossMatrix(tangent, T_skw);

        CurveVertex* v_prev = e_i->prevVert;
        CurveVertex* v_next = e_i->nextVert;

        // Derivatives of lengths and trangents
        VertJacobian ddx_t_prev_J = TPESC::edge_tangent_wrt_vert(e_i, v_prev);
        VertJacobian ddx_t_next_J = TPESC::edge_tangent_wrt_vert(e_i, v_next);
        Vector3 deriv_len_prev = TPESC::edge_length_wrt_vert(e_i, v_prev);
        Vector3 deriv_len_next = TPESC::edge_length_wrt_vert(e_i, v_next);

        // Transfer tangent Jacobians into Eigen 3x3 matrices
        Eigen::Matrix3d ddx_t_prev, ddx_t_next;
        ddx_t_prev_J.FillMatrix(ddx_t_prev);
        ddx_t_next_J.FillMatrix(ddx_t_next);

        // Transfer length gradients into Eigen 3-vectors
        Eigen::Vector3d vxt_vec, ddx_len_prev, ddx_len_next;
        vxt_vec(0) = vxt.x;
        vxt_vec(1) = vxt.y;
        vxt_vec(2) = vxt.z;
        ddx_len_prev(0) = deriv_len_prev.x;
        ddx_len_prev(1) = deriv_len_prev.y;
        ddx_len_prev(2) = deriv_len_prev.z;
        ddx_len_next(0) = deriv_len_next.x;
        ddx_len_next(1) = deriv_len_next.y;
        ddx_len_next(2) = deriv_len_next.z;

        // Derivative of (norm cross product)^2 = (cross product) * (cross product)
        Eigen::Vector3d deriv_prev = -2 * (Y_skw * ddx_t_prev - 0.5 * T_skw * deriv_Y) * vxt_vec;
        Eigen::Vector3d deriv_next = -2 * (Y_skw * ddx_t_next - 0.5 * T_skw * deriv_Y) * vxt_vec;
        double normcross2 = dot(vxt, vxt);
        double len = e_i->Length();

        gradient.row(v_prev->GlobalIndex()) += ((ddx_len_prev * normcross2) + (len * deriv_prev)) * weight;
        gradient.row(v_next->GlobalIndex()) += ((ddx_len_next * normcross2) + (len * deriv_next)) * weight;
    }
}
}