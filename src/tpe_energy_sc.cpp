#include "tpe_energy_sc.h"

namespace LWS {
TPEPointType TangentMassPoint::type() const {
    if (curvePt) {
        return TPEPointType::Point;
    } else {
        return TPEPointType::Cluster;
    }
}


void TPESC::FillGradientSingle(const PolyCurveNetwork* curveNetwork, Eigen::MatrixXd& gradients, int i, int j,
                               double alpha, double beta) {
    if (i == j) return;
    CurveVertex* i_pt = curveNetwork->GetVertex(i);
    CurveVertex* j_pt = curveNetwork->GetVertex(j);

    // Add i and neighbors of i
    std::vector<CurveVertex*> i_pts;
    i_pts.push_back(i_pt);
    for (int e = 0; e < i_pt->numEdges(); e++) {
        i_pts.push_back(i_pt->edge(e)->Opposite(i_pt));
    }

    // Add j and neighbors of j
    std::vector<CurveVertex*> j_pts;
    j_pts.push_back(j_pt);
    for (int e = 0; e < j_pt->numEdges(); e++) {
        j_pts.push_back(j_pt->edge(e)->Opposite(j_pt));
    }

    // Differentiate wrt neighbors of i
    for (CurveVertex* i_n : i_pts) {
        AddToRow(gradients, i_n->GlobalIndex(), TPESC::tpe_grad(i_pt, j_pt, alpha, beta, i_n));
    }
    // Differentiate wrt neighbors of j
    for (CurveVertex* j_n : j_pts) {
        bool noOverlap = true;
        // Only compute this derivative if j_n is not already included in one of the previous pairs
        for (CurveVertex* i_n : i_pts) {
            if (i_n == j_n) noOverlap = false;
        }
        if (noOverlap) {
            AddToRow(gradients, j_n->GlobalIndex(), TPESC::tpe_grad(i_pt, j_pt, alpha, beta, j_n));
        }
    }
}

void TPESC::FillGradientVectorDirect(const PolyCurveNetwork* curveNetwork, Eigen::MatrixXd& gradients,
                                     double alpha, double beta) {
    const int nVerts = curveNetwork->NumVertices();
    // Fill with zeros, so that the constraint entries are 0
    gradients.setZero();
    // Fill vertex entries with accumulated gradients
    for (int i = 0; i < nVerts; i++) {
        for (int j = 0; j < nVerts; j++) {
            if (i == j) continue;
            TPESC::FillGradientSingle(curveNetwork, gradients, i, j, alpha, beta);
        }
    }
}

double TPESC::tpe_Kf(const CurveVertex* i, const CurveVertex* j, double alpha, double beta) {
    if (i == j) return 0;
    else if (i->numEdges() > 2) return 0;

    const Vector3 disp = i->Position() - j->Position();
    const Vector3 T_i = i->Tangent();

    const Vector3 normal_proj = disp - dot(disp, T_i) * T_i;
    const double numer = pow(norm(normal_proj), alpha);
    const double denom = pow(norm(disp), beta);
    return numer / denom;
}

double TPESC::tpe_pair(const CurveVertex* i, const CurveVertex* j, double alpha, double beta) {
    const double kfxy = tpe_Kf(i, j, alpha, beta);
    const double l_x = i->DualLength();
    const double l_y = j->DualLength();
    return l_x * l_y * kfxy;
}

double TPESC::tpe_pair_pts(const Vector3& p_x, const Vector3& p_y, const Vector3& tangent_x,
                           double l_x, double l_y, double alpha, double beta) {
    const double kfxy = tpe_Kf_pts(p_x, p_y, tangent_x, alpha, beta);
    return l_x * l_y * kfxy;
}

double TPESC::tpe_total(const PolyCurveNetwork* curves, double alpha, double beta) {
    const int nVerts = curves->NumVertices();
    double sumEnergy = 0;

    for (int i = 0; i < nVerts; i++) {
        for (int j = 0; j < nVerts; j++) {
            CurveVertex* pt_i = curves->GetVertex(i);
            CurveVertex* pt_j = curves->GetVertex(j);
            sumEnergy += tpe_pair(pt_i, pt_j, alpha, beta);
        }
    }
    return sumEnergy;
}

Vector3 TPESC::proj_normal_plane(const CurveVertex* i, const CurveVertex* j) {
    const Vector3 disp = i->Position() - j->Position();
    const Vector3 T_i = i->Tangent();
    return disp - dot(disp, T_i) * T_i;
}

Vector3 TPESC::tpe_grad_Kf(const CurveVertex* i, const CurveVertex* j, double alpha, double beta,
                           const CurveVertex* wrt) {
    if (i == j) return Vector3::zero();
    else if (i->numEdges() > 2) return Vector3::zero();

    // Get positions and displacement vectors
    const Vector3 disp = i->Position() - j->Position();
    const Vector3 T_i = i->Tangent();
    // Normalized displacement direction
    const Vector3 unit_disp = disp.normalize();

    // Evaluate projection onto normal plane, v - <v,T> * T
    const Vector3 normal_proj = disp - dot(disp, T_i) * T_i;
    // Numerator of energy is norm of projection ^ alpha
    const double A = pow(norm(normal_proj), alpha);
    // Denominator of energy is distance between points ^ beta
    const double B = pow(norm(disp), beta);

    // Derivative of numerator
    const Vector3 deriv_A = grad_norm_proj_alpha(i, j, alpha, beta, wrt);
    // Derivative of denominator
    Vector3 deriv_B = Vector3::zero();
    if (wrt == i) {
        deriv_B = beta * pow(norm(disp), beta - 1) * unit_disp;
    } else if (wrt == j) {
        deriv_B = -beta * pow(norm(disp), beta - 1) * unit_disp;
    }
    // Quotient rule for A / B
    Vector3 total = (deriv_A * B - A * deriv_B) / (B * B);
    return total;
}

Vector3 TPESC::tpe_grad_Kf(const TangentMassPoint& i, const CurveVertex* j, double alpha, double beta,
                           const CurveVertex* wrt) {
    // Get positions and displacement vectors
    const Vector3 disp = i.point - j->Position();
    const Vector3 T_i = i.tangent;
    // Normalized displacement direction
    const Vector3 unit_disp = disp.normalize();

    // Evaluate projection onto normal plane, v - <v,T> * T
    const Vector3 normal_proj = disp - dot(disp, T_i) * T_i;
    // Numerator of energy is norm of projection ^ alpha
    const double A = pow(norm(normal_proj), alpha);
    // Denominator of energy is distance between points ^ beta
    const double B = pow(norm(disp), beta);

    // Derivative of numerator
    const Vector3 deriv_A = grad_norm_proj_alpha(i, j, alpha, beta, wrt);

    // Derivative of denominator
    Vector3 deriv_B = Vector3::zero();
    if (i.type() == TPEPointType::Point && wrt == i.curvePt) {
        deriv_B = beta * pow(norm(disp), beta - 1) * unit_disp;
    } else if (i.type() == TPEPointType::Edge && (wrt == i.curvePt || wrt == i.curvePt2)) {
        deriv_B = (beta * pow(norm(disp), beta - 1) * unit_disp) / 2;
    }
    if (wrt == j) {
        deriv_B += -beta * pow(norm(disp), beta - 1) * unit_disp;
    }
    // Quotient rule for A / B
    return (deriv_A * B - A * deriv_B) / (B * B);
}

Vector3 TPESC::tpe_grad_Kf(const CurveVertex* i, const TangentMassPoint& j, double alpha, double beta,
                           const CurveVertex* wrt) {
    // Get positions and displacement vectors
    const Vector3 disp = i->Position() - j.point;
    const Vector3 T_i = i->Tangent();
    // Normalized displacement direction
    const Vector3 unit_disp = disp.normalize();

    // Evaluate projection onto normal plane, v - <v,T> * T
    const Vector3 normal_proj = disp - dot(disp, T_i) * T_i;
    // Numerator of energy is norm of projection ^ alpha
    const double A = pow(norm(normal_proj), alpha);
    // Denominator of energy is distance between points ^ beta
    const double B = pow(norm(disp), beta);

    // Derivative of numerator
    const Vector3 deriv_A = grad_norm_proj_alpha(i, j, alpha, beta, wrt);
    // Derivative of denominator
    Vector3 deriv_B = Vector3::zero();
    if (wrt == i) {
        deriv_B = beta * pow(norm(disp), beta - 1) * unit_disp;
    }
    if (j.type() == TPEPointType::Point && wrt == j.curvePt) {
        deriv_B += -beta * pow(norm(disp), beta - 1) * unit_disp;
    } else if (j.type() == TPEPointType::Edge && (wrt == j.curvePt || wrt == j.curvePt2)) {
        deriv_B += (-beta * pow(norm(disp), beta - 1) * unit_disp) / 2;
    }

    // Quotient rule for A / B
    return (deriv_A * B - A * deriv_B) / (B * B);
}

Vector3 TPESC::tpe_grad(const CurveVertex* x, const CurveVertex* y, double alpha, double beta,
                        const CurveVertex* wrt) {
    // Computes the gradient of the kernel (K_f(x, y) dx dy) with respect to
    // the position of the vertex "wrt".
    if (x->numEdges() > 2) return Vector3{0, 0, 0};

    // First get the gradient of K_f(x, y)
    const Vector3 grad_Kf = tpe_grad_Kf(x, y, alpha, beta, wrt);
    const double Kf = tpe_Kf(x, y, alpha, beta);

    const double l_x = x->DualLength();
    const double l_y = y->DualLength();

    // d/dy of area(x)
    const Vector3 grad_lx = length_wrt_vert(x, wrt);
    // d/dy of area(y)
    const Vector3 grad_ly = length_wrt_vert(y, wrt);
    // Evaluate the product rule for dx*dy
    const Vector3 prod_rule = grad_lx * l_y + l_x * grad_ly;
    // Evaluate the product rule for k dx dy
    const Vector3 total = grad_Kf * l_x * l_y + Kf * prod_rule;
    return total;
}

Vector3 TPESC::tpe_grad(const TangentMassPoint& x, const CurveVertex* y, double alpha, double beta,
                        const CurveVertex* wrt) {
    // Computes the gradient of the kernel (K_f(x, y) dx dy) with respect to
    // the position of the vertex "wrt".

    // First get the gradient of K_f(x, y)
    const Vector3 grad_Kf = tpe_grad_Kf(x, y, alpha, beta, wrt);
    const double Kf = tpe_Kf_pts(x.point, y->Position(), x.tangent, alpha, beta);

    const double l_x = x.mass;
    const double l_y = y->DualLength();

    // Area gradient for x depends on whether the mass point is distant
    Vector3 grad_lx{0, 0, 0};
    if (x.type() == TPEPointType::Point) {
        grad_lx = length_wrt_vert(x.curvePt, wrt);
    } else if (x.type() == TPEPointType::Edge) {
        if (wrt == x.curvePt) {
            grad_lx = (x.curvePt->Position() - x.curvePt2->Position());
            grad_lx = grad_lx.normalize();
        } else if (wrt == x.curvePt2) {
            grad_lx = (x.curvePt2->Position() - x.curvePt->Position());
            grad_lx = grad_lx.normalize();
        }
    }
    // d/dy of area(y)
    const Vector3 grad_ly = length_wrt_vert(y, wrt);
    // Evaluate the product rule for dx*dy
    const Vector3 prod_rule = grad_lx * l_y + l_x * grad_ly;

    // Evaluate the product rule for k dx dy
    return grad_Kf * l_x * l_y + Kf * prod_rule;
}

Vector3 TPESC::tpe_grad(const CurveVertex* x, const TangentMassPoint& y, double alpha, double beta,
                        const CurveVertex* wrt) {
    // Here, y is a mass point (e.g. from Barnes-Hut), so gradients of y are assumed to be zero.
    // First get the gradient of K_f(x, y)
    const Vector3 grad_Kf = tpe_grad_Kf(x, y, alpha, beta, wrt);
    const double Kf = tpe_Kf_pts(x->Position(), y.point, x->Tangent(), alpha, beta);

    const double l_x = x->DualLength();
    const double l_y = y.mass;

    // d/dy of area(x)
    const Vector3 grad_lx = length_wrt_vert(x, wrt);
    // Area gradient for y depends on whether the mass point is distant
    Vector3 grad_ly = Vector3::zero();
    if (y.type() == TPEPointType::Point) {
        grad_ly = length_wrt_vert(y.curvePt, wrt);
    } else if (y.type() == TPEPointType::Edge) {
        if (wrt == y.curvePt) {
            grad_ly = (y.curvePt->Position() - y.curvePt2->Position());
            grad_ly = grad_ly.normalize();
        } else if (wrt == y.curvePt2) {
            grad_ly = (y.curvePt2->Position() - y.curvePt->Position());
            grad_ly = grad_ly.normalize();
        }
    }

    // Evaluate the product rule for dx*dy
    const Vector3 prod_rule = grad_lx * l_y + l_x * grad_ly;

    // Evaluate the product rule for k dx dy
    return grad_Kf * l_x * l_y + Kf * prod_rule;
}

Vector3 TPESC::grad_norm_proj_alpha(const CurveVertex* i, const CurveVertex* j, double alpha, double beta,
                                    const CurveVertex* wrt) {
    static const Vector3 zero = Vector3::zero();
    const Vector3 disp = i->Position() - j->Position();
    const Vector3 T_i = i->Tangent();
    // Projection onto normal plane
    const Vector3 normal_proj = disp - dot(disp, T_i) * T_i;
    const double proj_len = norm(normal_proj);

    // If the displacement is actually exactly perpendicular to the tangent,
    // then the contribution is exactly 0.
    if (proj_len < 1e-10) return zero;

    // Derivative of |f(x) - ...|^alpha = alpha * |f(x) - ...|^(alpha - 1)
    const double alpha_deriv = alpha * pow(proj_len, alpha - 1);
    // Normalized vector of projection onto normal plane
    const Vector3 proj_normalized = normal_proj / proj_len;

    VertJacobian deriv_disp{zero, zero, zero};
    if (wrt == i) {
        deriv_disp.directional_x = Vector3{1, 0, 0};
        deriv_disp.directional_y = Vector3{0, 1, 0};
        deriv_disp.directional_z = Vector3{0, 0, 1};
    } else if (wrt == j) {
        deriv_disp.directional_x = Vector3{-1, 0, 0};
        deriv_disp.directional_y = Vector3{0, -1, 0};
        deriv_disp.directional_z = Vector3{0, 0, -1};
    }

    // Derivative of <f(x) - f(y), T> * T
    const VertJacobian deriv_T_inner = grad_tangent_proj(i, j, wrt);
    const VertJacobian deriv_N_proj = deriv_disp - deriv_T_inner;

    const Vector3 total = alpha_deriv * deriv_N_proj.LeftMultiply(proj_normalized);
    return total;
}

Vector3 TPESC::grad_norm_proj_alpha(const TangentMassPoint& i, const CurveVertex* j, double alpha,
                                    double beta,
                                    const CurveVertex* wrt) {
    static const Vector3 zero = Vector3::zero();
    const Vector3 disp = i.point - j->Position();
    const Vector3 T_i = i.tangent;
    // Projection onto normal plane
    const Vector3 normal_proj = disp - dot(disp, T_i) * T_i;
    const double proj_len = norm(normal_proj);

    // If the displacement is actually exactly perpendicular to the tangent,
    // then the contribution is exactly 0.
    if (proj_len < 1e-10) return zero;

    // Derivative of |f(x) - ...|^alpha = alpha * |f(x) - ...|^(alpha - 1)
    const double alpha_deriv = alpha * pow(proj_len, alpha - 1);
    // Normalized vector of projection onto normal plane
    const Vector3 proj_normalized = normal_proj / proj_len;

    VertJacobian deriv_disp{zero, zero, zero};
    if (i.type() == TPEPointType::Point && wrt == i.curvePt) {
        deriv_disp.directional_x = Vector3{1, 0, 0};
        deriv_disp.directional_y = Vector3{0, 1, 0};
        deriv_disp.directional_z = Vector3{0, 0, 1};
    }
    if (i.type() == TPEPointType::Edge && (wrt == i.curvePt || wrt == i.curvePt2)) {
        deriv_disp.directional_x = Vector3{0.5, 0, 0};
        deriv_disp.directional_y = Vector3{0, 0.5, 0};
        deriv_disp.directional_z = Vector3{0, 0, 0.5};
    } else if (wrt == j) {
        deriv_disp.directional_x = Vector3{-1, 0, 0};
        deriv_disp.directional_y = Vector3{0, -1, 0};
        deriv_disp.directional_z = Vector3{0, 0, -1};
    }

    // Derivative of <f(x) - f(y), T> * T
    const VertJacobian deriv_T_inner = grad_tangent_proj(i, j, wrt);
    const VertJacobian deriv_N_proj = deriv_disp - deriv_T_inner;

    return alpha_deriv * deriv_N_proj.LeftMultiply(proj_normalized);
}

Vector3 TPESC::grad_norm_proj_alpha(const CurveVertex* i, const TangentMassPoint& j, double alpha,
                                    double beta,
                                    const CurveVertex* wrt) {
    static const Vector3 zero = Vector3::zero();
    const Vector3 disp = i->Position() - j.point;
    const Vector3 T_i = i->Tangent();
    // Projection onto normal plane
    const Vector3 normal_proj = disp - dot(disp, T_i) * T_i;
    const double proj_len = norm(normal_proj);

    // If the displacement is actually exactly perpendicular to the tangent,
    // then the contribution is exactly 0.
    if (proj_len < 1e-10) return zero;

    // Derivative of |f(x) - ...|^alpha = alpha * |f(x) - ...|^(alpha - 1)
    const double alpha_deriv = alpha * pow(proj_len, alpha - 1);
    // Normalized vector of projection onto normal plane
    const Vector3 proj_normalized = normal_proj / proj_len;

    VertJacobian deriv_disp{zero, zero, zero};
    if (wrt == i) {
        deriv_disp.directional_x = Vector3{1, 0, 0};
        deriv_disp.directional_y = Vector3{0, 1, 0};
        deriv_disp.directional_z = Vector3{0, 0, 1};
    } else if (j.type() == TPEPointType::Point && wrt == j.curvePt) {
        deriv_disp.directional_x = Vector3{-1, 0, 0};
        deriv_disp.directional_y = Vector3{0, -1, 0};
        deriv_disp.directional_z = Vector3{0, 0, -1};
    } else if (j.type() == TPEPointType::Edge && (wrt == j.curvePt || wrt == j.curvePt2)) {
        deriv_disp.directional_x = Vector3{-0.5, 0, 0};
        deriv_disp.directional_y = Vector3{0, -0.5, 0};
        deriv_disp.directional_z = Vector3{0, 0, -0.5};
    }

    // Derivative of <f(x) - f(y), T> * T
    const VertJacobian deriv_T_inner = grad_tangent_proj(i, j, wrt);
    const VertJacobian deriv_N_proj = deriv_disp - deriv_T_inner;

    return alpha_deriv * deriv_N_proj.LeftMultiply(proj_normalized);
}

Vector3 TPESC::grad_norm_proj_num(const CurveVertex* i, const CurveVertex* j, double alpha, double beta,
                                  const CurveVertex* wrt, double h) {
    const Vector3 origPos = wrt->Position();
    const double orig = pow(norm(proj_normal_plane(i, j)), alpha);

    wrt->SetPosition(origPos + Vector3{h, 0, 0});
    const double xVal = pow(norm(proj_normal_plane(i, j)), alpha);

    wrt->SetPosition(origPos + Vector3{0, h, 0});
    const double yVal = pow(norm(proj_normal_plane(i, j)), alpha);

    wrt->SetPosition(origPos + Vector3{0, 0, h});
    const double zVal = pow(norm(proj_normal_plane(i, j)), alpha);

    wrt->SetPosition(origPos);

    const double xDeriv = (xVal - orig) / h;
    const double yDeriv = (yVal - orig) / h;
    const double zDeriv = (zVal - orig) / h;

    return Vector3{xDeriv, yDeriv, zDeriv};
}

VertJacobian TPESC::grad_tangent_proj(const CurveVertex* i, const CurveVertex* j, const CurveVertex* wrt) {
    // Differentiate the inner product
    const Vector3 disp = i->Position() - j->Position();
    const Vector3 T_i = i->Tangent();
    const double disp_dot_T = dot(disp, T_i);

    Vector3 inner_deriv_A_B = Vector3::zero();
    if (wrt == i) {
        inner_deriv_A_B = T_i;
    } else if (wrt == j) {
        inner_deriv_A_B = -T_i;
    }
    const VertJacobian deriv_T = vertex_tangent_wrt_vert(i, wrt);
    const Vector3 inner_A_deriv_B = deriv_T.LeftMultiply(disp);
    const Vector3 deriv_inner = inner_deriv_A_B + inner_A_deriv_B;

    // Now use product rule for <f(x) - f(y), T)> * T
    const VertJacobian deriv_A_B = outer_product_to_jacobian(T_i, deriv_inner);
    const VertJacobian A_deriv_B = disp_dot_T * deriv_T;

    return deriv_A_B + A_deriv_B;
}

VertJacobian TPESC::grad_tangent_proj(const TangentMassPoint& i, const CurveVertex* j,
                                      const CurveVertex* wrt) {
    // Differentiate the inner product
    const Vector3 disp = i.point - j->Position();
    const Vector3 T_i = i.tangent;
    const double disp_dot_T = dot(disp, T_i);

    Vector3 inner_deriv_A_B = Vector3::zero();
    if (wrt == j) {
        inner_deriv_A_B = -T_i;
    }
    const VertJacobian deriv_T = vertex_tangent_wrt_vert(i.curvePt, wrt);
    const Vector3 inner_A_deriv_B = deriv_T.LeftMultiply(disp);
    const Vector3 deriv_inner = inner_deriv_A_B + inner_A_deriv_B;

    // Now use product rule for <f(x) - f(y), T)> * T
    const VertJacobian deriv_A_B = outer_product_to_jacobian(T_i, deriv_inner);
    const VertJacobian A_deriv_B = disp_dot_T * deriv_T;

    return deriv_A_B + A_deriv_B;
}

VertJacobian TPESC::grad_tangent_proj(const CurveVertex* i, const TangentMassPoint& j,
                                      const CurveVertex* wrt) {
    // Differentiate the inner product
    const Vector3 disp = i->Position() - j.point;
    const Vector3 T_i = i->Tangent();
    const double disp_dot_T = dot(disp, T_i);

    Vector3 inner_deriv_A_B = Vector3::zero();
    if (wrt == i) {
        inner_deriv_A_B = T_i;
    }
    const VertJacobian deriv_T = vertex_tangent_wrt_vert(i, wrt);
    const Vector3 inner_A_deriv_B = deriv_T.LeftMultiply(disp);
    const Vector3 deriv_inner = inner_deriv_A_B + inner_A_deriv_B;

    // Now use product rule for <f(x) - f(y), T)> * T
    const VertJacobian deriv_A_B = outer_product_to_jacobian(T_i, deriv_inner);
    const VertJacobian A_deriv_B = disp_dot_T * deriv_T;

    return deriv_A_B + A_deriv_B;
}

VertJacobian TPESC::grad_tangent_proj_num(const CurveVertex* i, const CurveVertex* j, const CurveVertex* wrt,
                                          double h) {
    const Vector3 origPos = wrt->Position();
    const Vector3 origTangent = dot(i->Position() - j->Position(), i->Tangent()) * i->Tangent();

    wrt->SetPosition(origPos + Vector3{h, 0, 0});
    const Vector3 xTangent = dot(i->Position() - j->Position(), i->Tangent()) * i->Tangent();

    wrt->SetPosition(origPos + Vector3{0, h, 0});
    const Vector3 yTangent = dot(i->Position() - j->Position(), i->Tangent()) * i->Tangent();

    wrt->SetPosition(origPos + Vector3{0, 0, h});
    const Vector3 zTangent = dot(i->Position() - j->Position(), i->Tangent()) * i->Tangent();

    wrt->SetPosition(origPos);

    const Vector3 xDeriv = (xTangent - origTangent) / h;
    const Vector3 yDeriv = (yTangent - origTangent) / h;
    const Vector3 zDeriv = (zTangent - origTangent) / h;

    return VertJacobian{xDeriv, yDeriv, zDeriv};
}


VertJacobian TPESC::edge_tangent_wrt_vert(const CurveEdge* edge, const CurveVertex* wrtVert) {
    // get positions
    const CurveVertex* prevVert = edge->prevVert;
    const CurveVertex* nextVert = edge->nextVert;
    const Vector3 v_h = prevVert->Position();
    const Vector3 v_i = nextVert->Position();

    if (wrtVert != prevVert && wrtVert != nextVert) {
        return VertJacobian{Vector3::zero(), Vector3::zero(), Vector3::zero()};
    }

    const Vector3 v_tangent = v_i - v_h;
    const double v_norm = norm(v_tangent);
    const Vector3 v_normalized = v_tangent / v_norm;

    constexpr VertJacobian I{Vector3{1, 0, 0}, Vector3{0, 1, 0}, Vector3{0, 0, 1}};

    const VertJacobian deriv_A_B = I * v_norm;
    const VertJacobian A_deriv_B = outer_product_to_jacobian(v_tangent, v_normalized);
    const VertJacobian deriv = (deriv_A_B - A_deriv_B) * (1.0 / (v_norm * v_norm));

    // If we're differentiating the tail vertex, the derivative is negative
    if (wrtVert == prevVert) return -1 * deriv;
        // Otherwise we're differentiating the head vertex, and the derivative is positive
    else return deriv;
}

VertJacobian TPESC::vertex_tangent_wrt_vert(const CurveVertex* tangentVert, const CurveVertex* wrtVert) {
    static const auto zero_jac = VertJacobian{Vector3::zero(), Vector3::zero(), Vector3::zero()};
    if (!tangentVert || !wrtVert) {
        return zero_jac;
    }

    if (tangentVert->numEdges() != 2) {
        if (tangentVert->numEdges() == 1) {
            const CurveEdge* edge = tangentVert->edge(0);
            //Vector3 tangent = edge->Tangent();
            // Derivative of T
            return edge_tangent_wrt_vert(edge, wrtVert);
        } else {
            return zero_jac;
        }
    }

    const CurveEdge* prevEdge = tangentVert->edge(0);
    const CurveEdge* nextEdge = tangentVert->edge(1);

    const Vector3 prevTangent = prevEdge->Tangent();
    const Vector3 nextTangent = nextEdge->Tangent();

    const Vector3 sumTangents = prevTangent + nextTangent;
    const double normSum = norm(sumTangents);
    const Vector3 vertTangent = sumTangents.normalize();

    // Quotient rule on (T1 + T2) / |T1 + T2|
    const VertJacobian derivSumTs = edge_tangent_wrt_vert(prevEdge, wrtVert)
                                    + edge_tangent_wrt_vert(nextEdge, wrtVert);

    const Vector3 derivNorm = derivSumTs.LeftMultiply(vertTangent);
    const VertJacobian deriv_A_B = derivSumTs * normSum;
    const VertJacobian A_deriv_B = outer_product_to_jacobian(sumTangents, derivNorm);

    return (deriv_A_B - A_deriv_B) * (1.0 / (normSum * normSum));
}


Vector3 TPESC::edge_length_wrt_vert(const CurveEdge* edge, const CurveVertex* wrt) {
    if (edge->prevVert != wrt && edge->nextVert != wrt) {
        return Vector3::zero();
    }
    const Vector3 wrt_pos = wrt->Position();
    const Vector3 opp_pos = edge->Opposite(wrt)->Position();

    // To increase the edge length, we want to move away from the opposite vertex
    const Vector3 outward = wrt_pos - opp_pos;
    return outward.normalize();
}

Vector3 TPESC::length_wrt_vert(const CurveVertex* lengthVert, const CurveVertex* wrt) {
    // If differentiating wrt self, need to consider both side
    if (lengthVert == wrt) {
        Vector3 sumDirections = Vector3::zero();
        Vector3 center = lengthVert->Position();
        for (int e = 0; e < lengthVert->numEdges(); e++) {
            const Vector3 other = lengthVert->edge(e)->Opposite(lengthVert)->Position();
            const Vector3 outward = (other - center).normalize();
            sumDirections += outward;
        }
        // Gradient is half the sum of unit vectors along outgoing edges
        return -0.5 * sumDirections;
    }
    // Otherwise, only consider the one edge between other and vert
    else {
        for (int e = 0; e < lengthVert->numEdges(); e++) {
            const CurveVertex* other = lengthVert->edge(e)->Opposite(lengthVert);
            if (other == wrt) {
                const Vector3 outward = (other->Position() - lengthVert->Position()).normalize();
                return outward / 2;
            }
        }
        return Vector3::zero();
    }
}
}