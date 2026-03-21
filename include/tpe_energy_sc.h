#pragma once

#include "vert_jacobian.h"
#include "poly_curve_network.h"

namespace LWS {
enum class TPEPointType {
    Cluster, Point, Edge
};

struct TangentMassPoint {
    Vector3 tangent;
    double mass;
    Vector3 point;
    CurveVertex* curvePt;
    CurveVertex* curvePt2;

    TPEPointType type() const;
};

class TPESC {
public:
    static void FillGradientSingle(const PolyCurveNetwork* curveNetwork, Eigen::MatrixXd& gradients,
                                   int i, int j, double alpha, double beta);
    static void FillGradientVectorDirect(const PolyCurveNetwork* curveNetwork, Eigen::MatrixXd& gradients,
                                         double alpha, double beta);

    static double tpe_Kf(const CurveVertex* i, const CurveVertex* j, double alpha, double beta);
    static inline double tpe_Kf_pts(const Vector3& p_x, const Vector3& p_y, const Vector3& tangent_x,
                                    double alpha, double beta);
    static inline double tpe_Kf_pts_sym(const Vector3& p_x, const Vector3& p_y,
                                        const Vector3& tangent_x, const Vector3& tangent_y,
                                        double alpha, double beta);

    static double tpe_pair(const CurveVertex* i, const CurveVertex* j, double alpha, double beta);
    static double tpe_pair_pts(const Vector3& p_x, const Vector3& p_y, const Vector3& tangent_x,
                               double l_x, double l_y,
                               double alpha, double beta);
    static double tpe_total(const PolyCurveNetwork* curves, double alpha, double beta);

    static Vector3 tpe_grad_Kf(const CurveVertex* i, const CurveVertex* j, double alpha, double beta,
                               const CurveVertex* wrt);
    static Vector3 tpe_grad_Kf(const TangentMassPoint& i, const CurveVertex* j, double alpha, double beta,
                               const CurveVertex* wrt);
    static Vector3 tpe_grad_Kf(const CurveVertex* i, const TangentMassPoint& j, double alpha, double beta,
                               const CurveVertex* wrt);

    static Vector3 tpe_grad(const CurveVertex* x, const CurveVertex* y, double alpha, double beta,
                            const CurveVertex* wrt);
    static Vector3 tpe_grad(const TangentMassPoint& x, const CurveVertex* y, double alpha, double beta,
                            const CurveVertex* wrt);
    static Vector3 tpe_grad(const CurveVertex* x, const TangentMassPoint& y, double alpha, double beta,
                            const CurveVertex* wrt);

    // Projection of (f(i) - f(j)) onto the plane normal to the tangent at i
    static Vector3 proj_normal_plane(const CurveVertex* i, const CurveVertex* j);

    // Gradient of |proj_N|^alpha
    static Vector3 grad_norm_proj_alpha(const CurveVertex* i, const CurveVertex* j, double alpha, double beta,
                                        const CurveVertex* wrt);
    static Vector3 grad_norm_proj_alpha(const TangentMassPoint& i, const CurveVertex* j, double alpha,
                                        double beta,
                                        const CurveVertex* wrt);
    static Vector3 grad_norm_proj_alpha(const CurveVertex* i, const TangentMassPoint& j, double alpha,
                                        double beta,
                                        const CurveVertex* wrt);
    static Vector3 grad_norm_proj_num(const CurveVertex* i, const CurveVertex* j, double alpha, double beta,
                                      const CurveVertex* wrt, double h);

    // Jacobian of proj_N = ((f(x) - f(y)) - <f(x) - f(y), T> * T)
    static VertJacobian grad_tangent_proj(const CurveVertex* i, const CurveVertex* j, const CurveVertex* wrt);
    static VertJacobian grad_tangent_proj(const TangentMassPoint& i, const CurveVertex* j,
                                          const CurveVertex* wrt);
    static VertJacobian grad_tangent_proj(const CurveVertex* i, const TangentMassPoint& j,
                                          const CurveVertex* wrt);
    static VertJacobian grad_tangent_proj_num(const CurveVertex* i, const CurveVertex* j,
                                              const CurveVertex* wrt, double h);

    // Gradient of edge length wrt a vertex
    static Vector3 edge_length_wrt_vert(const CurveEdge* edge, const CurveVertex* wrt);
    // Gradient of vertex dual length wrt a vertex
    static Vector3 length_wrt_vert(const CurveVertex* lengthVert, const CurveVertex* wrt);

    // Jacobian of the edge "tangent" (i.e. the edge vector) wrt a vertex
    static VertJacobian edge_tangent_wrt_vert(const CurveEdge* edge, const CurveVertex* wrtVert);
    // Jacobian of the vertex tangent (average of surrounding edges) wrt a vertex
    static VertJacobian vertex_tangent_wrt_vert(const CurveVertex* tangentVert, const CurveVertex* wrtVert);
};

inline double TPESC::tpe_Kf_pts(const Vector3& p_x, const Vector3& p_y, const Vector3& tangent_x,
                                double alpha, double beta) {
    Vector3 disp = p_x - p_y;
    Vector3 n_proj = disp - dot(disp, tangent_x) * tangent_x;
    double numer = pow(norm(n_proj), alpha);
    double denom = pow(norm(disp), beta);
    return numer / denom;
}

inline double TPESC::tpe_Kf_pts_sym(const Vector3& p_x, const Vector3& p_y,
                                    const Vector3& tangent_x, const Vector3& tangent_y,
                                    double alpha, double beta) {
    Vector3 disp = p_x - p_y;
    Vector3 n_proj_x = disp - dot(disp, tangent_x) * tangent_x;
    Vector3 n_proj_y = disp - dot(disp, tangent_y) * tangent_y;

    double numer_x = pow(norm(n_proj_x), alpha);
    double numer_y = pow(norm(n_proj_y), alpha);
    double denom = pow(norm(disp), beta);

    return 0.5 * (numer_x + numer_y) / denom;
}
}
