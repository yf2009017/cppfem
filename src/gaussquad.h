#ifndef GAUSSQUAD_H_
#define GAUSSQUAD_H_

#include "typedefs.h"
#include "../eigen/Eigen/Core"

using namespace Eigen;

// Returns Gauss quadratures for 3 node triangle (fig. 3, p. 9)
// @param N_ip: Number of integration points
// Returns r_i as 1st col, s_1 as 2nd col, w_i as 3rd col,
// and one row per integration point.
Matrix<Float, Dynamic, 3> gausst3 (int N_ip) {

    if (N_ip == 1) {
        Matrix<Float, 1, 3> m;
        m << 1.0/6.0, 1.0/6.0, 1.0/6.0;
        return m;

    } else if (N_ip == 3) {
        Matrix<Float, 3, 3> m;
        m << 1.0/6.0, 1.0/6.0, 1.0/6.0,
          2.0/3.0, 1.0/6.0, 1.0/6.0,
          1.0/6.0, 2.0/3.0, 1.0/6.0;
        return m;
        
    } else if (N_ip == 4) {
        Matrix<Float, 4, 3> m;
        m << 1.0/3.0, 1.0/3.0, -9.0/32.0,
          3.0/5.0, 1.0/5.0, 25.0/96.0,
          1.0/5.0, 3.0/5.0, 25.0/96.0,
          1.0/5.0, 1.0/5.0, 25.0/96.0;
        return m;
        
    } else if (N_ip == 7) {
        Matrix<Float, 7, 3> m;
        m << 0.0, 0.0, 1.0/40.0,
          0.5, 0.0, 1.0/15.0,
          1.0, 0.0, 1.0/40.0,
          0.5, 0.5, 1.0/15.0,
          0.0, 1.0, 1.0/40.0,
          0.0, 0.5, 1.0/15.0,
          1.0/3.0, 1.0/3.0, 9.0/40.0;
        return m;

    } else {
        throw "Error: gausst3 can only return data for 1, 3, 4, or 7 integration points (N_ip)";
    }
}

// Returns Gauss quadratures for a 2d line segment
// @param N_ip: Number of integration points
Matrix<Float, Dynamic, 2> gauss1d (int N_ip) {

    if (N_ip == 1) {
        Matrix<Float, 1, 2> m;
        m << 0.0, 2.0;
        return m;

    } else if (N_ip == 3) {
        Matrix<Float, 3, 2> m;
        m << -0.774596669241483377035835, 0.55555555555555555555556,
          0.000000000000000000000000, 0.88888888888888888888889,
          0.774596669241483377035835, 0.55555555555555555555556;
        return m;
        
    } else if (N_ip == 4) {
        Matrix<Float, 4, 2> m;
        m << -0.861136311594052575223946, 0.34785484513745385737306,
          -0.339981043584856264802666, 0.65214515486254614262694,
          0.339981043584856264802666, 0.65214515486254614262694,
          0.861136311594052575223946, 0.34785484513745385737306;
        return m;
        
    } else if (N_ip == 7) {
        Matrix<Float, 7, 2> m;
        m << -0.949107912342758524526190, 0.12948496616886969327061,
          -0.741531185599394439863865, 0.27970539148927666790147,
          -0.405845151377397166906607, 0.38183005050511894495037,
          0.000000000000000000000000, 0.41795918367346938775510,
          0.405845151377397166906607, 0.38183005050511894495037,
          0.741531185599394439863865, 0.27970539148927666790147,
          0.949107912342758524526190, 0.12948496616886969327061;
        return m;

    } else {
        throw "Error: gauss1d can only return data for 1, 3, 4, or 7 integration points (N_ip)";
    }
}

#endif
