#ifndef GAUSSQUAD_H_
#define GAUSSQUAD_H_

#include "typedefs.h"
#include "../eigen/Eigen/Core"

using namespace Eigen;

// Returns Gauss quadratures for 3 node triangle (fig. 3, p. 9)
// @param N_ip: Number of integration points [int]
// Returns r_i as 1st col, s_1 as 2nd col, w_i as 3rd col,
// and one row per integration point.
Matrix gausst3 (int N_ip) {

    if (N_ip == 1) {
        Matrix <Float, 1, 3> m;
        return (m << 1.0/6.0, 1.0/6.0, 1.0/6.0);
    }



}



#endif
