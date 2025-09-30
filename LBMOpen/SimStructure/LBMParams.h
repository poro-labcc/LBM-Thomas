#ifndef LBMPARAMS_H
#define LBMPARAMS_H

#include <vector>
#include "Constants.h"

struct LBMParams {
    int Nx, Ny, K;
    std::vector<double> f, feq, f_last, f_temp, rho, u, v,dist;
    std::vector<double> cx, cy, w, u_old;
    double omega, uo, rhoo;
    std::vector<bool> isSolid;
    bool stableFlow;
    std::vector<double> f_neq;

    LBMParams(int nx, int ny, int k) ;
};

#endif // LBMPARAMS_H

