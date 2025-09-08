#ifndef LBMPARAMS_H
#define LBMPARAMS_H

#include <vector>
#include "../SimStructure/Constants.h"

struct LBMParams {
    int Nx, Ny, K;
    std::vector<double> f, feq, f_last, f_temp, rho, u, v,dist;
    std::vector<double> cx, cy, w;
    double omega, uo, rhoo;
    std::vector<bool> isSolid;
    bool stableFlow;

    LBMParams(int nx, int ny, int k) ;
};

#endif // LBMPARAMS_H

