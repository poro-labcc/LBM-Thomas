#include "./MacroRecover.h"
#include <omp.h>
#include "../SimStructure/Constants.h"

void MacroRecover(LBMParams &params) {
#pragma omp parallel for collapse(2)
    //Recovering macroscopic properties
    for (int j = 0; j < Ny; j++) {
        for (int i = 0; i < Nx; i++) {
            if (!params.isSolid[index2d(i, j)]) {
                double ssum = 0;
                double usum = 0;
                double vsum = 0;
                for (int k = 0; k < K; k++) {
                    ssum = ssum + params.f[index3D(k, i, j)];
                    usum = usum + params.f[index3D(k, i, j)] * params.cx[k];
                    vsum = vsum + params.f[index3D(k, i, j)] * params.cy[k];
                }
                params.rho[index2d(i, j)] = ssum;
                params.u[index2d(i, j)] = usum / params.rho[index2d(i, j)];
                params.v[index2d(i, j)] = vsum / params.rho[index2d(i, j)];
            }
        }
    }
}

