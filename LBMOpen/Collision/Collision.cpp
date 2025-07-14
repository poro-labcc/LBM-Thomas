#include "./Collision.h"
#include <omp.h>
#include "../SimStructure/Constants.h"

void Collision(LBMParams &params) {
#pragma omp parallel for collapse(2)
    for (int j = 0; j < Ny; j++) {
        for (int i = 0; i < Nx; i++) {
            if (!params.isSolid[index2d(i, j)]) {
                double t1 = params.u[index2d(i, j)] * params.u[index2d(i, j)] + params.v[index2d(i, j)] * params.v[index2d(i, j)];
                for (int k = 0; k < K; k++) {
                    double t2 = params.u[index2d(i, j)] * params.cx[k] + params.v[index2d(i, j)] * params.cy[k];
                    params.feq[index3D(k, i, j)] =
                            params.rho[index2d(i, j)] * params.w[k] * (1.0 + 3.0 * t2 + 4.5 * t2 * t2 - 1.5 * t1);
                    params.f[index3D(k, i, j)] =
                            params.omega * params.feq[index3D(k, i, j)] + (1.0 - params.omega) * params.f[index3D(k, i, j)]; //Relaxation Step
                }
            }
        }
    }
}

