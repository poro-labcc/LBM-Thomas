#include "./Boundary.h"
#include <cmath>
#include <omp.h>
#include "../SimStructure/Constants.h"
#include "./BoundaryConditions.h"

void Boundary(LBMParams &params, BoundaryConditionType type) {
#pragma omp parallel for
    /*
    for (int j = 1; j < Ny-1; j++) {
        double vx = (-params.uo / 25281.0) * (j - 0.5) * (j - (Ny - 1.5));
        double vy = 0.0;

        double rhow = 2*params.rho[index2d(1, j)] - params.rho[index2d(2,j)];

        double t1 = vx * vx + vy * vy;
        for (int k = 0; k < K; k++) {
            double t2 = vx * params.cx[k] + vy * params.cy[k]; // Corrected: use vy instead of params.v[index2d(i,j)]
            params.f[index3D(k, 0, j)] =
                    rhow * params.w[k] * (1.0 + 3.0 * t2 + 4.5 * t2 * t2 - 1.5 * t1);
        }
    }
    */
    for (int j = 1; j < Ny-1; j++) {
        double vx = (-params.uo / 25281.0) * (j - 0.5) * (j - (Ny - 1.5));
        double rhow = (params.f[index3D(0,0,j)] + params.f[index3D(2,0,j)] + params.f[index3D(4,0,j)]
            + 2*(params.f[index3D(3,0,j)] + params.f[index3D(6,0,j)] +
                params.f[index3D(7,0,j)]))/(1-vx);

        params.f[index3D(1,0,j)] = params.f[index3D(3,0,j)] + 2.0/3.0 * rhow * vx;
        params.f[index3D(8,0,j)] = params.f[index3D(6,0,j)] + 0.5 * (params.f[index3D(2,0,j)] -
            params.f[index3D(4,0,j)]) - 1.0/6.0 * rhow * vx;
        params.f[index3D(5,0,j)] = params.f[index3D(7,0,j)] - 0.5 * (params.f[index3D(2,0,j)] -
            params.f[index3D(4,0,j)]) - 1.0/6.0 * rhow * vx;
    }
    applyBoundaryCondition(params, type);
}

