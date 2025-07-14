#include "./BoundaryConditions.h"
#include <omp.h>
#include <ostream>
#include "../SimStructure/Constants.h"

void applyEmerichBoundary(LBMParams &params) {
#pragma omp parallel for collapse(2)
    for (int j = 0; j < Ny; j++) {
        for (int k : {3,6,7}) {
            if (params.isSolid[index2d(Nx-1, j)]) continue;
            params.f[index3D(k, Nx-1, j)] = (params.rho[index2d(Nx-1,j)]/params.rho[index2d(Nx-2,j)])*params.f[index3D(k, Nx-2, j)];
        }
    }
}

void applySecondOrderExtrapolationBoundary(LBMParams &params) {
#pragma omp parallel for collapse(2)
    for (int j = 0; j < Ny; j++) {
        for (int k : {3,6,7}) {
            if (params.isSolid[index2d(Nx-1, j)]) continue;
            params.f[index3D(k, Nx-1, j)] = 2*params.f[index3D(k, Nx-2, j)] - params.f[index3D(k, Nx-3, j)];
        }
    }
}

void applyConvectiveBoundary(LBMParams &params) {
#pragma omp parallel for collapse(2)
    for (int j = 0; j < Ny; j++) {
        for (int k = 0; k < K; k++) {
            if (params.isSolid[index2d(Nx-2,j)]) continue;
            params.f[index3D(k, Nx - 1, j)] = (params.f_last[index3D(k, Nx - 1, j)] + params.uo * params.f[index3D(k, Nx - 2, j)]) / (1 + params.uo);
        }
    }
}

void applyExitInletBoundary(LBMParams &params) {
#pragma omp parallel for
    for (int j = 1; j < Ny-1; j++) {
        double rhoe = 1.0;
        double vx = -1 + (params.f[index3D(0,Nx-1,j)] + params.f[index3D(2,Nx-1,j)] +
            params.f[index3D(4,Nx-1,j)] + 2*(params.f[index3D(1,Nx-1,j)] + params.f[index3D(5,Nx-1,j)] +
                params.f[index3D(8,Nx-1,j)])) /rhoe;

        params.f[index3D(3,Nx-1,j)] = params.f[index3D(1,Nx-1,j)] - 2.0/3.0 * rhoe * vx;
        params.f[index3D(7,Nx-1,j)] = params.f[index3D(5,Nx-1,j)] + 0.5 * (params.f[index3D(2,Nx-1,j)] -
            params.f[index3D(4,Nx-1,j)]) - 1.0/6.0 * rhoe * vx;
        params.f[index3D(6,Nx-1,j)] = params.f[index3D(8,Nx-1,j)] - 0.5 * (params.f[index3D(2,Nx-1,j)] -
            params.f[index3D(4,Nx-1,j)]) - 1.0/6.0 * rhoe * vx;
    }
}

void applyBoundaryCondition(LBMParams &params, BoundaryConditionType type) {
    switch (type) {
        case BoundaryConditionType::EMERICH:
            applyEmerichBoundary(params);
            break;
        case BoundaryConditionType::SECOND_ORDER_EXTRAPOLATION:
            applySecondOrderExtrapolationBoundary(params);
            break;
        case BoundaryConditionType::CONVECTIVE:
            applyConvectiveBoundary(params);
            break;
        case BoundaryConditionType::EXIT_INLET:
            applyExitInletBoundary(params);
            break;
        default:
            // Handle unknown type or do nothing
            break;
    }
}

std::string boundaryConditionToString(BoundaryConditionType type) {
    switch (type) {
        case BoundaryConditionType::EMERICH:
            return "EMERICH";
        case BoundaryConditionType::SECOND_ORDER_EXTRAPOLATION:
            return "SECOND_ORDER_EXTRAPOLATION";
        case BoundaryConditionType::CONVECTIVE:
            return "CONVECTIVE";
        case BoundaryConditionType::EXIT_INLET:
            return "Extrapolation density outlet";
        default:
            return "UNKNOWN";
    }
}