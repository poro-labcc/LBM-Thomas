#include "./BoundaryConditions.h"

#include <iostream>
#include <omp.h>
#include <ostream>
#include "../SimStructure/Constants.h"
#include "../PostProcess/MacroRecover.h"

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

void applyConvectiveBoundaryNardelli(LBMParams &params) {

    double U = 0.0;
#pragma omp parallel for
    for (int j = 0; j < Ny; j++) {
        if (params.isSolid[index2d(Nx-2,j)]) continue;
        double u, v;
        MacroRecoverPoint(params, Nx-2, j, u, v);
        U += u;
    }

    U /= Ny-2;
#pragma omp parallel for collapse(2)
    for (int j = 0; j < Ny; j++) {
        for (int k : {3,6,7}) {
            if (params.isSolid[index2d(Nx-1,j)]) continue;
            params.f[index3D(k, Nx - 1, j)] = (params.f_last[index3D(k, Nx - 1, j)] + U* params.f[index3D(k, Nx - 2, j)]) / (1 + U);
        }
    }
}

void applyConvectiveBoundaryWang(LBMParams &params) {
    #pragma omp parallel for collapse(2)
    for (int j = 0; j < Ny; j++) {
        for (int k : {3,6,7}) {
            if (params.isSolid[index2d(Nx-1,j)]) continue;
            params.f[index3D(k, Nx - 1, j)] = (params.f_last[index3D(k, Nx - 1, j)] + params.uo * params.f[index3D(k, Nx - 2, j)]) / (1.0 + params.uo);
        }
    }
}

void applyExitInletBoundary(LBMParams &params) {
    int i = Nx-1;
#pragma omp parallel for
    for (int j = 1; j < Ny; j++) {

        double vx = params.uo * 4.0 / (H*H) * (j-0.5) * (H - (j-0.5));
        double rhoe =  (params.f[index3D(0,Nx-1,j)] + params.f[index3D(2,Nx-1,j)] +
            params.f[index3D(4,Nx-1,j)] + 2*(params.f[index3D(1,Nx-1,j)] + params.f[index3D(5,Nx-1,j)] +
                params.f[index3D(8,Nx-1,j)])) /(1.0+vx);

        params.f[index3D(3,Nx-1,j)] = params.f[index3D(1,Nx-1,j)] - 2.0/3.0 * rhoe * vx;
        params.f[index3D(6,i,j)] = params.f[index3D(8,i,j)] - 0.5 * (params.f[index3D(2,i,j)] -
            params.f[index3D(4,i,j)]) - (1.0/6.0) * rhoe * vx;
        params.f[index3D(7,i,j)] = params.f[index3D(5,i,j)] + 0.5 * (params.f[index3D(2,i,j)] -
                params.f[index3D(4,i,j)]) - (1.0/6.0) * rhoe * vx;

    }
}

void applyReisBoundary(LBMParams &params) {
    double ux = params.uo;
    //Camada Oeste
    for (int j = 2; j < Ny - 2; j++) {
        double rho = (params.f[index3D(0,0,j)] + params.f[index3D(2,0,j)] + params.f[index3D(4,0,j)]
            + 2*(params.f[index3D(3,0,j)] + params.f[index3D(6,0,j)] + params.f[index3D(7,0,j)]))/(1-ux);

        params.f[index3D(1,0,j)] = params.f[index3D(2,0,j)] + params.f[index3D(3,0,j)] + params.f[index3D(4,0,j)] +
            2*params.f[index3D(6,0,j)] + 2*params.f[index3D(7,0,j)] + rho * ux - rho/3.0 ;
        params.f[index3D(5,0,j)] = -params.f[index3D(2,0,j)] - params.f[index3D(6,0,j)] + rho/6.0 ;
        params.f[index3D(8,0,j)] = -params.f[index3D(4,0,j)] - params.f[index3D(7,0,j)] + rho/6.0 ;
    }

    //Quina Down
    params.f[index3D(1,0,1)] = (3.0*params.f[index3D(0,0,1)]*ux - params.f[index3D(0,0,1)]
        + 3.0*params.f[index3D(3,0,1)]*ux + params.f[index3D(3,0,1)]
        + 4.0*params.f[index3D(4,0,1)]
        + 8.0*params.f[index3D(7,0,1)]) / (3.0 - 3.0*ux);

    params.f[index3D(2,0,1)] = (-3.0*params.f[index3D(0,0,1)]*ux*ux
        + 3.0*params.f[index3D(0,0,1)]*ux - params.f[index3D(0,0,1)] - 6.0*params.f[index3D(3,0,1)]*ux*ux
        + 4.0*params.f[index3D(3,0,1)] - 6.0*params.f[index3D(4,0,1)]*ux*ux + 3.0*params.f[index3D(4,0,1)]*ux
        + params.f[index3D(4,0,1)] - 12.0*params.f[index3D(7,0,1)]*ux*ux + 8.0*params.f[index3D(7,0,1)]) / (3.0 - 3.0*ux);


    params.f[index3D(5,0,1)] = (3.0*params.f[index3D(0,0,1)]*ux*ux - 3.0*params.f[index3D(0,0,1)]*ux
        + 2.0*params.f[index3D(0,0,1)] + 6.0*params.f[index3D(3,0,1)]*ux*ux - 2.0*params.f[index3D(3,0,1)]
        + 6.0*params.f[index3D(4,0,1)]*ux*ux - 2.0*params.f[index3D(4,0,1)] + 12.0*params.f[index3D(7,0,1)]*ux*ux
        + 6.0*params.f[index3D(7,0,1)]*ux - 10.0*params.f[index3D(7,0,1)]) / (6.0 - 6.0*ux);


    params.f[index3D(6,0,1)] = (3.0*params.f[index3D(0,0,1)]*ux*ux - 3.0*params.f[index3D(0,0,1)]*ux
        + params.f[index3D(0,0,1)] + 6.0*params.f[index3D(3,0,1)]*ux*ux - 4.0*params.f[index3D(3,0,1)]
        + 6.0*params.f[index3D(4,0,1)]*ux*ux - 6.0*params.f[index3D(4,0,1)]*ux + 2.0*params.f[index3D(4,0,1)]
        + 12.0*params.f[index3D(7,0,1)]*ux*ux - 6.0*params.f[index3D(7,0,1)]*ux - 2.0*params.f[index3D(7,0,1)]) / (6.0 - 6.0*ux);


    params.f[index3D(8,0,1)] = (params.f[index3D(0,0,1)] + 2.0*params.f[index3D(3,0,1)] + 6.0*params.f[index3D(4,0,1)]*ux
        - 4.0*params.f[index3D(4,0,1)] + 6.0*params.f[index3D(7,0,1)]*ux - 2.0*params.f[index3D(7,0,1)]) / (6.0 - 6.0*ux);


    params.rho[index2d(0,1)] = (params.f[index3D(0,0,1)] + 2.0*params.f[index3D(3,0,1)] + 2.0*params.f[index3D(4,0,1)]
        + 4.0*params.f[index3D(7,0,1)]) / (1.0 - ux);


    //Quina Up
    params.f[index3D(1,0,Ny-2)] = (-3.0*params.f[index3D(0,0,Ny-2)]*ux + params.f[index3D(0,0,Ny-2)]
        - 4.0*params.f[index3D(2,0,Ny-2)] - 3.0*params.f[index3D(3,0,Ny-2)]*ux - params.f[index3D(3,0,Ny-2)]
        - 8.0*params.f[index3D(6,0,Ny-2)]) /(3.0*ux - 3.0);

    params.f[index3D(4,0,Ny-2)] = ( 3.0*params.f[index3D(0,0,Ny-2)]*(ux*ux) - 3.0*params.f[index3D(0,0,Ny-2)]*ux
        + params.f[index3D(0,0,Ny-2)] + 6.0*params.f[index3D(2,0,Ny-2)]*(ux*ux) - 3.0*params.f[index3D(2,0,Ny-2)]*ux
        - params.f[index3D(2,0,Ny-2)] + 6.0*params.f[index3D(3,0,Ny-2)]*(ux*ux) - 4.0*params.f[index3D(3,0,Ny-2)]
        +12.0*params.f[index3D(6,0,Ny-2)]*(ux*ux) - 8.0*params.f[index3D(6,0,Ny-2)])/(3.0*ux - 3.0);

    params.f[index3D(5,0,Ny-2)] = (-     params.f[index3D(0,0,Ny-2)] - 6.0*params.f[index3D(2,0,Ny-2)]*ux
        + 4.0*params.f[index3D(2,0,Ny-2)] - 2.0*params.f[index3D(3,0,Ny-2)] - 6.0*params.f[index3D(6,0,Ny-2)]*ux
        + 2.0*params.f[index3D(6,0,Ny-2)]) /(6.0*ux - 6.0);

    params.f[index3D(7,0,Ny-2)] = (-3.0*params.f[index3D(0,0,Ny-2)]*(ux*ux) + 3.0*params.f[index3D(0,0,Ny-2)]*ux
        - params.f[index3D(0,0,Ny-2)] - 6.0*params.f[index3D(2,0,Ny-2)]*(ux*ux) + 6.0*params.f[index3D(2,0,Ny-2)]*ux
        - 2.0*params.f[index3D(2,0,Ny-2)] - 6.0*params.f[index3D(3,0,Ny-2)]*(ux*ux) + 4.0*params.f[index3D(3,0,Ny-2)]
        -12.0*params.f[index3D(6,0,Ny-2)]*(ux*ux) + 6.0*params.f[index3D(6,0,Ny-2)]*ux + 2.0*params.f[index3D(6,0,Ny-2)]) /(6.0*ux - 6.0);


    params.f[index3D(8,0,Ny-2)] = (-3.0*params.f[index3D(0,0,Ny-2)]*(ux*ux) + 3.0*params.f[index3D(0,0,Ny-2)]*ux
        - 2.0*params.f[index3D(0,0,Ny-2)]
        - 6.0*params.f[index3D(2,0,Ny-2)]*(ux*ux)+ 2.0*params.f[index3D(2,0,Ny-2)]
        - 6.0*params.f[index3D(3,0,Ny-2)]*(ux*ux) + 2.0*params.f[index3D(3,0,Ny-2)] -12.0*params.f[index3D(6,0,Ny-2)]*(ux*ux)
        - 6.0*params.f[index3D(6,0,Ny-2)]*ux +10.0*params.f[index3D(6,0,Ny-2)]) /(6.0*ux - 6.0);

    params.rho[index2d(0,Ny-2)] = (- params.f[index3D(0,0,Ny-2)] - 2.0*params.f[index3D(2,0,Ny-2)]
        - 2.0*params.f[index3D(3,0,Ny-2)] - 4.0*params.f[index3D(6,0,Ny-2)]) /(ux - 1.0);

}

void applyZouHeCondition(LBMParams &params) {
    int i = 0;
    double rhow = 0.0;
#pragma omp parallel for
    for (int j = 1; j < Ny; j++) {
        //double rhow = 1.01;
        //double vx = params.uo;
        double vx = params.uo * 4.0 / (H*H) * (j-0.5) * (H - (j-0.5));

        rhow =  (params.f[index3D(0,i,j)] + params.f[index3D(2,i,j)] + params.f[index3D(4,i,j)]
            + 2*(params.f[index3D(3,i,j)] + params.f[index3D(6,i,j)] + params.f[index3D(7,i,j)]))/(1.0-vx);

        params.f[index3D(1,i,j)] = params.f[index3D(3,i,j)] + (2.0/3.0) * rhow * vx;
        params.f[index3D(8,i,j)] = params.f[index3D(6,i,j)] + 0.5 * (params.f[index3D(2,i,j)] -
                params.f[index3D(4,i,j)]) + (1.0/6.0) * rhow * vx;
        params.f[index3D(5,i,j)] = params.f[index3D(7,i,j)] - 0.5 * (params.f[index3D(2,i,j)] -
                params.f[index3D(4,i,j)]) + (1.0/6.0) * rhow * vx;

    }
}

void applyBreuerCondition(LBMParams &params) {
    int i = 0;
#pragma omp parallel for
    for (int j = 0; j < Ny; j++) {
        if (params.isSolid[index2d(i,j)]) continue;
        //double vx = params.uo;
        double vx = params.uo * 4.0 / (H*H) * (j-0.5) * (H - (j-0.5));
        double vy = 0.0;

        double rhow = 2* params.rho[index2d(i+1,j)] - params.rho[index2d(i+2,j)];

        double t1 = vx * vx;
        for (int k : {1,5,8}) {
            double t2 = vx * params.cx[k];
            params.f[index3D(k, i, j)] =
                    rhow * params.w[k] * (1.0 + 3.0 * t2 + 4.5 * t2 * t2 - 1.5 * t1);
        }
    }
}

void applyBreuerConditionOUT(LBMParams &params) {
    int i = Nx-1;
#pragma omp parallel for
    for (int j = 0; j < Ny; j++) {
        if (params.isSolid[index2d(i,j)]) continue;
        //double vx = params.uo;
        double vx = params.uo * 4.0 / (H*H) * (j-0.5) * (H - (j-0.5));
        double vy = 0.0;

        double rhow = 2* params.rho[index2d(i-1,j)] - params.rho[index2d(i-2,j)];

        double t1 = vx * vx;
        for (int k : {3,6,7}) {
            double t2 = vx * params.cx[k];
            params.f[index3D(k, i, j)] =
                    rhow * params.w[k] * (1.0 + 3.0 * t2 + 4.5 * t2 * t2 - 1.5 * t1);
        }
    }
}

void applyFirstOrderExtrapolation(LBMParams &params) {
#pragma omp parallel for collapse(2)
    for (int j = 0; j < Ny; j++) {
        for (int k : {3,6,7}) {
            if (params.isSolid[index2d(Nx-1, j)]) continue;
            params.f[index3D(k, Nx-1, j)] = params.f[index3D(k, Nx-2, j)];
        }
    }
}


void applyBoundaryCondition(LBMParams &params, BoundaryConditionType type) {
    switch (type) {
        case BoundaryConditionType::EMERICH:
            applyEmerichBoundary(params);
            break;
        case BoundaryConditionType::FIRST_ORD:
            applyFirstOrderExtrapolation(params);
            break;
        case BoundaryConditionType::SECOND_ORD:
            applySecondOrderExtrapolationBoundary(params);
            break;
        case BoundaryConditionType::CONVECTIVE:
            applyConvectiveBoundaryNardelli(params);
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
        case BoundaryConditionType::FIRST_ORD:
            return "FIRST_ORDER_EXTRAPOLATION";
        case BoundaryConditionType::SECOND_ORD:
            return "SECOND_ORDER_EXTRAPOLATION";
        case BoundaryConditionType::CONVECTIVE:
            return "CONVECTIVE";
        case BoundaryConditionType::EXIT_INLET:
            return "Extrapolation density outlet";
        default:
            return "UNKNOWN";
    }
}
