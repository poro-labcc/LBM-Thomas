#include "LBMParams.h"

LBMParams::LBMParams(int nx, int ny, int k) : Nx(nx), Ny(ny), K(k),
                                              f(K * Nx * Ny, 0.0), feq(K * Nx * Ny, 0.0), f_last(K * Nx * Ny, 0.0),
                                              f_temp(K * Nx * Ny, 0.0), rho(Nx * Ny, 0.0), u(Nx * Ny, 0.0), v(Nx * Ny, 0.0),
                                              cx(K), cy(K), w(K), omega(), uo(), rhoo(), isSolid(Nx * Ny, false),
                                              stableFlow(false), activation(false), counter(0) , dist(Nx*Ny*K){
    cx = {0.0, 1.0, 0.0, -1.0, 0.0, 1.0, -1.0, -1.0, 1.0};
    cy = {0.0, 0.0, 1.0, 0.0, -1.0, 1.0, 1.0, -1.0, -1.0};
    w[0] = 4.0 / 9.0;
    for (int i = 1; i < 5; i++) w[i] = 1.0 / 9.0;
    for (int i = 5; i < 9; i++) w[i] = 1.0 / 36.0;
}

