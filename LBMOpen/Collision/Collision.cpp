#include "./Collision.h"
#include <omp.h>
#include "../SimStructure/Constants.h"


// Matrizes de transformação MRT (D2Q9)
const double M[9][9] = {
    {1, 1, 1, 1, 1, 1, 1, 1, 1},
    {-4, -1, -1, -1, -1, 2, 2, 2, 2},
    {4, -2, -2, -2, -2, 1, 1, 1, 1},
    {0, 1, 0, -1, 0, 1, -1, -1, 1},
    {0, -2, 0, 2, 0, 1, -1, -1, 1},
    {0, 0, 1, 0, -1, 1, 1, -1, -1},
    {0, 0, -2, 0, 2, 1, 1, -1, -1},
    {0, 1, -1, 1, -1, 0, 0, 0, 0},
    {0, 0, 0, 0, 0, 1, -1, 1, -1}
};

const double M_inv[9][9] = {
    {1.0/9, -1.0/9, 1.0/9, 0, 0, 0, 0, 0, 0},
    {1.0/9, -1.0/36, -1.0/18, 1.0/6, -1.0/6, 0, 0, 1.0/4, 0},
    {1.0/9, -1.0/36, -1.0/18, 0, 0, 1.0/6, -1.0/6, -1.0/4, 0},
    {1.0/9, -1.0/36, -1.0/18, -1.0/6, 1.0/6, 0, 0, 1.0/4, 0},
    {1.0/9, -1.0/36, -1.0/18, 0, 0, -1.0/6, 1.0/6, -1.0/4, 0},
    {1.0/9, 1.0/18, 1.0/36, 1.0/6, 1.0/12, 1.0/6, 1.0/12, 0, 1.0/4},
    {1.0/9, 1.0/18, 1.0/36, -1.0/6, -1.0/12, 1.0/6, 1.0/12, 0, -1.0/4},
    {1.0/9, 1.0/18, 1.0/36, -1.0/6, -1.0/12, -1.0/6, -1.0/12, 0, 1.0/4},
    {1.0/9, 1.0/18, 1.0/36, 1.0/6, 1.0/12, -1.0/6, -1.0/12, 0, -1.0/4}
};

// Vetor de relaxação (S): ajustável
void compute_S_vector(double omega, double S[9]) {
    S[0] = 0.0;       // Conservado (ρ)
    S[1] = 1.64;      // e
    S[2] = 1.54;      // ε
    S[3] = 0.0;       // jx
    S[4] = 1.0;       // qx
    S[5] = 0.0;       // jy
    S[6] = 1.0;       // qy
    S[7] = omega;     // pxx
    S[8] = omega;     // pxy
}

// Multiplicação M * f → m
void matVecMult(const double M[9][9], const double in[9], double out[9]) {
    for (int i = 0; i < 9; i++) {
        out[i] = 0.0;
        for (int j = 0; j < 9; j++) {
            out[i] += M[i][j] * in[j];
        }
    }
}

// Momentos de equilíbrio
void computeMeq(double rho, double ux, double uy, double meq[9]) {
    double u2 = ux*ux + uy*uy;
    meq[0] = rho;
    meq[1] = -2*rho + 3*rho*u2;
    meq[2] = rho - 3*rho*u2;
    meq[3] = rho * ux;
    meq[4] = -rho * ux;
    meq[5] = rho * uy;
    meq[6] = -rho * uy;
    meq[7] = rho * (ux*ux - uy*uy);
    meq[8] = rho * ux * uy;
}

void CollisionBGK(LBMParams &params) {
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

void CollisionMRT(LBMParams &params) {
    double S[9];
    compute_S_vector(params.omega, S);

#pragma omp parallel for collapse(2)
    for (int j = 0; j < Ny; j++) {
        for (int i = 0; i < Nx; i++) {
            if (!params.isSolid[index2d(i, j)]) {
                double ux = params.u[index2d(i, j)];
                double uy = params.v[index2d(i, j)];
                double rho = params.rho[index2d(i, j)];

                double f[9], m[9], meq[9];
                for (int k = 0; k < 9; k++)
                    f[k] = params.f[index3D(k, i, j)];

                matVecMult(M, f, m);
                computeMeq(rho, ux, uy, meq);

                for (int k = 0; k < 9; k++)
                    m[k] = m[k] - S[k] * (m[k] - meq[k]);

                matVecMult(M_inv, m, f);

                for (int k = 0; k < 9; k++)
                    params.f[index3D(k, i, j)] = f[k];
            }
        }
    }
}