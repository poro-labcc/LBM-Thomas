#include "./Initialize.h"
#include <cmath>
#include <omp.h>

void Initialize(LBMParams &params, const DomainParams &geometry, bool Parabolic, bool Cube) {
    // Paralelizar a inicialização do cubo
    /*if (Cube) {
        double centerX = 499.5;
        double centerY = 159.5;
        double radius = geometry.CubeD / 2.0;

        #pragma omp parallel for collapse(2)
        for (int i = 0; i < Nx; i++) {
            for (int j = 0; j < Ny; j++) {
                double dx = i - centerX;
                double dy = j - centerY;
                double dist = std::sqrt(dx * dx + dy * dy);

                if (dist <= radius) {
                    params.isSolid[index2d(i, j)] = true;
                }
            }
        }

        //Calculate the distance vector

        double r = geometry.CubeD / 2.0;

        for (int i = 0; i < Nx; i++) {
            for (int j = 0; j < Ny; j++) {

                if (params.isSolid[index2d(i,j)])
                    continue; // Só consideramos nós fluido

                double x0 = i;
                double y0 = j;

                for (int k = 1; k < K; k++) { // direção 0 (parada) não tem vizinho

                    int inb = i + params.cx[k];
                    int jnb = j + params.cy[k];

                    // Verifica se está dentro dos limites do domínio
                    if (inb < 0 || inb >= Nx || jnb < 0 || jnb >= Ny)
                        continue;

                    if (params.isSolid[index2d(inb, jnb)]) {
                        // O vizinho é sólido: calcular a distância até a borda da circunferência

                        // Resolve (x0 + s*dx - xc)^2 + (y0 + s*dy - yc)^2 = R^2
                        double a = params.cx[k]*params.cx[k] + params.cy[k]*params.cy[k];
                        double b = 2 * (params.cx[k]*(x0 - centerX) + params.cy[k]*(y0 - centerY));
                        double c = (x0 - centerX)*(x0 - centerX) + (y0 - centerY)*(y0 - centerY) - r*r;

                        double discriminant = b*b - 4*a*c;

                        if (discriminant >= 0) {
                            double s = (-b - std::sqrt(discriminant)) / (2*a); // menor raiz
                            double distance = s * std::sqrt(a); // distância física

                            params.dist[index3D(k, i, j)] = distance / std::sqrt((i - inb)*(i - inb) + (j - jnb)*(j - jnb)) ; // armazena
                            //params.dist[index3D(k, i, j)] = 0.5 ; // armazena
                        }
                    }
                }
            }
        }
    }*/
    if (Cube) {
        #pragma omp parallel for collapse(2)
        for (int i = geometry.L1 - (geometry.CubeD/2); i < geometry.L1 + (geometry.CubeD/2); i++) {
            for (int j = geometry.middle - (geometry.CubeD/2); j < geometry.middle + (geometry.CubeD/2); j++) {
                params.isSolid[index2d(i, j)] = true;
            }
        }
    }
    // Paralelizar a inicialização das paredes
#pragma omp parallel for
    for (int i = 0; i < Nx; i++) {
        params.isSolid[index2d(i, 0)] = true;
        params.isSolid[index2d(i, Ny - 1)] = true;
        params.dist[index3D(4,i,1)] = 0.5;
        params.dist[index3D(2,i,Ny-2)] = 0.5;
        params.dist[index3D(7,i,1)] = 0.5 * sqrt(2);
        params.dist[index3D(6,i,Ny-2)] = 0.5 * sqrt(2);
        params.dist[index3D(8,i,1)] = 0.5 * sqrt(2);
        params.dist[index3D(5,i,Ny-2)] = 0.5 * sqrt(2);
    }

    //Initializing arrays
    if (Parabolic) {
#pragma omp parallel for collapse(2)
        for (int j = 0; j < Ny; j++) {
            for (int i = 0; i < Nx; i++) {
                if (!params.isSolid[index2d(i, j)]) {
                    double uinit = (-params.uo / 25281.0) * (j - 0.5) * (j - (Ny - 1.5));
                    params.rho[index2d(i, j)] = params.rhoo;
                    params.u[index2d(i, j)] = uinit;
                    params.v[index2d(i, j)] = 0.0;

                    double t1 = params.u[index2d(i, j)] * params.u[index2d(i, j)] + params.v[index2d(i, j)] * params.v[index2d(i, j)];
                    
                    for (int k = 0; k < K; k++) {
                        double t2 = params.u[index2d(i, j)] * params.cx[k] + params.v[index2d(i, j)] * params.cy[k];
                        params.feq[index3D(k, i, j)] =
                                params.rho[index2d(i, j)] * params.w[k] * (1.0 + 3.0 * t2 + 4.5 * t2 * t2 - 1.5 * t1);
                        params.f[index3D(k, i, j)] = params.feq[index3D(k, i, j)];
                    }
                } else {
                    params.rho[index2d(i, j)] = params.rhoo;
                    params.u[index2d(i, j)] = 0.0;
                    params.v[index2d(i, j)] = 0.0;
                }
            }
        }
    } else {
#pragma omp parallel for collapse(2)
        for (int j = 0; j < Ny; j++) {
            for (int i = 0; i < Nx; i++) {
                if (!params.isSolid[index2d(i, j)]) {
                    params.rho[index2d(i, j)] = params.rhoo;
                    params.u[index2d(i, j)] = 0.0;
                    params.v[index2d(i, j)] = 0.0;
                    for (int k = 0; k < K; k++) {
                        params.feq[index3D(k, i, j)] = params.rho[index2d(i, j)] * params.w[k];
                        params.f[index3D(k, i, j)] = params.feq[index3D(k, i, j)];
                    }
                } else {
                    params.rho[index2d(i, j)] = params.rhoo;
                    params.u[index2d(i, j)] = 0.0;
                    params.v[index2d(i, j)] = 0.0;
                }
            }
        }
    }
}
