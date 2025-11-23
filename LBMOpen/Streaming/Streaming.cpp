#include "./Streaming.h"
#include <omp.h>
#include <algorithm>
#include <iostream>
#include <bits/ostream.tcc>

#include "../SimStructure/Constants.h"

void Streaming(LBMParams &params, SimulationStats &stats) {
    static const std::vector<int> ops_k = {0, 3, 4, 1, 2, 7, 8, 5, 6};

#pragma omp parallel for collapse(2)
    for (int i = 0; i < Nx; i++) {
        for (int j = 0; j < Ny; j++) {
            if (params.isSolid[index2d(i, j)]) continue;

            for (int k = 0; k < K; k++) {
                int ii = i + params.cx[k];
                int jj = j + params.cy[k];
                // Verifica se está tentando puxar de um nó sólido ou fora do domínio
                if (ii < 0 || ii >= Nx || jj < 0 || jj >= Ny ) continue;

                if (!params.isSolid[index2d(ii, jj)] ) {
                    // Streaming normal: push to neighbor
                    params.f_temp[index3D(k, ii, jj)] = params.f[index3D(k, i, j)];
                } else{
                    //BounceBack Padrão
                    params.f_temp[index3D(ops_k[k], i, j)] = params.f[index3D(k, i, j)];
                    /*BounceBack Bouzini
                    double d = params.dist[index3D(k, i, j)];

                    // Para obstáculos internos (usar distância calculada)
                    if (d < 0.5) {
                        //std::cout << "To aqui mas n " << i<< "  " << j << "  " << d<< "  " << k<<std::endl;
                        params.f_temp[index3D(ops_k[k], i, j)] =
                            2.0 * d * params.f[index3D(k, i, j)] +
                            (1.0 - 2.0 * d) * params.f[index3D(k, i + params.cx[k], j + params.cy[k])];
                    } else {
                        //std::cout << "To aqui" << std::endl;
                        params.f_temp[index3D(ops_k[k], i, j)] =
                            (1.0 / (2.0 * d)) * params.f[index3D(k, i, j)] +
                            ((2.0 * d - 1.0) / (2.0 * d)) * params.f[index3D(k, i, j)];
                    }
                    */
                }
            }
        }
    }

    // Swap buffers
    std::swap(params.f, params.f_temp);
    std::fill(params.f_temp.begin(), params.f_temp.end(), 0.0);
}

