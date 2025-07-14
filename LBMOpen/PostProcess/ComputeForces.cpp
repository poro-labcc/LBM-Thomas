#include "ComputeForces.h"
#include <omp.h>
#include "../SimStructure/Constants.h"

void ComputeForces(const LBMParams &params, SimulationStats &stats, const DomainParams &geometry) {
    double fx = 0.0;
    double fy = 0.0;
    double density = 0.0;
    double count = 0.0;

#pragma omp parallel for collapse(2) reduction(+:fx,fy,density,count)
    for (int i = geometry.L1 - (geometry.CubeD/2) -2 ; i < geometry.L1 + (geometry.CubeD/2) + 2; i++) {
        for (int j = geometry.middle - (geometry.CubeD/2) - 2; j < geometry.middle + (geometry.CubeD/2) + 2; j++){
            if (params.isSolid[index2d(i, j)]) continue;
            for (int k = 0; k < K; k++) {
                int xx = i + params.cx[k];
                int yy = j + params.cy[k];

                if (params.isSolid[index2d(xx, yy)]) {
                    fx += 2 * params.f[index3D(k, i, j)] * params.cx[k];
                    fy += 2 * params.f[index3D(k, i, j)] * params.cy[k];
                    density += params.rho[index2d(i, j)];
                    count += 1;
                }
            }
        }
    }
    
    if (count > 0) {
        density /= count;
    }
    
    stats.Fx = fx / (density > 0 ? density : 1.0);
    stats.Fy = fy / (density > 0 ? density : 1.0);
}

