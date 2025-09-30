#include "../PostProcess/SimulationStats.h"
#include <cmath>
#include <algorithm>
#include <cassert>
#include "../SimStructure/Constants.h"

SimulationStats::SimulationStats() : Fx(0), Fy(0), Cd(0.0), Cl(0.0), inletMassFlow(0.0), outletMassFlow(0.0),
                                     relativeDifference(1.0) {
}

void SimulationStats::calculateMassFlow(const LBMParams &params) {
    inletMassFlow = 0.0;

    // Soma de todas as densidades no domínio
#pragma omp parallel for collapse(2) reduction(+:inletMassFlow)
    for (int j = 0; j < Ny; j++) {
        for (int i = 0; i < Nx; i++) {
            inletMassFlow += params.rho[index2d(i,j)];
        }
    }
}

void SimulationStats::computeRelativeDifference(const LBMParams &params) {
    double norm_diff = 0.0;
    double norm_old  = 0.0;

#pragma omp parallel for reduction(+:norm_diff, norm_old)
    for (size_t i = 0; i < params.u.size(); i++) {
        double df = params.u[i] - params.u_old[i];
        norm_diff += df * df;
        norm_old  += params.u_old[i] * params.u_old[i];
    }

    relativeDifference = std::sqrt(norm_diff / norm_old);
}


void SimulationStats::calculateCoefficients(const LBMParams &params, const DomainParams &dim) {
    double denominator = 0.5 * params.uo * params.uo * dim.CubeD;
    if (std::abs(denominator) < 1e-10) {
        Cd = 0.0;
        Cl = 0.0;
    } else {
        Cd = Fx / denominator;
        Cl = Fy / denominator;
    }
}

void SimulationStats::reset() {
    Fx = 0.0;
    Fy = 0.0;
    Cd = 0.0;
    Cl = 0.0;
    inletMassFlow = 0.0;
    outletMassFlow = 0.0;
}

