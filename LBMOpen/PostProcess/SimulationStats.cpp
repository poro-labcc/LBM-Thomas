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
    outletMassFlow = 0.0;

    // Fluxo na entrada (i=0)
#pragma omp parallel for collapse(2) reduction(+:inletMassFlow,outletMassFlow)
    for (int j = 0; j < Ny; j++) {
        for (int k = 0; k < K; k++) {
            inletMassFlow += params.f[index3D(k, 0, j)] * params.cx[k] + params.f[index3D(k, 0, j)] * params.cy[k];
            outletMassFlow += params.f[index3D(k, Nx - 1, j)] * params.cx[k] + params.f[index3D(k, Nx - 1, j)] * params.cy[k];
        }
    }
}

void SimulationStats::computeRelativeDifference(const LBMParams &params) {
    double norm_diff = 0.0;
    double norm_old = 0.0;

#pragma omp parallel for reduction(+:norm_diff, norm_old)
    for (size_t i = 0; i < params.f.size(); i++) {
        norm_diff += std::abs(params.f[i] - params.f_last[i]);
        norm_old += std::abs(params.f_last[i]);
    }

    relativeDifference = (norm_old > 0) ? (norm_diff / norm_old) : norm_diff;
}

void SimulationStats::calculateCoefficients(const LBMParams &params, const DomainParams &dim) {
    double denominator = 0.5 * params.rhoo * params.uo * params.uo * dim.CubeD;
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

