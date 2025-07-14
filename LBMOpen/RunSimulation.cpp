#include "RunSimulation.h"
#include <iostream>
#include <vector>
#include <cmath>
#include <chrono>
#include <omp.h>
#include "./SimStructure/Constants.h"
#include "./SimStructure/LBMParams.h"
#include "SimStructure/DomainParams.h"
#include "./PostProcess/SimulationStats.h"
#include "./SimStructure/Initialize.h"
#include "./Collision/Collision.h"
#include "./Streaming/Streaming.h"
#include "./Boundary/Boundary.h"
#include "./PostProcess/MacroRecover.h"
#include "./PostProcess/ComputeForces.h"
#include "./PostProcess/SaveFunctions.h"
#include "./Boundary/BoundaryConditions.h"

void RunSimulation(LBMParams &params, SimulationStats &stats, const DomainParams &domainParams, BoundaryConditionType bc_type, std::optional<int> ReOnly) {
    std::vector<int> Reynolds;

    if (ReOnly.has_value()) {
        Reynolds = {ReOnly.value()};
    } else {
        // Reynolds numbers to simulate (em ordem decrescente)
        Reynolds = {300, 250, 200, 150, 100, 90, 80, 70};
    }

    const double uo_final = 0.1 / sqrt(3);
    params.uo = uo_final;
    params.rhoo = 1.0;

    const int steps_per_Re = 350000;
    const int ramp_steps = 3000;

    Initialize(params, domainParams, false, false);

    for (size_t re_idx = 0; re_idx < Reynolds.size(); re_idx++) {
        int current_Re = Reynolds[re_idx];

        const double alpha = uo_final * 320 / current_Re;
        const double tau = (3. * alpha + 0.5);
        params.omega = 1.0 / tau;

        std::cout << "\nStarting simulation for Re = " << current_Re << std::endl;
        std::cout << "\nBoundary Outlet: " << boundaryConditionToString(bc_type) << std::endl;
        std::cout << "Alpha = " << alpha << "\t\tMach = " << uo_final * sqrt(3) << std::endl;
        std::cout << "Omega = " << params.omega << "\t\tRelaxation Time = " << tau << std::endl;

        int mstep = 0;

        #pragma omp parallel for
        for (size_t i = 0; i < params.f.size(); i++) {
            params.f_last[i] = params.f[i];
        }

        while (mstep < steps_per_Re) {
            if (re_idx == 0 && mstep < ramp_steps) {
                double ramp_factor = 0.0001 + (0.9999 * mstep) / ramp_steps;
                params.uo = uo_final * ramp_factor;
            } else {
                params.uo = uo_final;
            }

            #pragma omp parallel for
            for (size_t i = 0; i < params.f.size(); i++) {
                params.f_last[i] = params.f[i];
            }

            Collision(params);
            if (mstep % 10 == 0) ComputeForces(params, stats, domainParams);
            Streaming(params, stats);
            Boundary(params, bc_type);
            MacroRecover(params);

            if (mstep % 10 == 0) {
                stats.calculateCoefficients(params, domainParams);
                saveForce(stats, mstep, 10, current_Re);

                if (mstep % 100 == 0) {
                    stats.calculateMassFlow(params);
                    stats.computeRelativeDifference(params);

                    std::cout << "Re:" << current_Re << " Step:" << mstep
                              << "\tDensity-> " << params.rho[index2d(Nx / 2, Ny / 2)]
                              << "\tVel-> " << params.u[index2d(Nx / 2, Ny / 2)] << "\t" << params.v[index2d(Nx / 2, Ny / 2)]
                              << "\tDelta Mass: " << stats.inletMassFlow - stats.outletMassFlow
                              << "\tDeriv. " << stats.relativeDifference
                              << "\tCd: " << stats.Cd << " Cl: " << stats.Cl << std::endl;
                }
            }

            if (mstep % 100 == 0) {
                SaveVTK(mstep + re_idx * steps_per_Re, params);
            }

            stats.reset();
            mstep++;
        }
    }
}

