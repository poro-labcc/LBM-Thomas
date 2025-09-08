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

void RunSimulation(LBMParams &params, SimulationStats &stats, const DomainParams &domainParams,
    BoundaryConditionType bc_type, OutputType print_type, std::optional<int> ReOnly) {
    std::vector<int> Reynolds;

    if (ReOnly.has_value()) {
        Reynolds = {ReOnly.value()};
    } else {
        // Reynolds numbers to simulate (em ordem decrescente)
        Reynolds = {300, 250, 200, 150, 100, 90, 80, 70, 60,50};
        //Reynolds = {90, 80, 70, 60};
    }

    const double uo_final = 0.1 / sqrt(3); //* 0.43;
    //const double uo_final = 0.75;
    params.uo = uo_final;
    params.rhoo = 1.5;

    const int steps_per_Re = 350000;
    //const int ramp_steps = 5000;

    Initialize(params, domainParams, false, true);

    for (size_t re_idx = 0; re_idx < Reynolds.size(); re_idx++) {
        int current_Re = Reynolds[re_idx];

        const double alpha = uo_final * domainParams.CubeD / current_Re;
        const double tau = (3. * alpha + 0.5);
        params.omega = 1.0 / tau;

        std::cout << "\nStarting simulation for Re = " << current_Re << std::endl;
        std::cout << "\nBoundary Outlet: " << boundaryConditionToString(bc_type) << std::endl;
        std::cout << "Alpha = " << alpha << "\t\tMach = " << uo_final * sqrt(3) << std::endl;
        std::cout << "Omega = " << params.omega << "\t\tRelaxation Time = " << tau << std::endl;
        std::cout << "MaxVelo = " << params.uo  << std::endl;

        int mstep = 0;

        #pragma omp parallel for
        for (size_t i = 0; i < params.f.size(); i++) {
            params.f_last[i] = params.f[i];
        }

        while (mstep < steps_per_Re) {
            /*if (mstep < ramp_steps) {
                double ramp_factor = (std::exp(static_cast<double>(mstep) / ramp_steps) - 1.0)
                                     / (std::exp(1.0) - 1.0); // 0 → 1
                params.uo = uo_final * ramp_factor;
            } else {
                params.uo = uo_final;
            }*/

            #pragma omp parallel for
            for (size_t i = 0; i < params.f.size(); i++) {
                params.f_last[i] = params.f[i];
            }

            CollisionBGK(params);
            if (mstep % 10 == 0) ComputeForces(params, stats, domainParams);
            Streaming(params, stats);
            Boundary(params, bc_type);
            MacroRecover(params);

            if (mstep % 10 == 0) {
                stats.calculateCoefficients(params, domainParams);
                saveForce(stats, mstep, 10, current_Re);

                if (mstep % 10 == 0) {
                    stats.calculateMassFlow(params);
                    stats.computeRelativeDifference(params);

                    std::cout << "Re:" << current_Re << " Step:" << mstep
                              << "\tDensity-> " << params.rho[index2d(Nx / 2, Ny / 2)]
                              << "\tVel-> " << params.u[index2d(Nx / 2, Ny / 2)] << "\t" << params.v[index2d(Nx / 2, Ny / 2)]
                              << "\tDelta Mass: " << stats.inletMassFlow - stats.outletMassFlow
                              << "\tDeriv. " << stats.relativeDifference
                              << "\tInflow " << stats.inletMassFlow
                              << "\tCd: " << stats.Cd << " Cl: " << stats.Cl << std::endl;
                    SaveFlow(mstep + re_idx * steps_per_Re,stats);
                }
            }

            if (mstep % 1000 == 0) {
                switch (print_type) {
                    case OutputType::All:
                        SaveVTK(mstep + re_idx * steps_per_Re, params);
                        break;
                    case OutputType::OnlyInOut:
                        SaveFlow(mstep + re_idx * steps_per_Re,stats);
                        break;
                    case OutputType::None:
                        break;
                }
            }

            stats.reset();
            mstep++;
        }
    }
    SaveVTK(steps_per_Re,params);
}

