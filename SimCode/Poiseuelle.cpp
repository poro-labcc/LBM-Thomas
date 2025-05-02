#include <iostream>
#include <string>
#include <sstream>
#include <cmath>
#include <omp.h>
#include <fstream>
#include <vector>
#include <thread> // Required for sleep_for
#include <chrono> // Required for duration literals
#include <iomanip>
#include <algorithm> // para std::max_element
#include <cassert>
#include <optional> // para std::optional


constexpr int CubeD = 40;
constexpr int Nx = 50*CubeD;
constexpr int Ny = 8*CubeD;
constexpr int L1 = 12.5*CubeD;
constexpr int K = 9;

inline int index2d(int i, int j) {
    return i * Ny + j;
}

inline int index3D(int k, int i, int j) {
    return k * Nx * Ny + i * Ny + j;
}

inline int oppositeDirection(int k) {
    static const int opp[9] = {0, 3, 4, 1, 2, 7, 8, 5, 6};
    return opp[k];
}

struct LBMParams {
    int Nx, Ny, K, counter;
    double omega, uo, rhoo;
    std::vector<double> f, feq, f_last, f_temp, rho, u, v;
    std::vector<double> cx, cy, w;
    std::vector<bool> isSolid;
    bool stableFlow, activation;

    LBMParams(int nx, int ny, int k) : Nx(nx), Ny(ny), K(k),
                                       f(K * Nx * Ny, 0.0), feq(K * Nx * Ny, 0.0), f_last(K * Nx * Ny, 0.0),
                                       f_temp(K * Nx * Ny, 0.0), rho(Nx * Ny, 0.0), u(Nx * Ny, 0.0), v(Nx * Ny, 0.0),
                                       cx(K), cy(K), w(K), omega(), uo(), rhoo(), isSolid(Nx * Ny, false),
                                       stableFlow(false), activation(false), counter(0) {
        cx = {0.0, 1.0, 0.0, -1.0, 0.0, 1.0, -1.0, -1.0, 1.0};
        cy = {0.0, 0.0, 1.0, 0.0, -1.0, 1.0, 1.0, -1.0, -1.0};
        w[0] = 4.0 / 9.0;
        for (int i = 1; i < 5; i++) w[i] = 1.0 / 9.0;
        for (int i = 5; i < 9; i++) w[i] = 1.0 / 36.0;
    }
};

struct DomainParams {
    int CubeD, L1, middle;

    DomainParams(int cubeD, int l1, int ny)
        : CubeD(cubeD), L1(l1), middle((ny) / 2) {
    }
};

struct SimulationStats {
    double Fx, Fy; // Forces in x and y directions
    double Cd, Cl; // Drag and lift coefficients
    double inletMassFlow; // Mass flow at the inlet
    double outletMassFlow; // Mass flow at the outlet
    double relativeDifference; // Relative difference in the distribution functions

    // Constructor to initialize members
    SimulationStats() : Fx(0), Fy(0), Cd(0.0), Cl(0.0), inletMassFlow(0.0), outletMassFlow(0.0),
                        relativeDifference(1.0) {
    }

    // Method to calculate mass flow at inlet and outlet
    void calculateMassFlow(const LBMParams &params) {
        inletMassFlow = 0.0;
        outletMassFlow = 0.0;

        // Fluxo na entrada (i=0)
#pragma omp parallel for collapse(2) reduction(+:inletMassFlow,outletMassFlow)
        for (int j = 0; j < Ny; j++) {
            for (int k = 0; k < K; k++) {
                inletMassFlow += params.f[index3D(k, 0, j)] * params.cx[k];
                outletMassFlow += params.f[index3D(k, Nx - 1, j)] * params.cx[k];
            }
        }
    }

    // Method to compute the relative difference in distribution functions
    void computeRelativeDifference(const LBMParams &params) {
        double norm_diff = 0.0;
        double norm_old = 0.0;

#pragma omp parallel for reduction(+:norm_diff, norm_old)
        for (size_t i = 0; i < params.f.size(); i++) {
            norm_diff += std::abs(params.f[i] - params.f_last[i]);
            norm_old += std::abs(params.f_last[i]);
        }

        relativeDifference = (norm_old > 0) ? (norm_diff / norm_old) : norm_diff;
    }

    // Method to calculate drag and lift coefficients
    void calculateCoefficients(const LBMParams &params, const DomainParams &dim) {
        double denominator = 0.5 * params.rhoo * params.uo * params.uo * dim.CubeD;
        if (std::abs(denominator) < 1e-10) {
            Cd = 0.0;
            Cl = 0.0;
        } else {
            Cd = Fx / denominator;
            Cl = Fy / denominator;
        }
    }

    // Method to reset statistics
    void reset() {
        Fx = 0.0;
        Fy = 0.0;
        Cd = 0.0;
        Cl = 0.0;
        inletMassFlow = 0.0;
        outletMassFlow = 0.0;
    }
};

void Initialize(LBMParams &params, const DomainParams &geometry, bool Parabolic, bool Cube) {
    // Paralelizar a inicialização do cubo
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
    }

    //Initializing arrays
    if (Parabolic) {
#pragma omp parallel for collapse(2)
        for (int j = 0; j < Ny; j++) {
            for (int i = 0; i < Nx; i++) {
                if (!params.isSolid[index2d(i, j)]) {
                    double uinit = params.uo * (1.0 - std::pow((j-geometry.middle)/geometry.middle,2));
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
                    params.u[index2d(i, j)] = 0;
                    params.v[index2d(i, j)] = 0;
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

void Collision(LBMParams &params) {
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

void Streaming(LBMParams &params, SimulationStats &stats) {
    static const std::vector<int> ops_k = {0, 3, 4, 1, 2, 7, 8, 5, 6};

#pragma omp parallel for collapse(2)
    for (int i = 0; i < Nx; i++) {
        for (int j = 0; j < Ny; j++) {
            if (params.isSolid[index2d(i, j)]) continue;

            for (int k = 0; k < K; k++) {
                int ii = i - params.cx[k];
                int jj = j - params.cy[k];

                // Skip if neighbor is out of bounds
                if (ii < 0 || ii >= Nx || jj < 0 || jj >= Ny) continue;

                if (!params.isSolid[index2d(ii, jj)]) {
                    // Normal streaming: pull from neighbor
                    params.f_temp[index3D(k, i, j)] = params.f[index3D(k, ii, jj)];
                } else {
                    // Bounce-back: reflect from current node
                    params.f_temp[index3D(k, i, j)] = params.f[index3D(ops_k[k], i, j)];
                }
            }
        }
    }

    // Swap buffers
    std::swap(params.f, params.f_temp);
}

void Boundary(LBMParams &params) {
#pragma omp parallel for
    for (int j = 1; j < Ny-1; j++) {
        double vx = (-params.uo / 25281.0) * (j - 0.5) * (j - (Ny - 1.5));
        double vy = 0.0;

        double rhow = 2*params.rho[index2d(1, j)] - params.rho[index2d(2,j)];

        double t1 = vx * vx + vy * vy;
        for (int k : {1, 5, 8}) {
            double t2 = vx * params.cx[k] + vy * params.cy[k];
            params.f[index3D(k, 0, j)] =
                    rhow * params.w[k] * (1.0 + 3.0 * t2 + 4.5 * t2 * t2 - 1.5 * t1);
        }
    }

    //Emerich boundary condition
#pragma omp parallel for collapse(2)
    for (int j = 0; j < Ny; j++) {
        for (int k = 0; k < K; k++) {
            params.f[index3D(k, Nx-1, j)] = (params.rho[index2d(Nx-1,j)]/params.rho[index2d(Nx-2,j)])*params.f[index3D(k, Nx-2, j)];
        }
    }

/*
    //Convective Boundary in the Outlet using Local velocity
#pragma omp parallel for collapse(2)
    for (int j = 0; j < Ny; j++) {
        for (int k = 0; k < K; k++) {
            if (params.isSolid[index2d(Nx-2,j)]) continue;
            params.f[index3D(k, Nx - 1, j)] = (params.f_last[index3D(k, Nx - 1, j)] + params.uo * params.f[index3D(k, Nx - 2, j)]) / (1 + params.uo);
        }
    }*/
}

void MacroRecover(LBMParams &params) {
#pragma omp parallel for collapse(2)
    //Recovering macroscopic properties
    for (int j = 0; j < Ny; j++) {
        for (int i = 0; i < Nx; i++) {
            if (!params.isSolid[index2d(i, j)]) {
                double ssum = 0;
                double usum = 0;
                double vsum = 0;
                for (int k = 0; k < K; k++) {
                    ssum = ssum + params.f[index3D(k, i, j)];
                    usum = usum + params.f[index3D(k, i, j)] * params.cx[k];
                    vsum = vsum + params.f[index3D(k, i, j)] * params.cy[k];
                }
                params.rho[index2d(i, j)] = ssum;
                params.u[index2d(i, j)] = usum / params.rho[index2d(i, j)];
                params.v[index2d(i, j)] = vsum / params.rho[index2d(i, j)];
            }
        }
    }
}

void ComputeForces(const LBMParams &params, SimulationStats &stats, const DomainParams &geometry) {
    double fx = 0.0;
    double fy = 0.0;
    double density = 0.0;
    double count = 0.0;

#pragma omp parallel for collapse(2) reduction(+:fx,fy,density,count)
    for (int i = geometry.L1 - (geometry.CubeD/2) -1 ; i < geometry.L1 + (geometry.CubeD/2) + 1; i++) {
        for (int j = geometry.middle - (geometry.CubeD/2) - 1; j < geometry.middle + (geometry.CubeD/2) + 1; j++){
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

void saveField(LBMParams &params, const int time, int step) {
    if (time % step == 0) {
        std::stringstream n;
        n << time;
        std::ofstream file("velocity_profile" + n.str() + ".txt");
        
        if (!file.is_open()) {
            std::cerr << "Error opening velocity_profile file for writing!" << std::endl;
            return;
        }
        
#pragma omp parallel
        {
            std::ostringstream localBuffer;
            
#pragma omp for collapse(2) ordered
            for (int j = 0; j < Ny; j++) {
                for (int i = 0; i < Nx; i++) {
                    localBuffer << params.u[index2d(i, j)] << "\t" << params.v[index2d(i, j)] << std::endl;
                    
#pragma omp ordered
                    {
                        file << localBuffer.str();
                        localBuffer.str("");
                        localBuffer.clear();
                    }
                }
            }
        }
        file.close();
        
        std::ofstream file2("density_profile" + n.str() + ".txt");
        if (!file2.is_open()) {
            std::cerr << "Error opening density_profile file for writing!" << std::endl;
            return;
        }
        
#pragma omp parallel
        {
            std::ostringstream localBuffer;
            
#pragma omp for collapse(2) ordered
            for (int j = 0; j < Ny; j++) {
                for (int i = 0; i < Nx; i++) {
                    localBuffer << params.rho[index2d(i, j)] << std::endl;
                    
#pragma omp ordered
                    {
                        file2 << localBuffer.str();
                        localBuffer.str("");
                        localBuffer.clear();
                    }
                }
            }
        }
        file2.close();
    }
}

void saveForce(SimulationStats &stats, const int time, int step, int Re) {
    if (time % step == 0) {
        std::ofstream forceFile("force.txt", std::ios::app);
        if (forceFile.is_open()) {
            forceFile << Re << "\t" << stats.Cd << "\t" << stats.Cl << std::endl;
            forceFile.close();
        } else {
            std::cerr << "Error opening force.txt for writing!" << std::endl;
        }
    }
}

void SaveVTK(int timestep, const LBMParams &params) {
    // Format the filename with timestep
    std::ostringstream filename;
    filename << "velocity/field" << std::setw(5) << std::setfill('0') << timestep << ".vtk";

    std::ofstream file(filename.str());
    if (!file) {
        std::cerr << "Error opening file: " << filename.str() << std::endl;
        return;
    }

    // Write the VTK file header
    file << "# vtk DataFile Version 3.0" << std::endl;
    file << "LBM Output" << std::endl;
    file << "ASCII" << std::endl;
    file << "DATASET STRUCTURED_POINTS" << std::endl;
    file << "DIMENSIONS " << params.Nx << " " << params.Ny << " 1" << std::endl;
    file << "ORIGIN 0 0 0" << std::endl;
    file << "SPACING 1 1 1" << std::endl;
    file << "POINT_DATA " << params.Nx * params.Ny << std::endl;

    // Save velocity field as vectors
    file << "VECTORS Velocity float" << std::endl;
    
#pragma omp parallel
    {
        std::ostringstream localBuffer;
        
#pragma omp for collapse(2) ordered
        for (int j = 0; j < params.Ny; j++) {
            for (int i = 0; i < params.Nx; i++) {
                localBuffer << params.u[index2d(i, j)] << " "
                         << params.v[index2d(i, j)] << " 0.0" << std::endl;
                
#pragma omp ordered
                {
                    file << localBuffer.str();
                    localBuffer.str("");
                    localBuffer.clear();
                }
            }
        }
    }

    // Save density field as a scalar field
    file << "SCALARS Density float 1" << std::endl;
    file << "LOOKUP_TABLE default" << std::endl;
    
#pragma omp parallel
    {
        std::ostringstream localBuffer;
        
#pragma omp for collapse(2) ordered
        for (int j = 0; j < params.Ny; j++) {
            for (int i = 0; i < params.Nx; i++) {
                localBuffer << params.rho[index2d(i, j)] << std::endl;
                
#pragma omp ordered
                {
                    file << localBuffer.str();
                    localBuffer.str("");
                    localBuffer.clear();
                }
            }
        }
    }

    file.close();
}

void RunSimulation(LBMParams &params, SimulationStats &stats, const DomainParams &domainParams, std::optional<int> ReOnly = std::nullopt) {
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

    Initialize(params, domainParams, false, true);

    for (size_t re_idx = 0; re_idx < Reynolds.size(); re_idx++) {
        int current_Re = Reynolds[re_idx];

        const double alpha = uo_final * CubeD / current_Re;
        const double tau = (3. * alpha + 0.5);
        params.omega = 1.0 / tau;

        std::cout << "\nStarting simulation for Re = " << current_Re << std::endl;
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
            ComputeForces(params, stats, domainParams);
            Streaming(params, stats);
            Boundary(params);
            MacroRecover(params);

            if (mstep % 100 == 0) {
                stats.calculateCoefficients(params, domainParams);
                saveForce(stats, mstep, 100, current_Re);

                if (mstep % 1000 == 0) {
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

            if (mstep % 10000 == 0) {
                SaveVTK(mstep + re_idx * steps_per_Re, params);
            }

            stats.reset();
            mstep++;
        }
    }
}

int main() {
    // Configurar o número de threads OpenMP
    int num_threads = omp_get_max_threads();
    omp_set_num_threads(num_threads);
    
    std::cout << "####Geometrical Problem####"<< std::endl;
    std::cout << "Domain dimension: "<< Nx << "x"<< Ny << std::endl;
    std::cout << "Cube: "<< CubeD << " Position: "<< L1 << "\n";
    std::cout << "####Simulation specs####"<< std::endl;
    std::cout << "Running with " << num_threads << " OpenMP threads" << std::endl;

    LBMParams params(Nx, Ny, K);
    DomainParams const domainParams(CubeD, L1, Ny);
    SimulationStats stats;


    // Criar diretório para arquivos VTK se não existir
    system("mkdir -p velocity");

    // Clear force file at start
    std::ofstream clearFile("force.txt");
    clearFile.close();

    //Start clock
    auto const start = std::chrono::high_resolution_clock::now();

    //Start simulation
    RunSimulation(params, stats, domainParams,300);

    // Calcular tempo de execução
    auto const end = std::chrono::high_resolution_clock::now();
    auto const duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);
    
    std::cout << "Simulation completed in " << duration.count() / 1000.0 << " seconds" << std::endl;
    
    return 0;
}
