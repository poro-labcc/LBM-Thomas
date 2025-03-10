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


constexpr int Nx = 2001;
constexpr int Ny = 321;
constexpr int K = 9;

inline int index2d(int i, int j) {
    return i * Ny + j;
}

inline int index3D(int k, int i, int j) {
    return k * Nx * Ny + i * Ny + j;
}

struct LBMParams {
    int Nx, Ny, K, counter;
    double omega, uo, rhoo;
    std::vector<long double> f, feq, f_last, f_temp, rho, u, v;
    std::vector<double> cx, cy, w;
    std::vector<bool> isSolid;
    bool stableFlow, activation;

    LBMParams(int nx, int ny, int k) : Nx(nx), Ny(ny), K(k),
        f(K * Nx * Ny, 0.0), feq(K * Nx * Ny, 0.0), f_last(K * Nx * Ny, 0.0),
        f_temp(K * Nx * Ny, 0.0), rho(Nx * Ny, 0.0), u(Nx * Ny, 0.0), v(Nx * Ny, 0.0),
        cx(K), cy(K), w(K), omega(), uo(), rhoo(), isSolid(Nx*Ny,false),stableFlow(false), activation(false),counter(0){

        cx = {0.0, 1.0, 0.0, -1.0, 0.0, 1.0, -1.0, -1.0, 1.0};
        cy = {0.0, 0.0, 1.0, 0.0, -1.0, 1.0, 1.0, -1.0, -1.0};
        w[0] = 4.0/9.0;
        for(int i=1; i < 5;i++) w[i] = 1.0/9.0;
        for(int i=5; i < 9;i++) w[i] = 1.0/36.0;
    }
};

struct DomainParams {
    int CubeD, L1, middle;

    DomainParams(int cubeD, int l1, int ny)
        : CubeD(cubeD), L1(l1), middle((ny) / 2) {}
};

struct SimulationStats {
    long double Fx, Fy;          // Forces in x and y directions
    long double Cd, Cl;           // Drag and lift coefficients
    long double inletMassFlow;    // Mass flow at the inlet
    long double outletMassFlow;   // Mass flow at the outlet
    long double relativeDifference; // Relative difference in the distribution functions

    // Constructor to initialize members
    SimulationStats() : Fx(0), Fy(0), Cd(0.0), Cl(0.0), inletMassFlow(0.0), outletMassFlow(0.0), relativeDifference(1.0) {}

    // Method to calculate mass flow at inlet and outlet
    void calculateMassFlow(const LBMParams& params) {
        inletMassFlow = 0.0;
        outletMassFlow = 0.0;

        // Fluxo na entrada (i=0)
        for (int j = 0; j < Ny; j++) {
            for (int k = 0; k < K; k++) {
                inletMassFlow += params.f[index3D(k, 0, j)] * params.cx[k];
            }
        }

        // Fluxo na saída (i=Nx-1)
        for (int j = 0; j < Ny; j++) {
            for (int k = 0; k < K; k++) {
                outletMassFlow += params.f[index3D(k, Nx - 1, j)] * params.cx[k];
            }
        }
    }

    // Method to compute the relative difference in distribution functions
    void computeRelativeDifference(const LBMParams& params) {
        long double norm_diff = 0.0;
        long double norm_old = 0.0;

        #pragma omp parallel for reduction(+:norm_diff, norm_old)
        for (size_t i = 0; i < params.f.size(); i++) {
            norm_diff += std::abs(params.f[i] - params.f_last[i]);
            norm_old += std::abs(params.f_last[i]);
        }

        relativeDifference = (norm_old > 0) ? (norm_diff / norm_old) : norm_diff;
    }

    // Method to calculate drag and lift coefficients
    void calculateCoefficients(const LBMParams& params, const DomainParams& dim) {
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

void Initialize(LBMParams& params, bool Parabolic) {
    for (int i = 480 ; i <= 520; i++) {
        for (int j = 140; j <= 180 ; j++) {
            params.isSolid[index2d(i,j)] = true;
        }
    }
    //Initializing arrays
    if (Parabolic) {
        for (int j = 0; j < Ny; j++) {
            long double uinit = (6.0 * params.uo * (2.0 / 3.0)) / ((Ny-1) * (Ny-1)) * j * ((Ny-1) - j);
            for (int i = 0; i < Nx; i++) {
                if (!params.isSolid[index2d(i,j)]) {
                    params.rho[index2d(i, j)] = params.rhoo;
                    params.u[index2d(i, j)] = uinit;
                    params.v[index2d(i, j)] = 0.0;

                    long double t1 = params.u[index2d(i,j)] * params.u[index2d(i,j)] + params.v[index2d(i,j)] * params.v[index2d(i,j)];
                    for (int k = 0; k < K; k++) {
                        long double t2 = params.u[index2d(i,j)] * params.cx[k] + params.v[index2d(i,j)] * params.cy[k];
                        params.feq[index3D(k, i, j)] = params.rho[index2d(i, j)] * params.w[k] * (1.0 + 3.0 * t2 + 4.5 * t2 * t2 - 1.5 * t1);
                        params.f[index3D(k, i, j)] = params.feq[index3D(k, i, j)];
                    }
                }else {
                    params.rho[index2d(i, j)] = params.rhoo;
                    params.u[index2d(i, j)] = 0.0;
                    params.v[index2d(i, j)] = 0.0;
                }
            }
        }
    }else {
#pragma omp parallel for collapse(2)
        for(int j=0; j < Ny; j++){
            for(int i=0; i < Nx; i++){
                params.rho[index2d(i,j)] = params.rhoo;
                params.u[index2d(i,j)] = 0;
                params.v[index2d(i,j)] = 0;
                for(int k=0; k < K; k++){
                    params.feq[index3D(k,i,j)] = params.rho[index2d(i,j)]*params.w[k];
                    params.f[index3D(k,i,j)] = params.feq[index3D(k,i,j)];
                }
            }
        }
    }
}

void Colision(LBMParams& params) {
    #pragma omp parallel for collapse(2)
    for(int j=0; j < Ny; j++){
        for(int i=0; i < Nx; i++){
            if (!params.isSolid[index2d(i,j)]) {
                long double t1 = params.u[index2d(i,j)] * params.u[index2d(i,j)] + params.v[index2d(i,j)] * params.v[index2d(i,j)];
                for(int k=0; k < K; k++){
                    long double t2 = params.u[index2d(i,j)]*params.cx[k]+params.v[index2d(i,j)]*params.cy[k];
                    params.feq[index3D(k,i,j)] = params.rho[index2d(i,j)]*params.w[k]*(1.0+3.0*t2+4.5*t2*t2-1.5*t1);
                    params.f[index3D(k,i,j)]=params.omega*params.feq[index3D(k,i,j)]+(1.0-params.omega)*params.f[index3D(k,i,j)]; //Relaxation Step
                }
            }
        }
    }
}

void Streaming(LBMParams& params) {
    static const std::vector<int> ops_k = {0,3,4,1,2,7,8,5,6};
#pragma omp parallel for collapse(2)
    for(int i = 0; i < Nx; i++) {
        for(int j = 0; j < Ny; j++) {
            for(int k = 0; k < K; k++) {
                int xx = i + params.cx[k];
                int yy = j + params.cy[k];

                if (xx < 0) xx = Nx-1;
                if (yy < 0) yy = Ny-1;
                if (xx > Nx-1) xx = 0;
                if (yy > Ny-1) yy = 0;

                if (!params.isSolid[index2d(xx,yy)] && !params.isSolid[index2d(i,j)]) {
                    params.f_temp[index3D(k,xx,yy)] = params.f[index3D(k,i,j)];
                }else {
                    params.f_temp[index3D(ops_k[k],i,j)] = params.f[index3D(k,i,j)];
                }
            }
        }
    }
    params.f.swap(params.f_temp);
}

void Boundary(LBMParams& params) {
    for (int j = 0; j < Ny; j++) {
        // Corrected velocity profile (positive flow to the right)
        long double vx =  (6.0 * params.uo * (2.0 / 3.0)) / ((Ny-1) * (Ny-1)) * j * ((Ny-1) - j);
        // Corrected density formula (1 - vx[j])
        long double rhow = (params.f[index3D(0,0,j)] + params.f[index3D(2,0,j)] + params.f[index3D(4,0,j)]
                + 2 * (params.f[index3D(3,0,j)] + params.f[index3D(6,0,j)] + params.f[index3D(7,0,j)])) / (1 - vx);

        // Set incoming distributions (f3, f6, f7)
        params.f[index3D(1,0,j)] = params.f[index3D(3,0,j)] + (2.0/3.0) * rhow * vx;
        params.f[index3D(8,0,j)] = params.f[index3D(6,0,j)] + 0.5*(params.f[index3D(2,0,j)] - params.f[index3D(4,0,j)]) + (1.0/6.0) * rhow * vx;
        params.f[index3D(5,0,j)] = params.f[index3D(7,0,j)] - 0.5*(params.f[index3D(2,0,j)] - params.f[index3D(4,0,j)]) + (1.0/6.0) * rhow * vx;
    }

    for (int i = 0; i < Nx; i++) { //Fullway
        // Bounce back on south boundary
        params.f[index3D(2,i,0)] = params.f[index3D(4,i,0)];
        params.f[index3D(5,i,0)] = params.f[index3D(7,i,0)];
        params.f[index3D(6,i,0)] = params.f[index3D(8,i,0)];

        // Bounce back on north boundary
        params.f[index3D(4,i,Ny-1)] = params.f[index3D(2,i,Ny-1)];
        params.f[index3D(8,i,Ny-1)] = params.f[index3D(6,i,Ny-1)];
        params.f[index3D(7,i,Ny-1)] = params.f[index3D(5,i,Ny-1)];
    }

    //Convective Boundary in the Outlet
    double U = 0.;

#pragma omp parallel for
    for (int j = 0; j < Ny; j++){
        double usum = 0;
        for(int k=0; k < params.K; k++){
            usum = usum + params.f[index3D(k,Nx-2,j)]*params.cx[k];
        }
        U += usum/params.rho[index2d(Nx-2,j)];
    }
    U = U/Ny;

#pragma omp parallel for
    for(int j = 0; j < Ny; j++){
        for(int k = 0;k<K;k++){
            params.f[index3D(k,Nx-1,j)] = (params.f_last[index3D(k,Nx-1,j)] + U*params.f[index3D(k,Nx-2,j)])/(1+U);
        }
    }

}

void CubeBDC(LBMParams& params, const DomainParams& domain) {
    int cube_start = domain.L1 - domain.CubeD/2; //480
    int cube_end = domain.L1 + domain.CubeD/2; //520
    int j_top = domain.middle + domain.CubeD/2; //180
    int j_bottom = domain.middle - domain.CubeD/2; //140

    // Top and bottom edges
    for(int i = cube_start + 1; i < cube_end; i++) { // Exclude corners
        // Top edge
        params.f[index3D(2,i,j_top)] = params.f[index3D(4, i, j_top)];
        params.f[index3D(5,i,j_top)] = params.f[index3D(7, i, j_top)];
        params.f[index3D(6,i,j_top)] = params.f[index3D(8, i, j_top)];

        // Bottom edge
        params.f[index3D(4,i,j_bottom)] = params.f[index3D(2, i, j_bottom)];
        params.f[index3D(8,i,j_bottom)] = params.f[index3D(6, i, j_bottom)];
        params.f[index3D(7,i,j_bottom)] = params.f[index3D(5, i, j_bottom)];
    }

    // Left and right edges
    for(int j = j_bottom + 1; j < j_top; j++) { // Exclude corners
        // Left edge
        params.f[index3D(3,cube_start,j)] = params.f[index3D(1, cube_start, j)];
        params.f[index3D(6,cube_start,j)] = params.f[index3D(8, cube_start, j)];
        params.f[index3D(7,cube_start,j)] = params.f[index3D(5, cube_start, j)];

        // Right edge
        params.f[index3D(1,cube_end,j)] = params.f[index3D(3, cube_end, j)];
        params.f[index3D(8,cube_end,j)] = params.f[index3D(6, cube_end, j)];
        params.f[index3D(5,cube_end,j)] = params.f[index3D(7, cube_end, j)];
    }

    // Handle corners explicitly
    // Top-left (480,180)
    params.f[index3D(6, cube_start, j_top)] = params.f[index3D(8, cube_start, j_top)];
    params.f[index3D(3, cube_start, j_top)] = params.f[index3D(1, cube_start, j_top)];
    params.f[index3D(2, cube_start, j_top)] = params.f[index3D(4, cube_start, j_top)];

    // Top-right (520,180)
    params.f[index3D(8, cube_end, j_top)] = params.f[index3D(6, cube_end, j_top)];
    params.f[index3D(1, cube_end, j_top)] = params.f[index3D(3, cube_end, j_top)];
    params.f[index3D(2, cube_end, j_top)] = params.f[index3D(4, cube_end, j_top)];

    // Bottom-left (480,140)
    params.f[index3D(7, cube_start, j_bottom)] = params.f[index3D(5, cube_start, j_bottom)];
    params.f[index3D(3, cube_start, j_bottom)] = params.f[index3D(1, cube_start, j_bottom)];
    params.f[index3D(4, cube_start, j_bottom)] = params.f[index3D(2, cube_start, j_bottom)];

    // Bottom-right (520,140)
    params.f[index3D(5, cube_end, j_bottom)] = params.f[index3D(7, cube_end, j_bottom)];
    params.f[index3D(1, cube_end, j_bottom)] = params.f[index3D(3, cube_end, j_bottom)];
    params.f[index3D(4, cube_end, j_bottom)] = params.f[index3D(2, cube_end, j_bottom)];
}

void MacroRecover(LBMParams& params, const DomainParams& dim, bool isCube) {
    long double ssum, usum, vsum;
    #pragma omp parallel for collapse(2)
    //Recovering macroscopic properties
    for (int j = 0; j < Ny; j++) {
        for(int i=0; i < Nx; i++) {
            if (!params.isSolid[index2d(i,j)]) {
                ssum = 0;
                usum = 0;
                vsum = 0;
                for(int k=0; k < K; k++){
                    ssum = ssum + params.f[index3D(k,i,j)];
                    usum = usum + params.f[index3D(k,i,j)]*params.cx[k];
                    vsum = vsum + params.f[index3D(k,i,j)]*params.cy[k];
                }
                params.rho[index2d(i,j)] = ssum;
                params.u[index2d(i,j)] = usum/params.rho[index2d(i,j)];
                params.v[index2d(i,j)] = vsum/params.rho[index2d(i,j)];
            }
        }
    }
}

void ComputeForces(LBMParams& params, const DomainParams& domain, SimulationStats& stats) {
    int cube_start = domain.L1 - domain.CubeD/2 - 1;
    int cube_end = domain.L1 + domain.CubeD/2 ;
    int j_top = domain.middle + domain.CubeD/2;
    int j_bottom = domain.middle - domain.CubeD/2 -1  ;

    stats.Fx = 0.0;
    stats.Fy = 0.0;

    // =================================================
    // Forças nas bordas superiores e inferiores
    // =================================================
    for(int i = 480; i <= 520 ; i++) {
        // Borda superior (j_top)
        long double f4 = params.f[index3D(4, i, 180)];
        long double f7 = params.f[index3D(7, i, 180)];
        long double f8 = params.f[index3D(8, i, 180)];
        stats.Fx += 2.0 * (f7*(-1.0) + f8*(1.0));  // ΔFx = 2*(f7*(-cx[7]) + f8*(cx[8]))
        stats.Fy += 2.0 * (f4*(-1.0) + f7*(-1.0) + f8*(-1.0)); // ΔFy = 2*(f4*cy[4] + f7*cy[7] + f8*cy[8])

        // Borda inferior (j_bottom)
        long double f2 = params.f[index3D(2, i, 140)];
        long double f5 = params.f[index3D(5, i, 140)];
        long double f6 = params.f[index3D(6, i, 140)];
        stats.Fx += 2.0 * (f6*(-1.0) + f5*(1.0));  // ΔFx = 2*(f6*(-cx[6]) + f5*(cx[5]))
        stats.Fy += 2.0 * (f2*(1.0) + f6*(1.0) + f5*(1.0)); // ΔFy = 2*(f2*cy[2] + f6*cy[6] + f5*cy[5])
    }

    // =================================================
    // Forças nas bordas laterais (esquerda e direita)
    // =================================================
    for(int j = 140; j <= 180; j++) {
        // Borda esquerda (cube_start)
        long double f1 = params.f[index3D(1, 480, j)];
        long double f5 = params.f[index3D(5, 480, j)];
        long double f8 = params.f[index3D(8, 480, j)];
        stats.Fx += 2.0 * (f1*(1.0) + f5*(1.0) + f8*(1.0)); // ΔFx = 2*(f1*cx[1] + f5*cx[5] + f8*cx[8])
        stats.Fy += 2.0 * (f8*(-1.0) + f5*(1.0)); // ΔFy = 2*(f8*cy[8] + f5*cy[5])

        // Borda direita (cube_end)
        long double f3 = params.f[index3D(3, 520, j)];
        long double f6 = params.f[index3D(6, 520, j)];
        long double f7 = params.f[index3D(7, 520, j)];
        stats.Fx += 2.0 * (f3*(-1.0) + f6*(-1.0) + f7*(-1.0)); // ΔFx = 2*(f3*cx[3] + f6*cx[6] + f7*cx[7])
        stats.Fy += 2.0 * (f6*(1.0) + f7*(-1.0)); // ΔFy = 2*(f6*cy[6] + f7*cy[7])
    }

    // =================================================
    // Forças nas quinas do cubo (contribuição adicional)
    // =================================================
    // Quina Inferior-Esquerda (cube_start, j_bottom)
    long double f5_corner = params.f[index3D(5,480,140)];
    stats.Fx += 2.0 * f5_corner * 1.0;  // cx[5] = 1
    stats.Fy += 2.0 * f5_corner * 1.0;  // cy[5] = 1

    // Quina Superior-Esquerda (cube_start, j_top)
    long double f8_corner = params.f[index3D(8,480,180)];
    stats.Fx += 2.0 * f8_corner * 1.0;   // cx[8] = 1
    stats.Fy += 2.0 * f8_corner * (-1.0); // cy[8] = -1

    // Quina Superior-Direita (cube_end, j_top)
    long double f7_corner = params.f[index3D(7,520,180)];
    stats.Fx += 2.0 * f7_corner * (-1.0); // cx[7] = -1
    stats.Fy += 2.0 * f7_corner * (-1.0); // cy[7] = -1

    // Quina Inferior-Direita (cube_end, j_bottom)
    long double f6_corner = params.f[index3D(6,520,140)];
    stats.Fx += 2.0 * f6_corner * (-1.0); // cx[6] = -1
    stats.Fy += 2.0 * f6_corner * 1.0;    // cy[6] = 1
}

void saveField(LBMParams& params, const int time, int step) {
    if(time % step == 0){
        std::stringstream n;
        n << time;
        std::ofstream file("velocity_profile"+n.str()+".txt");
        for(int j = 0; j < Ny; j++){
            for(int i = 0; i < Nx; i++) {
                file << params.u[index2d(i,j)] << "\t"<< params.v[index2d(i,j)] << std::endl;
            }
        }
        file.close();
        std::ofstream file2("density_profile"+n.str()+".txt");
        for(int j = 0; j < Ny; j++){
            for(int i = 0; i < Nx; i++) {
                file2 << params.rho[index2d(i,j)] << std::endl;
            }
        }
        file2.close();
    }
}

void saveForce(SimulationStats& stats, const int time, int step){
    if(time % step == 0) {
        std::ofstream forceFile("force.txt", std::ios::app);
        if (forceFile.is_open()) {
            forceFile << stats.Cd << "\t" << stats.Cl << std::endl;
            forceFile.close();
        } else {
            std::cerr << "Error opening force.txt for writing!" << std::endl;
        }
    }
}

void Pertubation(LBMParams& params) {
    for (int i = 0; i < Nx; i++) {
        for (int j = 0; j < Ny; j++) {
            if (!params.isSolid[index2d(i,j)]) {
                params.v[index2d(i,j)] = params.uo * 0.05 * sin(2.0 * M_PI * ((i) / (Nx-1) + 1.0 / 4.0));
            }
        }
    }
}

void SaveVTK(int timestep, const LBMParams& params) {

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
    for (int j = 0; j < params.Ny; j++) {
        for (int i = 0; i < params.Nx; i++) {
            file << params.u[index2d(i, j)] << " "
                 << params.v[index2d(i, j)] << " 0.0" << std::endl;
        }
    }

    file.close();
}

int main() {
    constexpr int CubeD = 40, L1 = 500;

    constexpr int ReSet = 60;
    constexpr double alpha = 0.1;

    LBMParams params(Nx, Ny, K);
    params.uo = ReSet * alpha / CubeD;
    params.rhoo = 1.0;
    params.omega = 1.0/(3.*alpha+0.5);

    DomainParams const domainParams(CubeD, L1, Ny);


    SimulationStats stats;

    std::cout << "Re = " << ReSet << " ; Max. Velocity = " << params.uo << std::endl;
    std::cout << "Omega = " << params.omega  << std::endl;



    Initialize(params,true);
    int mstep = 0;
    //Start clock
    auto const start = std::chrono::high_resolution_clock::now();

    while (mstep < 100001) {
        params.f_last = params.f;

        if (stats.relativeDifference < 1e-5 && !params.activation) {
            Pertubation(params);
            params.activation = true;
        }
        Colision(params);
        Streaming(params);
        Boundary(params);
        //ComputeForces(params,domainParams,stats);
        MacroRecover(params,domainParams,true);

        stats.calculateMassFlow(params);
        stats.computeRelativeDifference(params);
        stats.calculateCoefficients(params,domainParams);

        //saveField(params,mstep,500);
        if (mstep % 500 == 0) {
            SaveVTK(mstep, params);
        }
        //saveForce(stats,mstep,500);

        std::cout << "Step:" << mstep <<"\t\tDensity-> "<<params.rho[index2d(Nx/2,Ny/2)]<<
            "\t\tVel. Vectors -> " << params.u[index2d(Nx/2,Ny/2)] << "\t" << params.v[index2d(Nx/2,Ny/2)] << "\t" <<
                "Delta Mass in/out: " << stats.inletMassFlow - stats.outletMassFlow << "\t" << "Deriv. Field " << stats.relativeDifference
            << std::endl;
        std::cout  << "<<<<<<<<<<<< Cd: " << stats.Cd << " Cl: " << stats.Cl << " >>>>>>>>>>>>"<< std::endl;

        stats.reset();
        mstep += 1;
    }
    auto const end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> const duration = end - start;
    std::cout << "Execution time: " << duration.count() << " seconds" << std::endl;
    std::cout << "END OF SIMULATION!" << std::endl;
    return 0;
}
