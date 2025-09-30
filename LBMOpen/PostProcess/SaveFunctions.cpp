#include "./SaveFunctions.h"
#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <fstream>
#include <string>
#include <omp.h>
#include "../SimStructure/Constants.h"

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

void SaveVTK(int timestep, int Re, const LBMParams &params) {
    // Format the filename with timestep
    std::ostringstream filename;
    filename << "velocity/"<< Re <<"field" << std::setw(5) << std::setfill('0') << timestep << ".vtk";

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

void SaveFlow(int timestep, const SimulationStats &stats) {
    std::ofstream file("flow.txt", std::ios::app); // append mode
    if (!file) {
        std::cerr << "Error opening flow.txt" << std::endl;
        return;
    }

    // Save the mass flow difference with timestep
    file << timestep << " " << stats.inletMassFlow << std::endl;

    file.close();
}



// Função para salvar todos os f_k
void saveDistribution(int timestep, const LBMParams& params) {
    std::string filename = "distributions_" + std::to_string(timestep) + ".vtk";
    std::ofstream file(filename);

    int nx = params.Nx;
    int ny = params.Ny;

    file << "# vtk DataFile Version 3.0\n";
    file << "LBM distributions f_k\n";
    file << "ASCII\n";
    file << "DATASET STRUCTURED_POINTS\n";
    file << "DIMENSIONS " << nx << " " << ny << " 1\n";
    file << "ORIGIN 0 0 0\n";
    file << "SPACING 1 1 1\n";
    file << "POINT_DATA " << nx * ny << "\n";

    // Para cada k, salva como SCALARS separado
    for (int k = 0; k < 9; ++k) {
        file << "SCALARS f" << k << " float 1\n";
        file << "LOOKUP_TABLE default\n";

        for (int j = 0; j < ny; ++j) {
            for (int i = 0; i < nx; ++i) {
                file << params.f[index3D(k,i,j)] << "\n";
            }
        }
    }

    file.close();
}
void salvarFAppend(LBMParams &params, const std::string &filename) {
    static bool header_written = false;
    std::ofstream file;

    int i = 0;
    int j = 1;
    int K = 9; // D2Q9

    // Abre arquivo em modo append
    file.open(filename, std::ios::app);

    // Escrever header apenas na primeira chamada
    if (!header_written) {
        for (int k = 0; k < K; k++) {
            file << "f" << k;
            if (k != K - 1) file << ",";
        }
        file << "\n";
        header_written = true;
    }

    // Salva o vetor f do ponto (0,1)
    for (int k = 0; k < K; k++) {
        file << std::setprecision(10) << params.f[index3D(k, i, j)];
        if (k != K - 1) file << ",";
    }
    file << "\n";

    file.close();
}

void saveCentralLineFneq(const LBMParams &params, const std::string &filename, int timestep) {

    std::ofstream fout(filename, std::ios::app); // usa "app" para acumular no mesmo arquivo
    fout << timestep;
    for (int i = 0; i < Nx; i++) {
        fout << " " << params.f_neq[index2d(1,i)];
    }
    fout << "\n"; // separador entre timesteps

    fout.close();
}