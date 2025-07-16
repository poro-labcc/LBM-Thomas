#include "./SaveFunctions.h"
#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>
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

void SaveFlow(int timestep, const LBMParams &params) {
    // Format the filename with timestep
    std::ostringstream filename;
    filename << "velocity/fieldFlow" << std::setw(5) << std::setfill('0') << timestep << ".vtk";

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
    file << "DIMENSIONS " << 2 << " " << params.Ny << " 1" << std::endl;
    file << "ORIGIN 0 0 0" << std::endl;
    file << "SPACING 1 1 1" << std::endl;
    file << "POINT_DATA " << 2 * params.Ny << std::endl;

    // Save velocity field as vectors
    file << "VECTORS Velocity float" << std::endl;

#pragma omp parallel
    {
        std::ostringstream localBuffer;

#pragma omp for collapse(2) ordered
        for (int j = 0; j < params.Ny; j++) {
            for (int i :{0,params.Nx-1}) {
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
    file.close();
}