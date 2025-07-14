#include <iostream>
#include <string>
#include <sstream>
#include <cmath>
#include <omp.h>
#include <fstream>
#include <vector>
#include <thread>
#include <chrono>
#include <iomanip>
#include <algorithm>
#include <cassert>
#include <optional>

#include "./SimStructure/Constants.h"
#include "./SimStructure/LBMParams.h"
#include "./SimStructure/DomainParams.h"
#include "./PostProcess/SimulationStats.h"
#include "RunSimulation.h"
#include "./Boundary/BoundaryConditions.h"

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
    RunSimulation(params, stats, domainParams, BoundaryConditionType::EXIT_INLET, 300);

    // Calcular tempo de execução
    auto const end = std::chrono::high_resolution_clock::now();
    auto const duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);
    
    std::cout << "Simulation completed in " << duration.count() / 1000.0 << " seconds" << std::endl;
    
    return 0;
}

