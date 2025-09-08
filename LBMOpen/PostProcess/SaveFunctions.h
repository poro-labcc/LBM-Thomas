#ifndef SAVEFUNCTIONS_H
#define SAVEFUNCTIONS_H

#include "../SimStructure/LBMParams.h"
#include "./SimulationStats.h"
#include <fstream>
#include <string>
enum class OutputType {
    All,
    OnlyInOut,
    None
};

void saveField(LBMParams &params, const int time, int step);
void saveForce(SimulationStats &stats, const int time, int step, int Re);
void SaveVTK(int timestep, const LBMParams &params);
void SaveFlow(int timestep, const SimulationStats &stats);
void saveDistribution(int timestep, const LBMParams &params);
void salvarFAppend(LBMParams &params, const std::string &filename);
#endif // SAVEFUNCTIONS_H

