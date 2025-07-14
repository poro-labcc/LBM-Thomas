#ifndef SAVEFUNCTIONS_H
#define SAVEFUNCTIONS_H

#include "../SimStructure/LBMParams.h"
#include "./SimulationStats.h"

void saveField(LBMParams &params, const int time, int step);
void saveForce(SimulationStats &stats, const int time, int step, int Re);
void SaveVTK(int timestep, const LBMParams &params);

#endif // SAVEFUNCTIONS_H

