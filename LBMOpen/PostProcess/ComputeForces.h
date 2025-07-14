#ifndef COMPUTEFORCES_H
#define COMPUTEFORCES_H

#include "../SimStructure/LBMParams.h"
#include "SimulationStats.h"
#include "../SimStructure/DomainParams.h"

void ComputeForces(const LBMParams &params, SimulationStats &stats, const DomainParams &geometry);

#endif // COMPUTEFORCES_H

