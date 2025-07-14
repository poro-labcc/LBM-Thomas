#ifndef RUNSIMULATION_H
#define RUNSIMULATION_H

#include "./SimStructure/LBMParams.h"
#include "./PostProcess/SimulationStats.h"
#include "SimStructure/DomainParams.h"
#include "./Boundary/BoundaryConditions.h"
#include <optional>

void RunSimulation(LBMParams &params, SimulationStats &stats, const DomainParams &domainParams, BoundaryConditionType bc_type, std::optional<int> ReOnly = std::nullopt);

#endif // RUNSIMULATION_H

