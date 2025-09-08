#ifndef BOUNDARYCONDITIONS_H
#define BOUNDARYCONDITIONS_H

#include "../SimStructure/LBMParams.h"
#include <string>

enum class BoundaryConditionType {
    EMERICH,
    SECOND_ORDER_EXTRAPOLATION,
    CONVECTIVE,
    EXIT_INLET
};

// Interface for boundary conditions
void applyBoundaryCondition(LBMParams &params, BoundaryConditionType type);

// Individual boundary condition implementations
void applyEmerichBoundary(LBMParams &params);
void applySecondOrderExtrapolationBoundary(LBMParams &params);
void applyConvectiveBoundaryNardelli(LBMParams &params);
void applyConvectiveBoundaryWang(LBMParams &params);
void applyExitInletBoundary(LBMParams &params);
void applyReisBoundary(LBMParams &params);
void applyZouHeCondition(LBMParams &params);
void applyBreuerCondition(LBMParams &params);
std::string boundaryConditionToString(BoundaryConditionType type);

#endif // BOUNDARYCONDITIONS_H

