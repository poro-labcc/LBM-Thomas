#include "./Boundary.h"
#include <cmath>
#include <omp.h>
#include "../SimStructure/Constants.h"
#include "./BoundaryConditions.h"

void Boundary(LBMParams &params, BoundaryConditionType type) {
	//applyZouHeCondition(params);
	applyBreuerCondition(params); //Inlet
	applyBoundaryCondition(params, type);
	//applyBreuerConditionOUT(params);
}
