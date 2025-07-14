#ifndef INITIALIZE_H
#define INITIALIZE_H

#include "LBMParams.h"
#include "DomainParams.h"

void Initialize(LBMParams &params, const DomainParams &geometry, bool Parabolic, bool Cube);

#endif // INITIALIZE_H

