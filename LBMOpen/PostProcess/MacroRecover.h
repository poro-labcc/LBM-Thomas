#ifndef MACRORECOVER_H
#define MACRORECOVER_H

#include "../SimStructure/LBMParams.h"

void MacroRecover(LBMParams &params);
void MacroRecoverPoint(LBMParams &params, int i, int j, double &u, double &v);
#endif // MACRORECOVER_H

