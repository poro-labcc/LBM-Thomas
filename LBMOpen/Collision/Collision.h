#ifndef COLLISION_H
#define COLLISION_H

#include "../SimStructure/LBMParams.h"

void CollisionBGK(LBMParams &params);
void CollisionMRT(LBMParams &params);
void CollisionTRT(LBMParams &params);
#endif // COLLISION_H

