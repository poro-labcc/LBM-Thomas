#ifndef COLLISION_H
#define COLLISION_H

#include "../SimStructure/LBMParams.h"
#include <string>

enum class CollisionType {
    BGK,
    MRT,
    TRT
};
void Collision(LBMParams &params, CollisionType type);

void CollisionBGK(LBMParams &params);
void CollisionMRT(LBMParams &params);
void CollisionTRT(LBMParams &params);
std::string CollisionToString(CollisionType type);


#endif // COLLISION_H

