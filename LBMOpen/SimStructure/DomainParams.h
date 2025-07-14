#ifndef DOMAINPARAMS_H
#define DOMAINPARAMS_H

#include "./Constants.h"

struct DomainParams {
    int CubeD, L1, middle;

    DomainParams(int cubeD, int l1, int ny);
};

#endif // DOMAINPARAMS_H

