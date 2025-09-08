#ifndef CONSTANTS_H
#define CONSTANTS_H

constexpr int CubeD = 40;
constexpr int Nx = 2000;
constexpr int Ny = 320;
constexpr int L1 = 12.5*CubeD;
constexpr int K = 9;

inline int index2d(int i, int j) {
    return i * Ny + j;
}

inline int index3D(int k, int i, int j) {
    return k * Nx * Ny + i * Ny + j;
}

inline int oppositeDirection(int k) {
    static const int opp[9] = {0, 3, 4, 1, 2, 7, 8, 5, 6};
    return opp[k];
}

#endif // CONSTANTS_H

