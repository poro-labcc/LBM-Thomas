#ifndef CONSTANTS_H
#define CONSTANTS_H

constexpr int RE_MAIN = 300;
constexpr int RUN_TIME = 350000;
constexpr int PROBE_TIME = 118;

//constexpr float t = 1/(8/2.5);
constexpr int CubeD   = 40;
constexpr int Nx      = static_cast<int>(50 * CubeD);

// Altura física do canal (em unidades de malha LBM)
constexpr int H = static_cast<int>(8 * CubeD);

// Número de nós fluidos
constexpr int Ny_fluid = H;

// Número total de nós (fluido + sólidos para bounce-back half-way)
constexpr int Ny = Ny_fluid + 2;

// Posição física das paredes (em unidades de malha)
constexpr double y_wall_bottom = 0.0;
constexpr double y_wall_top    = H;

constexpr int L1 = static_cast<int>(12.5*CubeD);
constexpr int K  = 9;

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
