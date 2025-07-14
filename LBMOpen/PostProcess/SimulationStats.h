#ifndef SIMULATIONSTATS_H
#define SIMULATIONSTATS_H

#include "../SimStructure/LBMParams.h"
#include "../SimStructure/DomainParams.h"

struct SimulationStats {
    double Fx, Fy; // Forces in x and y directions
    double Cd, Cl; // Drag and lift coefficients
    double inletMassFlow; // Mass flow at the inlet
    double outletMassFlow; // Mass flow at the outlet
    double relativeDifference; // Relative difference in the distribution functions

    // Constructor to initialize members
    SimulationStats();

    // Method to calculate mass flow at inlet and outlet
    void calculateMassFlow(const LBMParams &params);

    // Method to compute the relative difference in distribution functions
    void computeRelativeDifference(const LBMParams &params);

    // Method to calculate drag and lift coefficients
    void calculateCoefficients(const LBMParams &params, const DomainParams &dim);

    // Method to reset statistics
    void reset();
};

#endif // SIMULATIONSTATS_H

