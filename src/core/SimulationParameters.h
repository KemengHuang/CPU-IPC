#pragma once
#ifndef SIMULATION_PARAMETERS_H
#define SIMULATION_PARAMETERS_H

// Centralized default simulation parameters.
// All values are intentionally identical to the hardcoded literals previously
// scattered across Simulator.cpp and IPC_FUNC.cpp so that numerical behavior
// is preserved.
struct SimulationParameters {
    // Time stepping
    static constexpr double defaultDt = 1.0 / 240.0;          // 240 fps

    // Material defaults
    static constexpr double defaultDensity = 1e3;
    static constexpr double defaultClothDensity = 1e2;
    static constexpr double defaultYoungModulus = 1e4;
    static constexpr double defaultPoissonRate = 0.49;
    static constexpr double defaultClothYoungModulus = 1e4;
    static constexpr double defaultBendYoungModulus = 1e6;
    static constexpr double defaultClothThickness = 1e-3;
    static constexpr double defaultStrainRate = 1e2;

    // Damping / friction
    static constexpr double defaultDragCoefficient = 0.00;
    static constexpr double defaultFrictionCoefficient = 0.5;

    // Barrier / contact defaults
    static constexpr double defaultHhat = 9e-8;
    static constexpr double defaultFhat = 1e-6;
    static constexpr double defaultDTol = 1e-18;
    static constexpr double defaultKappa = 0.0;
    static constexpr double defaultNewtonSolverThreshold = 1e-2;

    // External forces
    static constexpr double gravityMagnitude = -9.8;

    // IPC / CCD parameters
    static constexpr double CCDDistRatio = 0.2;               // tet injective step slackness
    static constexpr double environmentSlackness = 0.8;
    static constexpr double selfContactSlackness = 0.8;

    // Newton solver convergence: tolerance factor inside sqrt(...)
    static constexpr double newtonConvergenceTolFactor = 1e-4;

    // Line search: lower-bound fraction of the feasible step size
    static constexpr double lineSearchMinFraction = 1e-3;

    // Kappa initialization / bounding
    static constexpr double kappaMaxMassFactor = 1e13;
    static constexpr double kappaUpperBoundFactor = 100.0;

    // Newton iteration cap
    static constexpr int newtonIterCap = 10000;
};

#endif // SIMULATION_PARAMETERS_H
