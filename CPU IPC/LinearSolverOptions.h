#pragma once

enum class LinearSolverBackend {
    Cholmod,
    EigenConjugateGradient
};

struct LinearSolverOptions {
    LinearSolverBackend backend = LinearSolverBackend::Cholmod;
    double relativeTolerance = 1e-6;
    int maximumIterations = 10000;
};
