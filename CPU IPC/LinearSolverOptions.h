#pragma once

enum class LinearSolverBackend {
    Cholmod,
    SuiteSparseLDL,
    Pardiso,
    EigenConjugateGradient
};

struct LinearSolverOptions {
    LinearSolverBackend backend = LinearSolverBackend::SuiteSparseLDL;
    double relativeTolerance = 1e-6;
    int maximumIterations = 10000;
    int pardisoThreadCount = 16;
};
