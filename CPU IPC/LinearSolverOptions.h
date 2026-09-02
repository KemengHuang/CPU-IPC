#pragma once

enum class LinearSolverBackend {
    Cholmod,
    SuiteSparseLDL,
    Pardiso,
    EigenConjugateGradient
};

struct LinearSolverOptions {
#ifdef CIPC_HAS_PARDISO
    LinearSolverBackend backend = LinearSolverBackend::Pardiso;
#else
    LinearSolverBackend backend = LinearSolverBackend::SuiteSparseLDL;
#endif
    double relativeTolerance = 1e-6;
    int maximumIterations = 10000;
    int cholmodThreadCount = 0;
    int pardisoThreadCount = 16;
};
