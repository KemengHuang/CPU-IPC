#pragma once
#ifndef CIPC_IPC_SOLVER_H
#define CIPC_IPC_SOLVER_H

#include "Elasticity.h"
#include "LinearSolverOptions.h"
#include <cstddef>
#include <memory>

class NewtonLinearSystem;

struct IPCStepStats {
    int frame = 0;
    int newtonIterations = 0;
    int kappaIterations = 0;
    int lineSearchBacktracks = 0;
    int energyBacktracks = 0;
    int intersectionBacktracks = 0;
    int maximumEnergyBacktracksPerNewton = 0;
    int newtonStepsWithMoreThanTwoBacktracks = 0;
    double stepMilliseconds = 0.0;
    double assemblyMilliseconds = 0.0;
    double linearSolveMilliseconds = 0.0;
    double ccdMilliseconds = 0.0;
    double lineSearchMilliseconds = 0.0;
    double postLineSearchMilliseconds = 0.0;
    double pardisoAnalysisMilliseconds = 0.0;
    double pardisoFactorizationMilliseconds = 0.0;
    double pardisoSolveMilliseconds = 0.0;
    double collisions = 0.0;
    double kappa = 0.0;
    double minConstraintDistance2 = -1.0;
    double acceptedStepSum = 0.0;
    double minimumAcceptedStep = 1.0;
    double maximumAcceptedStep = 0.0;
    std::size_t groundContacts = 0;
    std::size_t selfContacts = 0;
    std::size_t mollifiedContacts = 0;
    std::size_t matrixNonZeros = 0;
    std::size_t symbolicAnalyses = 0;
    std::size_t numericFactorizations = 0;
    std::size_t factorNonZeros = 0;
    int linearSolverThreads = 0;
};

struct IPCSolverContext {
    bool checkpointLoadAttempted = false;
    bool writeRuntimeFiles = true;
    bool writeCheckpoints = true;
    bool verbose = true;
    bool diagnoseLineSearch = false;
    LinearSolverOptions linearSolver;
    std::shared_ptr<NewtonLinearSystem> linearSystem;
    bool metricsInitialized = false;
    int stepIndex = 0;
    int totalNewtonIterations = 0;
    double cumulativeStepMilliseconds = 0.0;
    double cumulativeStageMilliseconds[5] = { 0.0, 0.0, 0.0, 0.0, 0.0 };
    double totalCollisions = 0.0;
    IPCStepStats lastStep;
};

int solveIPCStep(
    int& stepId,
    mesh3D& mesh,
    SpatialHash& broadPhase,
    Ground& ground,
    IPCSolverContext& context);
void updateInertialTarget(mesh3D& mesh);
#endif // CIPC_IPC_SOLVER_H
