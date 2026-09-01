#pragma once

#include "CollisionBroadPhase.h"
#include "IPCSolver.h"

enum class SimulationScene {
    TwistingMat,
    ClothOverBunny,
    Bunny2
};

struct SimulationOptions {
    SimulationScene scene = SimulationScene::ClothOverBunny;
    bool resumeCheckpoint = false;
    bool writeRuntimeFiles = true;
    bool writeCheckpoints = false;
    bool verbose = true;
    bool diagnoseLineSearch = false;
    bool disableBarrier = false;
    BroadPhaseBackend broadPhaseBackend = BroadPhaseBackend::LinearBVH;
    LinearSolverOptions linearSolver;
};

class FEMSimulator {
public:
    bool buildModels(SimulationScene scene = SimulationScene::ClothOverBunny);
    bool buildModels(const SimulationOptions& options);

    SimulationModel& getModel() { return model_; }
    const SimulationModel& getModel() const { return model_; }
    const IPCStepStats& getLastStepStats() const { return solverContext_.lastStep; }
    int simulateStep(int& stepId);

private:
    void rebuildCollisionSets();

    SimulationModel model_;
    SpatialHash broadPhase_;
    Ground ground_;
    IPCSolverContext solverContext_;
};
