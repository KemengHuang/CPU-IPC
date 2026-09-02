#include "RuntimePaths.h"
#include "Simulator.h"

#include <Eigen/Core>

#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <stdexcept>
#include <string>

namespace {

struct CommandLineOptions {
    SimulationScene scene = SimulationScene::ClothOverBunny;
    int steps = 1;
    bool resume = false;
    bool writeOutput = true;
    bool writeCheckpoints = false;
    bool verbose = false;
    bool diagnoseLineSearch = false;
    bool disableBarrier = false;
    BroadPhaseBackend broadPhase = BroadPhaseBackend::LinearBVH;
    LinearSolverBackend linearSolver = LinearSolverOptions{}.backend;
    int cholmodThreadCount = LinearSolverOptions{}.cholmodThreadCount;
    int pardisoThreadCount = LinearSolverOptions{}.pardisoThreadCount;
    std::string outputDirectory;
};

void printUsage(const char* executable)
{
    std::cout
        << "Usage: " << executable << " [options]\n"
        << "  --scene cloth-bunny|twisting-mat|twisting-mat-soft|bunny2\n"
        << "  --steps N\n"
        << "  --output DIRECTORY\n"
        << "  --broad-phase spatial-hash|lbvh\n"
        << "  --linear-solver cholmod|pardiso|eigen-cg\n"
        << "    default: pardiso when available, otherwise optimized cholmod\n"
        << "  --cholmod-threads N (default 0: auto-select 4 or 8)\n"
        << "  --pardiso-threads N (default 16; 0 uses the oneMKL default)\n"
        << "  --disable-barrier\n"
        << "  --diagnose-line-search\n"
        << "  --resume\n"
        << "  --write-checkpoints\n"
        << "  --no-output\n"
        << "  --verbose\n"
        << "  --help\n";
}

int parseNonNegativeInt(const std::string& value, const char* option)
{
    size_t parsed = 0;
    const int result = std::stoi(value, &parsed);
    if (parsed != value.size() || result < 0) {
        throw std::invalid_argument(std::string(option) + " expects a non-negative integer");
    }
    return result;
}

CommandLineOptions parseCommandLine(int argc, char** argv)
{
    CommandLineOptions options;
    for (int i = 1; i < argc; ++i) {
        const std::string argument = argv[i];
        auto requireValue = [&](const char* option) -> std::string {
            if (++i >= argc) {
                throw std::invalid_argument(std::string(option) + " requires a value");
            }
            return argv[i];
        };

        if (argument == "--scene") {
            options.scene = parseSimulationScene(requireValue("--scene"));
        }
        else if (argument == "--steps") {
            options.steps = parseNonNegativeInt(requireValue("--steps"), "--steps");
        }
        else if (argument == "--output") {
            options.outputDirectory = requireValue("--output");
        }
        else if (argument == "--broad-phase") {
            const std::string value = requireValue("--broad-phase");
            if (value == "spatial-hash") {
                options.broadPhase = BroadPhaseBackend::SpatialHash;
            }
            else if (value == "lbvh") {
                options.broadPhase = BroadPhaseBackend::LinearBVH;
            }
            else {
                throw std::invalid_argument("unknown broad phase: " + value);
            }
        }
        else if (argument == "--linear-solver") {
            const std::string value = requireValue("--linear-solver");
            if (value == "cholmod") {
                options.linearSolver = LinearSolverBackend::Cholmod;
            }
            else if (value == "pardiso") {
                options.linearSolver = LinearSolverBackend::Pardiso;
            }
            else if (value == "eigen-cg") {
                options.linearSolver = LinearSolverBackend::EigenConjugateGradient;
            }
            else {
                throw std::invalid_argument("unknown linear solver: " + value);
            }
        }
        else if (argument == "--pardiso-threads") {
            options.pardisoThreadCount = parseNonNegativeInt(
                requireValue("--pardiso-threads"), "--pardiso-threads");
        }
        else if (argument == "--cholmod-threads") {
            options.cholmodThreadCount = parseNonNegativeInt(
                requireValue("--cholmod-threads"), "--cholmod-threads");
        }
        else if (argument == "--disable-barrier") {
            options.disableBarrier = true;
        }
        else if (argument == "--diagnose-line-search") {
            options.diagnoseLineSearch = true;
        }
        else if (argument == "--resume") {
            options.resume = true;
        }
        else if (argument == "--write-checkpoints") {
            options.writeCheckpoints = true;
        }
        else if (argument == "--no-output") {
            options.writeOutput = false;
        }
        else if (argument == "--verbose") {
            options.verbose = true;
        }
        else if (argument == "--help" || argument == "-h") {
            printUsage(argv[0]);
            std::exit(0);
        }
        else {
            throw std::invalid_argument("unknown option: " + argument);
        }
    }
    return options;
}

void printFinalSummary(const FEMSimulator& simulator, int completedSteps)
{
    const SimulationModel& model = simulator.getModel();
    if (model.meshes.empty()) {
        return;
    }

    const mesh3D& mesh = model.meshes.front();
    Eigen::Vector3d coordinateSum = Eigen::Vector3d::Zero();
    double squaredNormSum = 0.0;
    for (const Eigen::Vector3d& position : mesh.vertexes) {
        coordinateSum += position;
        squaredNormSum += position.squaredNorm();
    }

    const IPCStepStats& stats = simulator.getLastStepStats();
    std::cout << std::setprecision(17)
        << "RESULT steps=" << completedSteps
        << " vertices=" << mesh.vertexes.size()
        << " sum_x=" << coordinateSum.x()
        << " sum_y=" << coordinateSum.y()
        << " sum_z=" << coordinateSum.z()
        << " squared_norm_sum=" << squaredNormSum;
    if (completedSteps > 0) {
        std::cout << " last_step_ms=" << stats.stepMilliseconds
            << " last_newton=" << stats.newtonIterations
            << " last_min_distance2=" << stats.minConstraintDistance2
            << " last_matrix_nnz=" << stats.matrixNonZeros;
    }
    std::cout << std::endl;
}

} // namespace

int main(int argc, char** argv)
{
    try {
        const CommandLineOptions commandLine = parseCommandLine(argc, argv);
        if (!commandLine.outputDirectory.empty()
            && !RuntimePaths::setOutputDirectory(commandLine.outputDirectory)) {
            throw std::runtime_error("invalid output directory");
        }
        else if (commandLine.outputDirectory.empty()) {
            const std::string defaultHeadlessOutput = RuntimePaths::outputFile("headless");
            RuntimePaths::setOutputDirectory(defaultHeadlessOutput);
        }

        SimulationOptions simulationOptions;
        simulationOptions.scene = commandLine.scene;
        simulationOptions.resumeCheckpoint = commandLine.resume;
        simulationOptions.writeRuntimeFiles = commandLine.writeOutput;
        simulationOptions.writeCheckpoints = commandLine.writeCheckpoints;
        simulationOptions.verbose = commandLine.verbose;
        simulationOptions.broadPhaseBackend = commandLine.broadPhase;
        simulationOptions.linearSolver.backend = commandLine.linearSolver;
        simulationOptions.linearSolver.cholmodThreadCount = commandLine.cholmodThreadCount;
        simulationOptions.linearSolver.pardisoThreadCount = commandLine.pardisoThreadCount;
        simulationOptions.diagnoseLineSearch = commandLine.diagnoseLineSearch;
        simulationOptions.disableBarrier = commandLine.disableBarrier;

        FEMSimulator simulator;
        if (!simulator.buildModels(simulationOptions)) {
            throw std::runtime_error("failed to build simulation scene");
        }

        int stepId = 0;
        for (int step = 0; step < commandLine.steps; ++step) {
            simulator.simulateStep(stepId);
        }
        printFinalSummary(simulator, commandLine.steps);
    }
    catch (const std::exception& error) {
        std::cerr << "cipc_headless: " << error.what() << std::endl;
        return 1;
    }
    return 0;
}
