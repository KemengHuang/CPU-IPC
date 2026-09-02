#include "Simulator.h"
#include "BoundaryConditions.h"
#include "iostream"
#include "Eigen/Eigen"
#include "RuntimePaths.h"
#include <fstream>
#include <cmath>
#include <stdexcept>
#include <unordered_map>
#include <unordered_set>

using namespace Eigen;

namespace {

constexpr double pi = 3.14159265358979323846;
constexpr double twistingSoftConstraintWeight = 100.0;

enum class TwistingBoundaryMode {
    HardDirichlet,
    SoftTarget
};

Vector3d twistingBoundaryDirection(
    const Vector3d& vertex,
    const Vector3d& restVertex,
    double stepFraction,
    double timeStep)
{
    const double angleX = 3.14 / 5 * timeStep * stepFraction;
    Matrix3d rotationL, rotationR;
    rotationL << 1, 0, 0,
        0, cos(angleX), sin(angleX),
        0, -sin(angleX), cos(angleX);
    rotationR << 1, 0, 0,
        0, cos(angleX), -sin(angleX),
        0, sin(angleX), cos(angleX);

    if (restVertex[0] < 0.0) {
        return rotationL * vertex - vertex;
    }
    if (restVertex[0] > 0.0) {
        return rotationR * vertex - vertex;
    }
    return Vector3d::Zero();
}

void prepareQuadraticBending(mesh3D& mesh)
{
    mesh.quadBendingInfo.clear();
    std::unordered_map<int64_t, int> edgeMap;

    auto edgeHash = [](const Vector2i& edge) -> int64_t {
        return (int64_t(edge[0]) << 32) | edge[1];
        };

    for (const auto& triangle : mesh.triangles)
    {
        for (int i = 0; i != 3; ++i)
        {
            Vector2i edge = Vector2i(triangle[i], triangle[(i + 1) % 3]);
            int otherVert = triangle[(i + 2) % 3];
            if (edge[0] > edge[1])
                std::swap(edge[0], edge[1]);

            if (edgeMap.find(edgeHash(edge)) == edgeMap.end())
                edgeMap[edgeHash(edge)] = otherVert;
            else  // found a quad
            {
                int trig1OuterIndex = edgeMap[edgeHash(edge)];
                int trig2OuterIndex = otherVert;

                using Vec3 = Vector3d;
                Vec3 v0 = mesh.vertexes[trig1OuterIndex],
                    v1 = mesh.vertexes[edge[0]],
                    v2 = mesh.vertexes[edge[1]],
                    v3 = mesh.vertexes[trig2OuterIndex];

                Vec3 e0 = v2 - v1,
                    e1 = v0 - v1,
                    e2 = v3 - v1,
                    e3 = v0 - v2,
                    e4 = v3 - v2;

                double a0 = 0.5 * e0.cross(e1).norm(),
                    a1 = 0.5 * e0.cross(e2).norm();

                auto cot = [](const Vec3 e0, const Vec3 e1) -> double {
                    return e0.dot(e1) / e0.cross(e1).norm();
                    };

                double c01 = cot(e0, e1),
                    c02 = cot(e0, e2),
                    c03 = cot(-e0, e3),
                    c04 = cot(-e0, e4);

                Vector4d K(
                    c03 + c04,
                    c01 + c02,
                    -c01 - c03,
                    -c02 - c04
                );

                Matrix4d Q = 3.0 / (a0 + a1) * K * K.transpose();
                mesh.quadBendingInfo.push_back(QuadBendingInfo(edge[0], edge[1], trig1OuterIndex, trig2OuterIndex, Q));
            }
        }
    }
}

void prepareHingeBending(mesh3D& mesh)
{
    mesh.hingeBendingInfo.clear();
    mesh.hingeBendingInfo.reserve(mesh.tri_edges.size());
    for (int edgeIndex = 0;
         edgeIndex < static_cast<int>(mesh.tri_edges.size());
         ++edgeIndex) {
        const Eigen::Vector2i& edge = mesh.tri_edges[edgeIndex];
        const Eigen::Vector2i& opposite = mesh.tri_edges_adj_points[edgeIndex];
        if (opposite[0] < 0 || opposite[1] < 0) {
            continue;
        }
        mesh.hingeBendingInfo.emplace_back(makeHingeBendingInfo(
            edge[0], edge[1], opposite[0], opposite[1], mesh.vertexes));
    }
}

void updateMaterial(mesh3D& mesh)
{
    mesh.lengthRateLame = mesh.YoungModulus / (2 * (1 + mesh.PoissonRate));
    mesh.volumeRateLame = mesh.YoungModulus * mesh.PoissonRate
        / ((1 + mesh.PoissonRate) * (1 - 2 * mesh.PoissonRate));
    mesh.lengthRate = 4 * mesh.lengthRateLame / 3;
    mesh.volumeRate = mesh.volumeRateLame + 5 * mesh.lengthRateLame / 6;


    mesh.stretchStiffness = mesh.clothYoungModulus / (1 - mesh.PoissonRate * mesh.PoissonRate);
    mesh.shearStiffness = mesh.shearYoungModulus / (2 * (1 + mesh.PoissonRate)) / mesh.clothThicness;
    mesh.plateRigidity = mesh.bendYoungModulus * std::pow(mesh.clothThicness, 3)
        / (12 * (1 - mesh.PoissonRate * mesh.PoissonRate));
}

void applyDefaultSettings(mesh3D& mesh3d) {

    mesh3d.density = 1e3;
    mesh3d.cloth_density = 1e2;
    mesh3d.clothThicness = 1e-3;

    mesh3d.use_barrier = 1;
    mesh3d.IPC_dt = 1 / 240.; // 240 fps

    mesh3d.drag_coeff = 0.00;
    mesh3d.Hhat = 9e-8;

    mesh3d.Fhat = 1e-6;
    mesh3d.Kappa = 0;
    mesh3d.dTol = 1e-18;

    mesh3d.YoungModulus = 1e4;
    mesh3d.PoissonRate = 0.49;
    mesh3d.friction = 0.5;

    mesh3d.clothYoungModulus = 1e4;
    mesh3d.shearYoungModulus = 1e1;
    mesh3d.bendYoungModulus = 1e4;
    mesh3d.strainRate = 1e2;

    mesh3d.Newton_Solver_Threshold = 1e-2;

    updateMaterial(mesh3d);
}


void loadSettings(mesh3D& mesh3d) {
    const std::string configFile = std::string{ CIPC_ASSETS_DIR } + "scene/parameterSetting.txt";
    std::ifstream infile(configFile);
    if (!infile) {
        std::cerr << "Warning: failed loading settings; using defaults." << std::endl;
        applyDefaultSettings(mesh3d);
        return;
    }

    const std::unordered_set<std::string> requiredKeys = {
        "volume_mesh_density",
        "volume_mesh_YoungsModulus",
        "poisson_rate",
        "friction_rate",
        "triangle(cloth)_mesh_thickness",
        "triangle(cloth)_mesh_YoungsModulus",
        "triangle(cloth)_mesh_shearYoungsModulus",
        "triangle(cloth)_mesh_bendingYoungsModulus",
        "triangle(cloth)_mesh_density",
        "strainRate",
        "enable_collision_handling",
        "drag_coeff",
        "ipc_time_step",
        "Newton_solver_threshold",
        "IPC_ralative_dHat"
    };
    std::unordered_map<std::string, double> values;
    std::string key;
    double value = 0.0;
    while (infile >> key) {
        if (!(infile >> value)) {
            throw std::runtime_error("invalid value for setting '" + key + "'");
        }
        if (requiredKeys.count(key) == 0) {
            throw std::runtime_error("unknown setting '" + key + "'");
        }
        if (!values.emplace(key, value).second) {
            throw std::runtime_error("duplicate setting '" + key + "'");
        }
    }
    for (const std::string& required : requiredKeys) {
        if (values.count(required) == 0) {
            throw std::runtime_error("missing setting '" + required + "'");
        }
    }

    auto requirePositive = [&](const char* setting) {
        const double settingValue = values.at(setting);
        if (!std::isfinite(settingValue) || settingValue <= 0.0) {
            throw std::runtime_error(std::string(setting) + " must be finite and positive");
        }
        return settingValue;
    };
    auto requireNonNegative = [&](const char* setting) {
        const double settingValue = values.at(setting);
        if (!std::isfinite(settingValue) || settingValue < 0.0) {
            throw std::runtime_error(std::string(setting) + " must be finite and non-negative");
        }
        return settingValue;
    };

    mesh3d.density = requirePositive("volume_mesh_density");
    mesh3d.YoungModulus = requirePositive("volume_mesh_YoungsModulus");
    mesh3d.PoissonRate = values.at("poisson_rate");
    if (!std::isfinite(mesh3d.PoissonRate)
        || mesh3d.PoissonRate <= -1.0 || mesh3d.PoissonRate >= 0.5) {
        throw std::runtime_error("poisson_rate must be in (-1, 0.5)");
    }
    mesh3d.friction = requireNonNegative("friction_rate");
    mesh3d.clothThicness = requirePositive("triangle(cloth)_mesh_thickness");
    mesh3d.clothYoungModulus = requirePositive("triangle(cloth)_mesh_YoungsModulus");
    mesh3d.shearYoungModulus = requirePositive("triangle(cloth)_mesh_shearYoungsModulus");
    mesh3d.bendYoungModulus = requireNonNegative("triangle(cloth)_mesh_bendingYoungsModulus");
    mesh3d.cloth_density = requirePositive("triangle(cloth)_mesh_density");
    mesh3d.strainRate = requireNonNegative("strainRate");
    const double collisionFlag = values.at("enable_collision_handling");
    if (collisionFlag != 0.0 && collisionFlag != 1.0) {
        throw std::runtime_error("enable_collision_handling must be 0 or 1");
    }
    mesh3d.use_barrier = collisionFlag != 0.0;
    mesh3d.drag_coeff = requireNonNegative("drag_coeff");
    mesh3d.IPC_dt = requirePositive("ipc_time_step");
    mesh3d.Newton_Solver_Threshold = requirePositive("Newton_solver_threshold");
    const double relativeDHat = requirePositive("IPC_ralative_dHat");
    mesh3d.Hhat = relativeDHat * relativeDHat;
    mesh3d.Fhat = 1e-6;
    mesh3d.Kappa = 0.0;
    mesh3d.dTol = 1e-18;
    updateMaterial(mesh3d);
}

void buildTwistingMatScene(
    mesh3D& mesh3d,
    const std::string& assetDirectory,
    TwistingBoundaryMode boundaryMode) {

    mesh3d.YoungModulus = 1e6;
    updateMaterial(mesh3d);
    mesh3d.objMaxConer = Vector3d(0, 0, 0);
    mesh3d.objMinConer = Vector3d(0, 0, 0);
	mesh3d.is_quasi_static = true;
	mesh3d.apply_gravity = false;
    if (!mesh3d.load_tetrahedraMesh(
            assetDirectory + "tetrahedraMesh/ipcmesh/mat40x40.msh",
            1,
            Vector3d(0, 0, 0))) {
        throw std::runtime_error("failed to load twisting-mat mesh");
    }

    const double eps = 1e-4;
    for (int i = 0; i < mesh3d.vertexes.size(); i++)
    {
        if (mesh3d.vertexes[i][0] < mesh3d.objMinConer[0] + eps || mesh3d.vertexes[i][0] > mesh3d.objMaxConer[0] - eps)
        {
            if (boundaryMode == TwistingBoundaryMode::HardDirichlet) {
                mesh3d.boundaryConditions.dirichlet.vertexIndices.push_back(i);
                mesh3d.boundaryTypes[i] =
                    boundaryTypeCode(VertexBoundaryType::Dirichlet);
            }
            else {
                mesh3d.boundaryConditions.soft.vertexIndices.push_back(i);
            }
        }
    }

    if (boundaryMode == TwistingBoundaryMode::HardDirichlet) {
        mesh3d.boundaryConditions.dirichlet.updateDirection =
            [](const Vector3d& vertex,
                const Vector3d& restVertex,
                int,
                double alpha,
                double ipc_dt) -> Vector3d
            {
                return twistingBoundaryDirection(
                    vertex, restVertex, alpha, ipc_dt);
            };
    }
    else {
        mesh3d.boundaryConditions.soft.weight = twistingSoftConstraintWeight;
        mesh3d.boundaryConditions.soft.updateTarget =
            [](const Vector3d& vertex,
                const Vector3d& restVertex,
                int,
                double ipc_dt) -> Vector3d
            {
                // Use the same endpoint as the hard x_new = x - p update,
                // but attract the free vertex to it with a finite penalty.
                return vertex - twistingBoundaryDirection(
                    vertex, restVertex, 1.0, ipc_dt);
            };
    }
}

void buildClothOverBunnyScene(mesh3D& mesh3d, const std::string& assetDirectory) {

    mesh3d.YoungModulus = 1e4;
    updateMaterial(mesh3d);
    mesh3d.objMaxConer = Vector3d(0, 0, 0);
    mesh3d.objMinConer = Vector3d(0, 0, 0);
    //mesh3d.is_quasi_static = true;
    {
        if (!mesh3d.load_triangleMesh(
                assetDirectory + "triangleMesh/planes/plane1024.obj",
                1,
                Vector3d(0, 0, 0))) {
            throw std::runtime_error("failed to load cloth mesh");
        }

#if defined USE_QUADRATIC_BENDING
        prepareQuadraticBending(mesh3d);
#else
        prepareHingeBending(mesh3d);
#endif

        Matrix3d rotate;
        const double angleX = pi / 2;
        rotate << 1, 0, 0, 0, cos(angleX), -sin(angleX), 0, sin(angleX), cos(angleX);
        for (int j = 0; j < mesh3d.vertexNum; j++) {
            mesh3d.vertexes[j] = (rotate * mesh3d.vertexes[j]);

        }
        double xmin = 1e32, ymin = 1e32, zmin = 1e32;
        double xmax = -1e32, ymax = -1e32, zmax = -1e32;
        for (int j = 0; j < mesh3d.vertexNum; j++) {

            Vector3d pos = mesh3d.vertexes[j];
            if (xmin > pos[0]) xmin = pos[0];
            if (ymin > pos[1]) ymin = pos[1];
            if (zmin > pos[2]) zmin = pos[2];
            if (xmax < pos[0]) xmax = pos[0];
            if (ymax < pos[1]) ymax = pos[1];
            if (zmax < pos[2]) zmax = pos[2];


        }
        mesh3d.maxConer = Vector3d(xmax, ymax, zmax);
        mesh3d.minConer = Vector3d(xmin, ymin, zmin);
    }

    if (!mesh3d.load_tetrahedraMesh(
            assetDirectory + "tetrahedraMesh/bunny.msh",
            0.5,
            Vector3d(0, -0.5, 0))) {
        throw std::runtime_error("failed to load bunny mesh");
    }

}

void buildBunny2Scene(mesh3D& mesh3d, const std::string& assetDirectory)
{
    mesh3d.YoungModulus = 1e5;
    updateMaterial(mesh3d);
    mesh3d.objMaxConer = Vector3d::Zero();
    mesh3d.objMinConer = Vector3d::Zero();

    constexpr double scale = 0.2;
    const std::string meshPath = assetDirectory + "tetrahedraMesh/bunny2.msh";
    if (!mesh3d.load_tetrahedraMesh(
            meshPath, scale, Vector3d(0.0, 0.65, 0.0))) {
        throw std::runtime_error("failed to load upper bunny2 mesh");
    }
    if (!mesh3d.load_tetrahedraMesh(
            meshPath, scale, Vector3d::Zero())) {
        throw std::runtime_error("failed to load lower bunny2 mesh");
    }
}

} // namespace

SimulationScene parseSimulationScene(const std::string& name)
{
    if (name == "cloth-bunny") {
        return SimulationScene::ClothOverBunny;
    }
    if (name == "twisting-mat") {
        return SimulationScene::TwistingMat;
    }
    if (name == "twisting-mat-soft") {
        return SimulationScene::TwistingMatSoft;
    }
    if (name == "bunny2") {
        return SimulationScene::Bunny2;
    }
    throw std::invalid_argument("unknown scene: " + name);
}

bool FEMSimulator::buildModels(SimulationScene scene) {
    SimulationOptions options;
    options.scene = scene;
    return buildModels(options);
}

bool FEMSimulator::buildModels(const SimulationOptions& options) {
    if (!RuntimePaths::initialize()) {
        std::cerr << "Warning: failed to create runtime output directories under "
                  << RuntimePaths::outputFile("") << std::endl;
    }

    model_.meshes.clear();

    if (!std::isfinite(options.linearSolver.relativeTolerance)
        || options.linearSolver.relativeTolerance <= 0.0
        || options.linearSolver.maximumIterations <= 0
        || options.linearSolver.cholmodThreadCount < 0
        || options.linearSolver.pardisoThreadCount < 0) {
        throw std::invalid_argument("linear solver options are invalid");
    }

    mesh3D mesh3d;
    loadSettings(mesh3d);
    if (options.disableBarrier) {
        mesh3d.use_barrier = false;
    }
    const std::string assetDirectory{ CIPC_ASSETS_DIR };
    mesh3d.maxConer = Vector3d(-1e32, -1e32, -1e32);
    mesh3d.minConer = Vector3d(1e32, 1e32, 1e32);

    switch (options.scene) {
    case SimulationScene::TwistingMat:
        buildTwistingMatScene(
            mesh3d, assetDirectory, TwistingBoundaryMode::HardDirichlet);
        break;
    case SimulationScene::TwistingMatSoft:
        buildTwistingMatScene(
            mesh3d, assetDirectory, TwistingBoundaryMode::SoftTarget);
        break;
    case SimulationScene::ClothOverBunny:
        buildClothOverBunnyScene(mesh3d, assetDirectory);
        break;
    case SimulationScene::Bunny2:
        buildBunny2Scene(mesh3d, assetDirectory);
        break;
    }

    mesh3d.vertexNum = mesh3d.vertexes.size();
    mesh3d.tetrahedraNum = mesh3d.tetrahedras.size();
    mesh3d.triangleNum = mesh3d.triangles.size();

    if (options.verbose) {
        printf("tets num: %d\n", mesh3d.tetrahedraNum);
        printf("verts num: %d\n", mesh3d.vertexNum);
    }

    initMesh3D(mesh3d, 1, 0.2);

    mesh3d.v_rest = mesh3d.vertexes;
    mesh3d.V_prev = mesh3d.vertexes;
    BoundaryConditionOps::initialize(mesh3d);

    updateInertialTarget(mesh3d);


    mesh3d.bboxDiagSize2 = (mesh3d.maxConer - mesh3d.minConer).squaredNorm();
    mesh3d.Hhat *= mesh3d.bboxDiagSize2;
    mesh3d.Fhat *= mesh3d.bboxDiagSize2;;
    mesh3d.dTol *= mesh3d.bboxDiagSize2;


    model_.meshes.push_back(mesh3d);
    model_.calculateSurface();
    mesh3D& simulationMesh = model_.meshes.front();

    if (simulationMesh.surfEdges.empty()) {
        throw std::runtime_error("simulation surface has no edges");
    }

    double length = 0;
    for (const auto& edge : simulationMesh.surfEdges) {
        length += (simulationMesh.vertexes[edge.first]
            - simulationMesh.vertexes[edge.second]).norm();
    }
    simulationMesh.averageEdgeLenth =
        length / (simulationMesh.surfEdges.size() * 3);
    if (options.verbose) {
        printf("triangles num: %zu\n", simulationMesh.surface.size());
    }
    simulationMesh.resumedFromCheckpoint = options.resumeCheckpoint
        && simulationMesh.load_tetTempData();
    if (options.verbose) {
        printf(simulationMesh.resumedFromCheckpoint
            ? "load temp data\n"
            : "no temp data\n");
    }

    broadPhase_.setBackend(options.broadPhaseBackend);
    rebuildCollisionSets();

    solverContext_ = IPCSolverContext{};
    solverContext_.writeRuntimeFiles = options.writeRuntimeFiles;
    solverContext_.writeCheckpoints = options.writeCheckpoints;
    solverContext_.verbose = options.verbose;
    solverContext_.diagnoseLineSearch = options.diagnoseLineSearch;
    solverContext_.linearSolver = options.linearSolver;

    return true;
}


void FEMSimulator::rebuildCollisionSets() {
    mesh3D& mesh = model_.meshes.front();
    if (mesh.use_barrier) {
        broadPhase_.build(mesh, mesh.averageEdgeLenth);
        broadPhase_.calculateActivateSet(mesh);
    }

    ground_.calculateActivateSet(mesh);
}

int FEMSimulator::simulateStep(int &stepId) {
    return solveIPCStep(
        stepId,
        model_.meshes.front(),
        broadPhase_,
        ground_,
        solverContext_);
}
