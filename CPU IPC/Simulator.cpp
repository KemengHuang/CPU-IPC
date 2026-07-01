#include "fem_parameters.h"
#include "Simulator.h"
#include "iostream"
#include "Eigen/Eigen"
#include <fstream>

using namespace Eigen;
using namespace FEM;
using namespace std;

void buildSpecialPoints(mesh3D &mesh) {
    Vector3d tp[6];
    tp[0] = Vector3d(-0.26554,0.67295,-0.23324);
    tp[1] = Vector3d(-0.05724, 0.69143, -0.10306);
    tp[2] = Vector3d(0.04884, 0.47191, 0.04952);
    tp[3] = Vector3d(-0.26145, -0.12829, -0.23713);
    tp[4] = Vector3d(-0.06139, -0.10531, -0.09838);
    tp[5] = Vector3d(0.015541, -0.32489, 0.02275);
    mesh.specialPointsArray.resize(6);
    for(int i=0;i<mesh.specialPointsArray.size();i++){
        double mindist = 1e32;
        for(int j=0;j<mesh.surfVerts.size();j++){
            double dist = (mesh.vertexes[mesh.surfVerts[j]]-tp[i]).norm();
            if(dist<mindist){
                mindist = dist;
                mesh.specialPointsArray[i] = mesh.surfVerts[j];
            }
        }
    }
}

void PrepQuadBending(mesh3D& mesh)
{
    std::unordered_map<int64_t, int> edgeMap;

    auto edgeHash = [](const Vector2i& edge) -> int64_t {
        return (int64_t(edge[0]) << 32) | edge[1];
        };

    for (const auto trig : mesh.triangles)
    {
        for (int i = 0; i != 3; ++i)
        {
            Vector2i edge = Vector2i(trig[i], trig[(i + 1) % 3]);
            int otherVert = trig[(i + 2) % 3];
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

void DefaultSettings(mesh3D& mesh3d) {

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
    mesh3d.bendYoungModulus = 1e6;
    mesh3d.strainRate = 1e2;

    mesh3d.Newton_Solver_Threshold = 1e-2;

    mesh3d.lengthRateLame = mesh3d.YoungModulus / (2 * (1 + mesh3d.PoissonRate));
    mesh3d.volumeRateLame = mesh3d.YoungModulus * mesh3d.PoissonRate / ((1 + mesh3d.PoissonRate) * (1 - 2 * mesh3d.PoissonRate));
    mesh3d.lengthRate = 4 * mesh3d.lengthRateLame / 3;
    mesh3d.volumeRate = mesh3d.volumeRateLame + 5 * mesh3d.lengthRateLame / 6;
    mesh3d.stretchStiffness = mesh3d.clothYoungModulus / (2 * (1 + mesh3d.PoissonRate));
    mesh3d.shearStiffness = mesh3d.stretchStiffness * 1;
    mesh3d.bendingStiffness = mesh3d.bendYoungModulus * pow(mesh3d.clothThicness, 3) / (24 * (1 - mesh3d.PoissonRate * mesh3d.PoissonRate));

}


void LoadSettings(mesh3D& mesh3d) {
    bool successfulRead = false;

    //read file
    std::ifstream infile;
    string asset_dir = string{ CIPC_ASSETS_DIR };

    string DEFAULT_CONFIG_FILE = asset_dir+"scene/parameterSetting.txt";

    infile.open(DEFAULT_CONFIG_FILE, std::ifstream::in);
    if (successfulRead = infile.is_open())
    {
        char ignoreToken[256];
        mesh3d.Fhat = 1e-6;
        mesh3d.Kappa = 0;
        mesh3d.dTol = 1e-18;

        // global settings:
        infile >> ignoreToken >> mesh3d.density;
        infile >> ignoreToken >> mesh3d.YoungModulus;
        infile >> ignoreToken >> mesh3d.PoissonRate;
        infile >> ignoreToken >> mesh3d.friction;
        infile >> ignoreToken >> mesh3d.clothThicness;
        infile >> ignoreToken >> mesh3d.clothYoungModulus;
        infile >> ignoreToken >> mesh3d.bendYoungModulus;
        infile >> ignoreToken >> mesh3d.cloth_density;
        infile >> ignoreToken >> mesh3d.strainRate;
        infile >> ignoreToken >> mesh3d.use_barrier;
        infile >> ignoreToken >> mesh3d.drag_coeff;
        infile >> ignoreToken >> mesh3d.IPC_dt;
        infile >> ignoreToken >> mesh3d.Newton_Solver_Threshold;
        infile >> ignoreToken >> mesh3d.Hhat;

        mesh3d.Hhat *= mesh3d.Hhat;
        mesh3d.lengthRateLame = mesh3d.YoungModulus / (2 * (1 + mesh3d.PoissonRate));
        mesh3d.volumeRateLame = mesh3d.YoungModulus * mesh3d.PoissonRate / ((1 + mesh3d.PoissonRate) * (1 - 2 * mesh3d.PoissonRate));
        mesh3d.lengthRate = 4 * mesh3d.lengthRateLame / 3;
        mesh3d.volumeRate = mesh3d.volumeRateLame + 5 * mesh3d.lengthRateLame / 6;
        mesh3d.stretchStiffness = mesh3d.clothYoungModulus / (2 * (1 + mesh3d.PoissonRate));
        mesh3d.shearStiffness = mesh3d.stretchStiffness * 0.3;
        mesh3d.bendingStiffness = mesh3d.bendYoungModulus * pow(mesh3d.clothThicness, 3) / (36 * (1 - mesh3d.PoissonRate * mesh3d.PoissonRate));

        cout << mesh3d.bendingStiffness << endl;

        infile.close();
    }

    if (!successfulRead)
    {
        std::cerr << "Waning: failed loading settings, set to defaults." << std::endl;
        DefaultSettings(mesh3d);
    }
}

void update_material(mesh3D& mesh3d) {
    mesh3d.lengthRateLame = mesh3d.YoungModulus / (2 * (1 + mesh3d.PoissonRate));
    mesh3d.volumeRateLame = mesh3d.YoungModulus * mesh3d.PoissonRate / ((1 + mesh3d.PoissonRate) * (1 - 2 * mesh3d.PoissonRate));
    mesh3d.lengthRate = 4 * mesh3d.lengthRateLame / 3;
    mesh3d.volumeRate = mesh3d.volumeRateLame + 5 * mesh3d.lengthRateLame / 6;
    mesh3d.stretchStiffness = mesh3d.clothYoungModulus / (2 * (1 + mesh3d.PoissonRate));
    mesh3d.shearStiffness = mesh3d.stretchStiffness * 0.3;
    mesh3d.bendingStiffness = mesh3d.bendYoungModulus * pow(mesh3d.clothThicness, 3) / (36 * (1 - mesh3d.PoissonRate * mesh3d.PoissonRate));

}

void case1(mesh3D& mesh3d, string asset_dir) {

    mesh3d.YoungModulus = 1e6;
    update_material(mesh3d);
    mesh3d.objMaxConer = Vector3d(0, 0, 0);
    mesh3d.objMinConer = Vector3d(0, 0, 0);
	mesh3d.is_quasi_static = true;
	mesh3d.apply_gravity = false;
    mesh3d.load_tetrahedraMesh(asset_dir + "tetrahedraMesh/ipcmesh/mat40x40.msh", 1, Vector3d(0, -0, 0));

    const double eps = 1e-4;
    for (int i = 0; i < mesh3d.vertexes.size(); i++)
    {
        if (mesh3d.vertexes[i][0] < mesh3d.objMinConer[0] + eps || mesh3d.vertexes[i][0] > mesh3d.objMaxConer[0] - eps)
        {
            mesh3d.boundary_vertexes_indices.push_back(i);
            mesh3d.boundaryTypes[i] = 1;
        }
    }

    mesh3d.update_hard_constraint_functor =
        [](Vector3d vertex, double alpha, double ipc_dt) -> Vector3d
        {
            double angleX = 3.14 / 5 * ipc_dt * alpha;
            Matrix3d rotationL, rotationR;
            rotationL << 1, 0, 0, 0, cos(angleX), sin(angleX), 0, -sin(angleX), cos(angleX);
            rotationR << 1, 0, 0, 0, cos(angleX), -sin(angleX), 0, sin(angleX), cos(angleX);
            double mvl = -1 * ipc_dt * alpha;
            Vector3d moveDir = Vector3d(0, 0, 0);
            if (vertex[0] < 0)
            {
                moveDir = rotationL * vertex - vertex;
            }
            if (vertex[0] > 0)
            {
                // rotate along x axis counterclockwise
                moveDir = rotationR * vertex - vertex;
            }
            return moveDir;
        };
}

void case2(mesh3D& mesh3d, string asset_dir) {

    mesh3d.YoungModulus = 1e4;
    update_material(mesh3d);
    mesh3d.objMaxConer = Vector3d(0, 0, 0);
    mesh3d.objMinConer = Vector3d(0, 0, 0);
    //mesh3d.is_quasi_static = true;
    if (true)
    {
        mesh3d.load_triangleMesh(asset_dir + "triangleMesh/CMU/plane100.obj", 1, Vector3d(0, 0, 0));

        PrepQuadBending(mesh3d);

        Matrix3d rotate, rotatey;
        float angleX = FEM::PI / 2, angleY = FEM::PI / 4, angleZ = FEM::PI / 2;
        rotate << 1, 0, 0, 0, cos(angleX), -sin(angleX), 0, sin(angleX), cos(angleX);
        rotatey << cos(angleY), -sin(angleY), 0, sin(angleY), cos(angleY), 0, 0, 0, 1;
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
        mesh3d.maxCorner = Vector3d(xmax, ymax, zmax);
        mesh3d.minCorner = Vector3d(xmin, ymin, zmin);

        //mesh3d.boundaryTypes[0] = 1;
        //mesh3d.Constraints[0].setZero();
    }

    mesh3d.load_tetrahedraMesh(asset_dir + "tetrahedraMesh/bunny.msh", 0.5, Vector3d(0, -0.5, 0));

}

bool FEMSimulator::buildModels(unsigned int buildType, unsigned int sceneType) {
    mesh3D mesh3d;
    LoadSettings(mesh3d);
    string asset_dir = string{ CIPC_ASSETS_DIR };
    mesh3d.maxCorner = Vector3d(-1e32, -1e32, -1e32);
    mesh3d.minCorner = Vector3d(1e32, 1e32, 1e32);

	case2(mesh3d, asset_dir);

    mesh3d.vertexNum = mesh3d.vertexes.size();
    mesh3d.tetrahedraNum = mesh3d.tetrahedras.size();
    mesh3d.triangleNum = mesh3d.triangles.size();

    printf("tets num: %d\n", mesh3d.tetrahedraNum);

    printf("verts num: %d\n", mesh3d.vertexNum);

    initMesh3D(mesh3d, 1, 0.2);

    mesh3d.v_rest = mesh3d.vertexes;
    mesh3d.V_prev = mesh3d.vertexes;

    computeXTilta(mesh3d);


    mesh3d.bboxDiagSize2 = (mesh3d.maxCorner - mesh3d.minCorner).squaredNorm();
    mesh3d.Hhat *= mesh3d.bboxDiagSize2;
    mesh3d.Fhat *= mesh3d.bboxDiagSize2;;
    mesh3d.dTol *= mesh3d.bboxDiagSize2;


    tetrahedra_meshes.mesh3Ds.push_back(mesh3d);
    tetrahedra_meshes.calculate_surface();

    double length = 0;
    for (const auto& edg : (tetrahedra_meshes.mesh3Ds[0]).surfEdges) {
        length += (tetrahedra_meshes.mesh3Ds[0].vertexes[edg.first] -
            tetrahedra_meshes.mesh3Ds[0].vertexes[edg.second]).norm();
    }
    tetrahedra_meshes.mesh3Ds[0].averageEdgeLength = length / (tetrahedra_meshes.mesh3Ds[0].surfEdges.size() * 3);
    printf("triangles num: %d\n", tetrahedra_meshes.mesh3Ds[0].surface.size());
    buildIntegrator(buildType, sceneType);

    if (!tetrahedra_meshes.mesh3Ds[0].load_tetTempData()) {
        printf("no temp data\n");
    }
    else {
        printf("load temp data\n");
    }

    buildCollisionSets();

    return true;
}


void FEMSimulator::buildIntegrator(const int integratorType, unsigned int sceneType) {
    integrator = std::make_unique<ImplicitFEMIntegrator>(&tetrahedra_meshes, sceneType);
}

void FEMSimulator::buildCollisionSets() {
    if (tetrahedra_meshes.mesh3Ds[0].use_barrier) {
        sh.build(tetrahedra_meshes.mesh3Ds[0], tetrahedra_meshes.mesh3Ds[0].averageEdgeLength);
        sh.calculateActivateSet(tetrahedra_meshes.mesh3Ds[0]);
    }
    
    gd.calculateActivateSet(tetrahedra_meshes.mesh3Ds[0]);
}

int FEMSimulator::simulateStick(int &stepId) {
    int cg_loops = 0;
    int newTon_loops = 0;
    int k = integrator->integrate(stepId, cg_loops, newTon_loops, sh, gd);
    //buildCollisionSets();
    return k;
}