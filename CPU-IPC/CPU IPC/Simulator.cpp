#include "fem_parameters.h"
#include "Simulator.h"
using namespace FEM;
bool FEMSimulator::buildModels(unsigned int buildType, unsigned int sceneType) {
    mesh3D mesh3d;
    mesh3d.maxConer = Vector3d(-DBL_MAX, -DBL_MAX, -DBL_MAX);
    mesh3d.minConer = Vector3d(DBL_MAX, DBL_MAX, DBL_MAX);

    mesh3d.objMaxConer = Vector3d(0, 0, 0);
    mesh3d.objMinConer = Vector3d(0, 0, 0);


    //mesh3d.load_tetrahedraMesh_IPC_TetMesh("tetrahedraMesh/ipcmesh/cube.msh", 0.5, Vector3d(0.25, 0.0, 0.25));
    //mesh3d.load_tetrahedraMesh_IPC_TetMesh("tetrahedraMesh/ipcmesh/cube.msh", 0.5, Vector3d(0.25, 0.6, 0.25));
    //mesh3d.load_tetrahedraMesh_IPC_TetMesh("tetrahedraMesh/ipcmesh/cube.msh", 0.5, Vector3d(0.25, 0.8, 0.25));
    //mesh3d.load_tetrahedraMesh("tetrahedraMesh/cube15.msh", 0.5, Vector3d(-0.25, -0.95, -0.25));
    //mesh3d.load_tetrahedraMesh("tetrahedraMesh/ipcmesh/sqballTet_.msh", 1, Vector3d(0, 0.75, 0));
    //mesh3d.load_tetrahedraMesh("tetrahedraMesh/ipcmesh/sqballTet_.msh", 1, Vector3d(0, -0.3, 0));
    //mesh3d.load_tetrahedraMesh_IPC_TetMesh("tetrahedraMesh/ipcmesh/sqballTet_.msh", 1, Vector3d(0, -0.75, 0));

    //if(true)
    //{
        //mesh3d.load_tetrahedraMesh("tetrahedraMesh/bunny2.msh", 0.2, Vector3d(0, 0, 0));
        //mesh3d.load_tetrahedraMesh("tetrahedraMesh/bunny2.msh", 0.2, Vector3d(0.6, 0, 0));
    //    mesh3d.load_tetrahedraMesh("tetrahedraMesh/bunny2.msh", 0.2, Vector3d(-0.6, 0, 0));
    //    mesh3d.load_tetrahedraMesh("tetrahedraMesh/bunny2.msh", 0.2, Vector3d(0, 0, 0.4));
    //    mesh3d.load_tetrahedraMesh("tetrahedraMesh/bunny2.msh", 0.2, Vector3d(0, 0, -0.4));
    //    mesh3d.load_tetrahedraMesh("tetrahedraMesh/bunny2.msh", 0.2, Vector3d(0.6, 0, 0.4));
    //    mesh3d.load_tetrahedraMesh("tetrahedraMesh/bunny2.msh", 0.2, Vector3d(-0.6, 0, 0.4));
    //    mesh3d.load_tetrahedraMesh("tetrahedraMesh/bunny2.msh", 0.2, Vector3d(0.6, 0, -0.4));
    //    mesh3d.load_tetrahedraMesh("tetrahedraMesh/bunny2.msh", 0.2, Vector3d(-0.6, 0, -0.4));

    //    mesh3d.load_tetrahedraMesh("tetrahedraMesh/bunny2.msh", 0.2, Vector3d(0, 0.65, 0));
    //    mesh3d.load_tetrahedraMesh("tetrahedraMesh/bunny2.msh", 0.2, Vector3d(0.6, 0.65, 0));
    //    mesh3d.load_tetrahedraMesh("tetrahedraMesh/bunny2.msh", 0.2, Vector3d(-0.6, 0.65, 0));
    //    mesh3d.load_tetrahedraMesh("tetrahedraMesh/bunny2.msh", 0.2, Vector3d(0, 0.65, 0.4));
    //    mesh3d.load_tetrahedraMesh("tetrahedraMesh/bunny2.msh", 0.2, Vector3d(0, 0.65, -0.4));
    //    mesh3d.load_tetrahedraMesh("tetrahedraMesh/bunny2.msh", 0.2, Vector3d(0.6, 0.65, 0.4));
    //    mesh3d.load_tetrahedraMesh("tetrahedraMesh/bunny2.msh", 0.2, Vector3d(-0.6, 0.65, 0.4));
    //    mesh3d.load_tetrahedraMesh("tetrahedraMesh/bunny2.msh", 0.2, Vector3d(0.6, 0.65, -0.4));
    //    mesh3d.load_tetrahedraMesh("tetrahedraMesh/bunny2.msh", 0.2, Vector3d(-0.6, 0.65, -0.4));


    //    mesh3d.load_tetrahedraMesh("tetrahedraMesh/bunny2.msh", 0.2, Vector3d(0, -0.65, 0));
    //    mesh3d.load_tetrahedraMesh("tetrahedraMesh/bunny2.msh", 0.2, Vector3d(0.6, -0.65, 0));
    //    mesh3d.load_tetrahedraMesh("tetrahedraMesh/bunny2.msh", 0.2, Vector3d(-0.6, -0.65, 0));
    //    mesh3d.load_tetrahedraMesh("tetrahedraMesh/bunny2.msh", 0.2, Vector3d(0, -0.65, 0.4));
    //    mesh3d.load_tetrahedraMesh("tetrahedraMesh/bunny2.msh", 0.2, Vector3d(0, -0.65, -0.4));
    //    mesh3d.load_tetrahedraMesh("tetrahedraMesh/bunny2.msh", 0.2, Vector3d(0.6, -0.65, 0.4));
    //    mesh3d.load_tetrahedraMesh("tetrahedraMesh/bunny2.msh", 0.2, Vector3d(-0.6, -0.65, 0.4));
    //    mesh3d.load_tetrahedraMesh("tetrahedraMesh/bunny2.msh", 0.2, Vector3d(0.6, -0.65, -0.4));
    //    mesh3d.load_tetrahedraMesh("tetrahedraMesh/bunny2.msh", 0.2, Vector3d(-0.6, -0.65, -0.4));
    //}
    mesh3d.load_tetrahedraMesh("tetrahedraMesh/twoBunny.msh", 2, Vector3d(-0, -0, -0));
    //mesh3d.load_tetrahedraMesh("tetrahedraMesh/bunny.msh", 0.5, Vector3d(-0, 0.4, -0));
    //mesh3d.load_tetrahedraMesh("tetrahedraMesh/bunny.msh", 0.5, Vector3d(-0, -0.4, -0));
    //{
    //    //mesh3d.load_tetrahedraMesh_IPC_TetMesh("tetrahedraMesh/ipcmesh/rod300x33.msh", 1, Vector3d(0, 0.1, 0));
    //    //mesh3d.load_tetrahedraMesh_IPC_TetMesh("tetrahedraMesh/ipcmesh/rod300x33.msh", 1, Vector3d(0, -0.1, 0));
    //    //mesh3d.load_tetrahedraMesh_IPC_TetMesh("tetrahedraMesh/ipcmesh/rod300x33.msh", 1, Vector3d(0, 0, 0.1));
    //    //mesh3d.load_tetrahedraMesh_IPC_TetMesh("tetrahedraMesh/ipcmesh/rod300x33.msh", 1, Vector3d(0, 0, -0.1));

    //    mesh3d.load_tetrahedraMesh_IPC_TetMesh("tetrahedraMesh/ipcmesh/mat150x150t40.msh", 1, Vector3d(0, 0, 0));
    //    for (int j = 0; j < mesh3d.vertexNum; j++) {
    //        if ((mesh3d.vertexes[j].x()) > 0.5 - 1e-4) {
    //            mesh3d.boundaryTypes[j] = 1;
    //            mesh3d.Constraints[j].setZero();
    //        }
    //        if (((mesh3d.vertexes[j].x()) < -0.5 + 1e-4)) {
    //            mesh3d.boundaryTypes[j] = -1;
    //            mesh3d.Constraints[j].setZero();
    //        }
    //    }
    //}
    //for (int j = mesh3d.vertexNum/2; j < mesh3d.vertexNum; j++) {
    //    mesh3d.boundaryTypes[j] = 1;
    //    mesh3d.Constraints[j].setZero();
    //}
    mesh3d.vertexNum = mesh3d.vertexes.size();
    mesh3d.tetrahedraNum = mesh3d.tetrahedras.size();
    printf("tets num: %d\n", mesh3d.tetrahedraNum);
    printf("verts num: %d\n", mesh3d.vertexNum);

    initMesh3D(mesh3d, 1, 0.2);
    //for (int i = 0; i < mesh3d.tetrahedraNum; i++) {
    //    Vector4i tet = Vector4i(mesh3d.tetrahedras[i][0], mesh3d.tetrahedras[i][1], mesh3d.tetrahedras[i][2], mesh3d.tetrahedras[i][3]);
    //    mesh3d.tetrahedrasV.push_back(tet);
    //}

    //mesh3d.restSNKE = getObjRestEnergy_StableNHK2_3D(mesh3d.vertexes, mesh3d, FEM::lengthRate, FEM::volumeRate) * FEM::IPC_dt * FEM::IPC_dt;

    mesh3d.v_rest = mesh3d.vertexes;
    mesh3d.V_prev = mesh3d.vertexes;

    computeXTilta(mesh3d);


    mesh3d.bboxDiagSize2 = (mesh3d.maxConer - mesh3d.minConer).squaredNorm();
    mesh3d.Hhat = 1e-6 * mesh3d.bboxDiagSize2;// (mesh3d.objMaxConer - mesh3d.objMinConer).squaredNorm();
    mesh3d.Kappa = 0;
    mesh3d.dTol = 1e-18 * mesh3d.bboxDiagSize2;
    tetrahedra_meshes.mesh3Ds.push_back(mesh3d);
    tetrahedra_meshes.calculate_surface();

    double length = 0;
    for (const auto& edg : (tetrahedra_meshes.mesh3Ds[0]).surfEdges) {
        length += (tetrahedra_meshes.mesh3Ds[0].vertexes[edg.first] - tetrahedra_meshes.mesh3Ds[0].vertexes[edg.second]).norm();
    }
    tetrahedra_meshes.mesh3Ds[0].averageEdgeLenth = length / (tetrahedra_meshes.mesh3Ds[0].surfEdges.size() * 3);

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
    sh.build(tetrahedra_meshes.mesh3Ds[0], tetrahedra_meshes.mesh3Ds[0].averageEdgeLenth);
    gd.calculateActivateSet(tetrahedra_meshes.mesh3Ds[0]);
    sh.calculateActivateSet(tetrahedra_meshes.mesh3Ds[0]);
}

int FEMSimulator::simulateStick(int& stepId) {
    int cg_loops = 0;
    int newTon_loops = 0;
    int k =  integrator->integrate(stepId, cg_loops, newTon_loops, sh, gd);
    //buildCollisionSets();
    return k;
}