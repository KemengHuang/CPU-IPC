#include "fem_parameters.h"
#include "Simulator.h"
using namespace FEM;
bool FEMSimulator::buildModels(unsigned int buildType, unsigned int sceneType) {
    if (sceneType == 0) {
        mesh3D mesh3d;
        mesh3d.maxConer = Vector3d(-DBL_MAX, -DBL_MAX, -DBL_MAX);
        mesh3d.minConer = Vector3d(DBL_MAX, DBL_MAX, DBL_MAX);

        mesh3d.objMaxConer = Vector3d(0, 0, 0);
        mesh3d.objMinConer = Vector3d(0, 0, 0);

        //mesh3d.load_tetrahedraMesh_IPC_TetMesh("tetrahedraMesh/ipcmesh/ANSYS1e-3Corner_coarse.msh", Vector3d(10, 10, 10), Vector3d(0, 0.0, 0));
        ////mesh3d.load_tetrahedraMesh_IPC_TetMesh("tetrahedraMesh/ipcmesh/cube.msh", Vector3d(0.1, 0.1, 0.1), Vector3d(0, 0.5, 0));

        //Matrix3d rotate;
        //double theta = PI / 4;
        //rotate << cos(theta), -sin(theta), 0, sin(theta), cos(theta), 0, 0, 0, 1;

        //for (int j = 0; j < mesh3d.vertexNum; j++) {
        //    if (mesh3d.vertexes[j].y() < 0.105) {
        //        mesh3d.boundaryTypes[j] = 1;
        //        mesh3d.Constraints[j].setZero();
        //    }
        //    else {
        //        mesh3d.vertexes[j].y() += 0.3;
        //        //mesh3d.boundaryTypes[j] = 1;
        //        //mesh3d.Constraints[j].setZero();
        //    }
        //}

       

        //mesh3d.load_tetrahedraMesh_IPC_TetMesh("tetrahedraMesh/ipcmesh/cube.msh", Vector3d(0.3, 0.3, 0.3), Vector3d(-3.7, -5.0, -4.5));
        //mesh3d.load_tetrahedraMesh_IPC_TetMesh("tetrahedraMesh/ipcmesh/cube.msh", Vector3d(10.20, 0.5, 10.20), Vector3d(4.25, 0.6, 4.25));
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
        //mesh3d.load_tetrahedraMesh("tetrahedraMesh/cubes0.msh", 1, Vector3d(-0, -0, -0));
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
        mesh3d.IPC_dt = 0.01;


        mesh3d.v_rest = mesh3d.vertexes;
        mesh3d.V_prev = mesh3d.vertexes;

        computeXTilta(mesh3d);


        mesh3d.bboxDiagSize2 = (mesh3d.maxConer - mesh3d.minConer).squaredNorm();
        mesh3d.Hhat = 1e-6 * mesh3d.bboxDiagSize2;
        mesh3d.Kappa = 0;
        mesh3d.dTol = 1e-18 * mesh3d.bboxDiagSize2;
        tetrahedra_meshes.mesh3Ds.push_back(mesh3d);
        tetrahedra_meshes.calculate_surface();

        mesh3d.restSNKE = getObjRestEnergy_StableNHK2_3D(mesh3d.vertexes, mesh3d, FEM::lengthRate, FEM::volumeRate) * mesh3d.IPC_dt * mesh3d.IPC_dt;

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
    }
    else if (sceneType == 1) {
        mesh3D mesh3d;
        mesh3d.maxConer = Vector3d(-DBL_MAX, -DBL_MAX, -DBL_MAX);
        mesh3d.minConer = Vector3d(DBL_MAX, DBL_MAX, DBL_MAX);

        mesh3d.objMaxConer = Vector3d(0, 0, 0);
        mesh3d.objMinConer = Vector3d(0, 0, 0);


        {
            //mesh3d.load_tetrahedraMesh_IPC_TetMesh("tetrahedraMesh/ipcmesh/ANSYS1e-3Corner_coarse.msh", Vector3d(10, 10, 10), Vector3d(0, 0.0, 0));
            //mesh3d.load_tetrahedraMesh_IPC_TetMesh("tetrahedraMesh/ipcmesh/cube.msh", Vector3d(0.1, 0.1, 0.1), Vector3d(0, 0.5, 0));

            Matrix3d rotate;
            double theta = PI / 1;
            rotate << cos(theta), -sin(theta), 0, sin(theta), cos(theta), 0, 0, 0, 1;

            Vector3d P0 = Vector3d(-sqrt(3), 1, -1);
            Vector3d P1 = Vector3d(sqrt(3), 1, -1);
            Vector3d P2 = Vector3d(0, 1, 2);
            Vector3d P3 = Vector3d(0, -1, 0);

            //Vector3d P0 = Vector3d(0, -1, 1);
            //Vector3d P1 = Vector3d(0, -1, -1);
            //Vector3d P2 = Vector3d(-1, 1, 0);
            //Vector3d P3 = Vector3d(1, 1, 0);

            vector<Vector3d> Parray;
            Parray.push_back(P0);
            Parray.push_back(P1);
            Parray.push_back(P2);
            Parray.push_back(P3);


            mesh3d.load_UnitTest(Parray, Vector3d(0.3, 0.3, 0.3), Vector3d(0, -0.6, 0));
            Vector3d center = Vector3d(0, 0, 0);//0.5 * (mesh3d.maxConer - mesh3d.minConer);

            //mesh3d.load_UnitTest2(Vector3d(0.6, 0.6, 0.6), Vector3d(0, 0.3, 0));
            //for (int j = 4; j < mesh3d.vertexNum; j++) {
            //    //if (mesh3d.vertexes[j].y() > 0.105) {
            //    //mesh3d.vertexes[j] -= center;
            //    mesh3d.boundaryTypes[j] = 1;
            //    mesh3d.Constraints[j].setZero();
            //    //}
            //}

            //P0 = Vector3d(0, 1, 1);
            //P1 = Vector3d(0, 1, -1);
            //P2 = Vector3d(-1, -1, 0);
            //P3 = Vector3d(1, -1, 0);

            P0 = Vector3d(-sqrt(3), 1, -1);
            P1 = Vector3d(sqrt(3), 1, -1);
            P2 = Vector3d(0, 1, 2);
            P3 = Vector3d(0, -1, 0);

            Parray.resize(0);
            Parray.push_back(P0);
            Parray.push_back(P1);
            Parray.push_back(P2);
            Parray.push_back(P3);

            Matrix3d rotate2;
            theta = PI / 2;
            rotate2 << cos(theta), 0, -sin(theta), 0, 1, 0, sin(theta), 0, cos(theta);
            rotate = rotate2 * rotate;
            mesh3d.load_UnitTest(Parray, Vector3d(0.3, 0.3, 0.3), Vector3d(0, 0.4, 0));
            //mesh3d.load_tetrahedraMesh_IPC_TetMesh("tetrahedraMesh/ipcmesh/cube.msh", Vector3d(0.3, 0.3, 0.3), Vector3d(0, 0, 0));
            for (int j = mesh3d.vertexNum / 2; j < mesh3d.vertexNum; j++) {
                //if (mesh3d.vertexes[j].y() > 0.105) {
                //mesh3d.vertexes[j] -= center;
                //mesh3d.vertexes[j] = rotate * mesh3d.vertexes[j];
                //mesh3d.vertexes[j] += center;
                double dx, dy, dz;
                dx = dy = dz = 0;
                int conner = 1;
                if (conner == 1) {
                    dx = 0.5 * sqrt(3) * 0.5;
                    dz = 0.25;
                }
                else if (conner == -1) {
                    dx = -0.5 * sqrt(3) * 0.5;
                    dz = 0.25;
                }
                else if (conner == 2) {
                    dx = 0;
                    dz = -0.5;
                }
                //mesh3d.vertexes[j].y() -= 0.3;
                mesh3d.vertexes[j].z() += dz;
                mesh3d.vertexes[j].x() += dx;
                mesh3d.boundaryTypes[j] = 1;
                mesh3d.Constraints[j].setZero();
                //}
            }
        }


        mesh3d.vertexNum = mesh3d.vertexes.size();
        mesh3d.tetrahedraNum = mesh3d.tetrahedras.size();
        printf("tets num: %d\n", mesh3d.tetrahedraNum);
        printf("verts num: %d\n", mesh3d.vertexNum);

        initMesh3D(mesh3d, 1, 0.2);

        mesh3d.IPC_dt = 0.01;


        mesh3d.v_rest = mesh3d.vertexes;
        mesh3d.V_prev = mesh3d.vertexes;

        computeXTilta(mesh3d);


        mesh3d.bboxDiagSize2 = (mesh3d.maxConer - mesh3d.minConer).squaredNorm();
        mesh3d.Hhat = 1e-6 * mesh3d.bboxDiagSize2;///*22.5e-8*/ */* mesh3d.bboxDiagSize2;*/ (mesh3d.objMaxConer - mesh3d.objMinConer).squaredNorm();
        mesh3d.Kappa = 0;
        mesh3d.dTol = 1e-18 * mesh3d.bboxDiagSize2;
        tetrahedra_meshes.mesh3Ds.push_back(mesh3d);
        tetrahedra_meshes.calculate_surface();

        mesh3d.restSNKE = getObjRestEnergy_StableNHK2_3D(mesh3d.vertexes, mesh3d, FEM::lengthRate, FEM::volumeRate) * mesh3d.IPC_dt * mesh3d.IPC_dt;

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
    }
    else {
        Vector3d A = Vector3d((1.0 / 200.0) * 1e-3, 1.0 / 2.0 * 1e-3, 0);
        Vector3d B = Vector3d(0, 0, 0);
        Vector3d C = Vector3d(1e-3, 0, 0);
        Vector3d D = Vector3d(5e-4, 0, sqrt(3) * 0.5e-3);
        mesh3D mesh3d;
        mesh3d.vertexes.push_back(A);
        mesh3d.vertexes.push_back(B);
        mesh3d.vertexes.push_back(C);
        mesh3d.vertexes.push_back(D);
        mesh3d.vertexes.push_back(Vector3d(0, 0, 0));
        mesh3d.vertexes.push_back(Vector3d(0, 0, 0));
        mesh3d.vertexes.push_back(Vector3d(0, 0, 0));
        mesh3d.vertexes.push_back(Vector3d(0, 0, 0));



        mesh3d.boundaryTypes.push_back(0);
        mesh3d.boundaryTypes.push_back(1);
        mesh3d.boundaryTypes.push_back(1);
        mesh3d.boundaryTypes.push_back(1);
        mesh3d.masses.push_back(1e0);
        mesh3d.masses.push_back(1e0);
        mesh3d.masses.push_back(1e0);
        mesh3d.masses.push_back(1e0);
        mesh3d.Hhat = 1e-6;
        mesh3d.Kappa = 2e6;

        Matrix3d Constraint1; Constraint1.setIdentity();
        Matrix3d Constraint2; Constraint2.setIdentity();//setIdentity();
        Matrix3d Constraint3; Constraint3.setIdentity();//setIdentity();
        Matrix3d Constraint4; Constraint4.setIdentity();
        mesh3d.Constraints.push_back(Constraint1);
        mesh3d.Constraints.push_back(Constraint2);
        mesh3d.Constraints.push_back(Constraint3);
        mesh3d.Constraints.push_back(Constraint4);



        tetrahedra_meshes.mesh3Ds.push_back(mesh3d);
        buildIntegrator(buildType, sceneType);
    }
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