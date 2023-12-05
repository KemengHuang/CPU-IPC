#include "fem_parameters.h"
#include "Simulator.h"
#include "iostream"
#include "Eigen/Eigen"
using namespace Eigen;
using namespace FEM;

void buildSpecialPoints(mesh3D &mesh) {
    Vector3d tp[6];
    tp[0] = Vector3d(-0.26554,0.67295,-0.23324);
    tp[1] = Vector3d(-0.05724, 0.69143, -0.10306);
    tp[2] = Vector3d(0.04884, 0.47191, 0.04952);
    tp[3] = Vector3d(-0.26145, -0.12829, -0.23713);
    tp[4] = Vector3d(-0.06139, -0.10531, -0.09838);
    tp[5] = Vector3d(0.015541, -0.32489, 0.02275);
    mesh.spectialPontsArray.resize(6);
    for(int i=0;i<mesh.spectialPontsArray.size();i++){
        double mindist = 1e32;
        for(int j=0;j<mesh.surfVerts.size();j++){
            double dist = (mesh.vertexes[mesh.surfVerts[j]]-tp[i]).norm();
            //printf("dist: %f\n", dist);
            if(dist<mindist){
                mindist = dist;
                mesh.spectialPontsArray[i] = mesh.surfVerts[j];
                //printf("special Point: %d\n", mesh.surfVerts[j]);
            }
        }
        //printf("special Point: %d\n", mesh.spectialPontsArray[i]);
        //cout<<mesh.vertexes[mesh.spectialPontsArray[i]]<<endl;
    }


}


bool FEMSimulator::buildModels(unsigned int buildType, unsigned int sceneType) {
    mesh3D mesh3d;
    mesh3d.maxConer = Vector3d(-1e32, -1e32, -1e32);
    mesh3d.minConer = Vector3d(1e32, 1e32, 1e32);

    mesh3d.objMaxConer = Vector3d(0, 0, 0);
    mesh3d.objMinConer = Vector3d(0, 0, 0);

    if(false)
    {
        mesh3d.load_tetrahedraMesh_IPC_TetMesh("../CPU IPC/tetrahedraMesh/ipcmesh/dolphin5K.msh", 0.01, Vector3d(0, 0, 0));

        Matrix3d rotate;
        float angleZ = FEM::PI / 4;
        rotate << cos(angleZ), -sin(angleZ), 0,  sin(angleZ), cos(angleZ), 0,0,0,1;

        for (int j = 0; j < mesh3d.vertexNum; j++) {
            mesh3d.vertexes[j] = (rotate * mesh3d.vertexes[j]);
            mesh3d.vertexes[j] = (mesh3d.vertexes[j] + Vector3d(-0.7, -0.45, -0.525));
        }
        int tnum = mesh3d.vertexNum;

        mesh3d.load_triangleMesh("../CPU IPC/triangleMesh/newtubing.obj", 0.05, Vector3d(0, 0, 0), 2);

        angleZ = FEM::PI / 2;
        Matrix3d rotate2;// << cos(angleZ), -sin(angleZ), 0,  sin(angleZ), cos(angleZ), 0,0,0,1;
        rotate2 << cos(angleZ), -sin(angleZ), 0,  sin(angleZ), cos(angleZ), 0,0,0,1;

        for (int j = tnum; j < mesh3d.vertexNum; j++) {
            mesh3d.vertexes[j] = (rotate2 * mesh3d.vertexes[j]);
            mesh3d.vertexes[j] = (mesh3d.vertexes[j] + Vector3d(0.5, .0, 0));
        }


        double xmin = 1e32, ymin = 1e32, zmin = 1e32;
        double xmax = -1e32, ymax = -1e32, zmax = -1e32;
        for (int j = 0; j < mesh3d.vertexNum; j++) {
            double dis = pow(mesh3d.vertexes[j].y(),2)+pow(mesh3d.vertexes[j].z(),2);
            if (dis<0.012*0.012&&mesh3d.vertexes[j].x()>-0.25&&mesh3d.vertexes[j].x()<-0.07) {
                mesh3d.boundaryTypes[j] = 1;
                mesh3d.Constraints[j].setZero();
            }

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


    if (false)
    {
        mesh3d.load_tetrahedraMesh_IPC_TetMesh("../CPU IPC/tetrahedraMesh/ipcmesh/Armadillo13K.msh", 0.006, Vector3d(0, 0, 0));

        Matrix3d rotate;
        float angleX = FEM::PI / 2, angleY = -FEM::PI / 2, angleZ = FEM::PI / 2;
        rotate << cos(angleY), 0, -sin(angleY), 0, 1, 0, sin(angleY), 0, cos(angleY);

        for (int j = 0; j < mesh3d.vertexNum; j++) {
            mesh3d.vertexes[j] = (rotate * mesh3d.vertexes[j]);

            mesh3d.vertexes[j] = (mesh3d.vertexes[j] + Vector3d(-0.1, 0.4, 0));
        }
        mesh3d.load_tetrahedraMesh("../CPU IPC/tetrahedraMesh/bunny2.msh", 0.2, Vector3d(0, -0.6, 0));
        mesh3d.load_triangleMesh("../CPU IPC/triangleMesh/tricloth.obj", 0.7, Vector3d(0, 0, 0));
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

    if(false)
    {
        //mesh3d.load_triangleMesh("../CPU IPC/triangleMesh/tricloth.obj", 0.3, Vector3d(0, 0, 0));
        mesh3d.load_triangleMesh("../CPU IPC/triangleMesh/CMU/plane1024.obj", 0.3, Vector3d(0, 0, 0));

        Matrix3d rotate,rotatey;
        float angleX = FEM::PI / 2, angleY = FEM::PI / 4, angleZ = FEM::PI / 2;
        rotate << 1,0,0,0,cos(angleX), -sin(angleX), 0, sin(angleX), cos(angleX);
        rotatey<<cos(angleY), -sin(angleY), 0,sin(angleY), cos(angleY),0,0,0,1;
        for (int j = 0; j < mesh3d.vertexNum; j++) {
            //mesh3d.vertexes[j] = (rotate * mesh3d.vertexes[j]);

            //mesh3d.vertexes[j] = (rotatey * mesh3d.vertexes[j]);
        }
        //mesh3d.load_tetrahedraMesh("../CPU IPC/tetrahedraMesh/bunny2.msh", 0.2, Vector3d(0, 0, 0));
        //mesh3d.load_triangleMesh("../CPU IPC/triangleMesh/tricloth.obj", 0.7, Vector3d(0, 0, 0));
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



    
    //mesh3d.load_tetrahedraMesh_IPC_TetMesh("../CPU IPC/tetrahedraMesh/ipcmesh/cube.msh", 0.5, Vector3d(0.25, -0.6, 0.25));
    //mesh3d.load_tetrahedraMesh("../CPU IPC/tetrahedraMesh/bunny2.msh", 0.2, Vector3d(0, -0.65, 0));
    //mesh3d.load_tetrahedraMesh("../CPU IPC/tetrahedraMesh/bunny2.msh", 0.2, Vector3d(0, 0.15, 0));
    //mesh3d.load_tetrahedraMesh("../CPU IPC/tetrahedraMesh/bunny2.msh", 0.2, Vector3d(0, 0.2, 0));
    //mesh3d.load_tetrahedraMesh_IPC_TetMesh("tetrahedraMesh/ipcmesh/cube.msh", 0.5, Vector3d(0.25, 0.6, 0.25));
    //mesh3d.load_tetrahedraMesh_IPC_TetMesh("tetrahedraMesh/ipcmesh/cube.msh", 0.5, Vector3d(0.25, 0.8, 0.25));
    //mesh3d.load_tetrahedraMesh("tetrahedraMesh/cube15.msh", 0.5, Vector3d(-0.25, -0.95, -0.25));
    //mesh3d.load_tetrahedraMesh("tetrahedraMesh/ipcmesh/sqballTet_.msh", 1, Vector3d(0, 0.75, 0));
    //mesh3d.load_tetrahedraMesh("tetrahedraMesh/ipcmesh/sqballTet_.msh", 1, Vector3d(0, -0.3, 0));
    //mesh3d.load_tetrahedraMesh_IPC_TetMesh("tetrahedraMesh/ipcmesh/sqballTet_.msh", 1, Vector3d(0, -0.75, 0));

    if(false)
    {
        mesh3d.load_tetrahedraMesh("../CPU IPC/tetrahedraMesh/bunny2.msh", 0.2, Vector3d(0, 0, 0));
    //mesh3d.load_triangleMesh("../CPU IPC/triangleMesh/tricloth.obj", 0.7, Vector3d(0, 0, 0));
        mesh3d.load_tetrahedraMesh("../CPU IPC/tetrahedraMesh/bunny2.msh", 0.2, Vector3d(0.6, 0, 0));
       mesh3d.load_tetrahedraMesh("../CPU IPC/tetrahedraMesh/bunny2.msh", 0.2, Vector3d(-0.6, 0, 0));
       mesh3d.load_tetrahedraMesh("../CPU IPC/tetrahedraMesh/bunny2.msh", 0.2, Vector3d(0, 0, 0.4));
       mesh3d.load_tetrahedraMesh("../CPU IPC/tetrahedraMesh/bunny2.msh", 0.2, Vector3d(0, 0, -0.4));
       mesh3d.load_tetrahedraMesh("../CPU IPC/tetrahedraMesh/bunny2.msh", 0.2, Vector3d(0.6, 0, 0.4));
       mesh3d.load_tetrahedraMesh("../CPU IPC/tetrahedraMesh/bunny2.msh", 0.2, Vector3d(-0.6, 0, 0.4));
       mesh3d.load_tetrahedraMesh("../CPU IPC/tetrahedraMesh/bunny2.msh", 0.2, Vector3d(0.6, 0, -0.4));
       mesh3d.load_tetrahedraMesh("../CPU IPC/tetrahedraMesh/bunny2.msh", 0.2, Vector3d(-0.6, 0, -0.4));

     //mesh3d.load_tetrahedraMesh("/home/mine/Codes/CPU_IPC_BENCHMARK/CPU-IPC/tetrahedraMesh/twoBunny.msh", 0.2, Vector3d(0, 0.65, 0));
        mesh3d.load_tetrahedraMesh("../CPU IPC/tetrahedraMesh/bunny2.msh", 0.2, Vector3d(0, 0.65, 0));
    //mesh3d.load_tetrahedraMesh_IPC_TetMesh("/home/mine/Codes/CPU_IPC_BENCHMARK/CPU-IPC/tetrahedraMesh/sphere1K.vtk", 0.2, Vector3d(0, 0.65, 0));
       mesh3d.load_tetrahedraMesh("../CPU IPC/tetrahedraMesh/bunny2.msh", 0.2, Vector3d(0.6, 0.65, 0));
       mesh3d.load_tetrahedraMesh("../CPU IPC/tetrahedraMesh/bunny2.msh", 0.2, Vector3d(-0.6, 0.65, 0));
       mesh3d.load_tetrahedraMesh("../CPU IPC/tetrahedraMesh/bunny2.msh", 0.2, Vector3d(0, 0.65, 0.4));
       mesh3d.load_tetrahedraMesh("../CPU IPC/tetrahedraMesh/bunny2.msh", 0.2, Vector3d(0, 0.65, -0.4));
       mesh3d.load_tetrahedraMesh("../CPU IPC/tetrahedraMesh/bunny2.msh", 0.2, Vector3d(0.6, 0.65, 0.4));
       mesh3d.load_tetrahedraMesh("../CPU IPC/tetrahedraMesh/bunny2.msh", 0.2, Vector3d(-0.6, 0.65, 0.4));
       mesh3d.load_tetrahedraMesh("../CPU IPC/tetrahedraMesh/bunny2.msh", 0.2, Vector3d(0.6, 0.65, -0.4));
       mesh3d.load_tetrahedraMesh("../CPU IPC/tetrahedraMesh/bunny2.msh", 0.2, Vector3d(-0.6, 0.65, -0.4));


       mesh3d.load_tetrahedraMesh("../CPU IPC/tetrahedraMesh/bunny2.msh", 0.2, Vector3d(0, -0.65, 0));
       mesh3d.load_tetrahedraMesh("../CPU IPC/tetrahedraMesh/bunny2.msh", 0.2, Vector3d(0.6, -0.65, 0));
       mesh3d.load_tetrahedraMesh("../CPU IPC/tetrahedraMesh/bunny2.msh", 0.2, Vector3d(-0.6, -0.65, 0));
       mesh3d.load_tetrahedraMesh("../CPU IPC/tetrahedraMesh/bunny2.msh", 0.2, Vector3d(0, -0.65, 0.4));
       mesh3d.load_tetrahedraMesh("../CPU IPC/tetrahedraMesh/bunny2.msh", 0.2, Vector3d(0, -0.65, -0.4));
       mesh3d.load_tetrahedraMesh("../CPU IPC/tetrahedraMesh/bunny2.msh", 0.2, Vector3d(0.6, -0.65, 0.4));
       mesh3d.load_tetrahedraMesh("../CPU IPC/tetrahedraMesh/bunny2.msh", 0.2, Vector3d(-0.6, -0.65, 0.4));
       mesh3d.load_tetrahedraMesh("../CPU IPC/tetrahedraMesh/bunny2.msh", 0.2, Vector3d(0.6, -0.65, -0.4));
       mesh3d.load_tetrahedraMesh("../CPU IPC/tetrahedraMesh/bunny2.msh", 0.2, Vector3d(-0.6, -0.65, -0.4));
    }
    //mesh3d.load_tetrahedraMesh("tetrahedraMesh/twoBunny.msh", 2, Vector3d(-0, -0, -0));
    mesh3d.load_tetrahedraMesh("tetrahedraMesh/bunny.msh", 0.5, Vector3d(-0, 0.4, -0));
    mesh3d.load_tetrahedraMesh("tetrahedraMesh/bunny.msh", 0.5, Vector3d(-0, -0.4, -0));
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
    mesh3d.triangleNum = mesh3d.triangles.size();
    //mesh3d.triangleNum = 0;
    printf("tets num: %d\n", mesh3d.tetrahedraNum);
    //printf("triangles num: %d\n", mesh3d.triangleNum);
    printf("verts num: %d\n", mesh3d.vertexNum);

    initMesh3D(mesh3d, 1, 0.2);
    //mesh3d.triangleNum = 0;

    //for (int i = 0; i < mesh3d.tetrahedraNum; i++) {
    //    Vector4i tet = Vector4i(mesh3d.tetrahedras[i][0], mesh3d.tetrahedras[i][1], mesh3d.tetrahedras[i][2], mesh3d.tetrahedras[i][3]);
    //    mesh3d.tetrahedrasV.push_back(tet);
    //}

    //mesh3d.restSNKE = getObjRestEnergy_StableNHK2_3D(mesh3d.vertexes, mesh3d, FEM::lengthRate, FEM::volumeRate) * FEM::IPC_dt * FEM::IPC_dt;

    mesh3d.v_rest = mesh3d.vertexes;
    mesh3d.V_prev = mesh3d.vertexes;

    computeXTilta(mesh3d);


    mesh3d.bboxDiagSize2 = (mesh3d.maxConer - mesh3d.minConer).squaredNorm();
    mesh3d.Hhat = 9e-8 * mesh3d.bboxDiagSize2;// (mesh3d.objMaxConer - mesh3d.objMinConer).squaredNorm();
    mesh3d.Fhat = 1e-6 * mesh3d.bboxDiagSize2;;
    mesh3d.Kappa = 0;
    mesh3d.dTol = 1e-18 * mesh3d.bboxDiagSize2;



    tetrahedra_meshes.mesh3Ds.push_back(mesh3d);
    tetrahedra_meshes.calculate_surface();
    //tetrahedra_meshes.mesh3Ds[0].triangles.clear();
    //buildSpecialPoints(tetrahedra_meshes.mesh3Ds[0]);

    double length = 0;
    for (const auto &edg: (tetrahedra_meshes.mesh3Ds[0]).surfEdges) {
        length += (tetrahedra_meshes.mesh3Ds[0].vertexes[edg.first] -
                   tetrahedra_meshes.mesh3Ds[0].vertexes[edg.second]).norm();
    }
    tetrahedra_meshes.mesh3Ds[0].averageEdgeLenth = length / (tetrahedra_meshes.mesh3Ds[0].surfEdges.size() * 3);
    printf("triangles num: %d\n", tetrahedra_meshes.mesh3Ds[0].surface.size());
    buildIntegrator(buildType, sceneType);

    if (!tetrahedra_meshes.mesh3Ds[0].load_tetTempData()) {
        printf("no temp data\n");
    } else {
        printf("load temp data\n");
    }

    buildCollisionSets();

    //printf("finish build cp\n");

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

int FEMSimulator::simulateStick(int &stepId) {
    int cg_loops = 0;
    int newTon_loops = 0;
    int k = integrator->integrate(stepId, cg_loops, newTon_loops, sh, gd);
    //buildCollisionSets();
    return k;
}