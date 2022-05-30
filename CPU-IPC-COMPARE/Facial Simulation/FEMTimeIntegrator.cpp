#include "FEMTimeIntegrator.h"
#include "fem_parameters.h"
//#include "fem3D_Device.cuh"



ImplicitFEMIntegrator::ImplicitFEMIntegrator(model_tet* tetra_mesh3d, unsigned int m_sceneType) {
    meshTetes = tetra_mesh3d;
    sceneType = m_sceneType;
}



int ImplicitFEMIntegrator::integrate(double& mfsum, int& total_cg_iterations, int& total_newton_iterations, SpatialHash& sh, Ground& gd) {
    return IPC_Solver(meshTetes, sh, gd);
}






