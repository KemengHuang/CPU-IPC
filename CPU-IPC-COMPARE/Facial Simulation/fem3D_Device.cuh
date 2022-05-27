#pragma once
#ifndef __FEM3D_DEVICE_H
#define __FEM3D_DEVICE_H

#include "FEMMeshes.cuh"
#include "Device_PCG.cuh"

void Modify_Geometry(Tetrahedra_Data* mesh, const __GEIGEN__::Matrix3x3d* rotateMat, int modifiedNum);
void GPU_Projected_Newton3D(Tetrahedra_Data* mesh, PCG_Data* pcg_data, int vertexNum, int tetrahedraNum, double* fsum);

#endif //__FEM3D_DEVICE_H