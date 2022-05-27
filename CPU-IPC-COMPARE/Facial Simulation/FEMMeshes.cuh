#ifndef  __FEM_MESHES_CUH__
#define __FEM_MESHES_CUH__

#include <cuda_runtime.h>
#include "Eigen_Data.h"

class Tetrahedra_Data{
public:
	double3* vertexes;
	double3* temp_vertexes;
	double3* tetra_fiberDir;
	uint4* tetrahedras;
	double* squeue;
	double* volum;
	int* isInside;
	int* modified_ids;
	bool* isMuscle;
	__GEIGEN__::Matrix3x3d* DmInverses;
	__GEIGEN__::Matrix3x3d* Constraints;
	__GEIGEN__::Matrix12x12d* Hessians;
public:
	Tetrahedra_Data() {}
	~Tetrahedra_Data();
	void Malloc_DEVICE_MEM(const int& vertex_num, const int& tetradedra_num, const int& modified_num);
	void FREE_DEVICE_MEM();

};


#endif // ! __FEM_MESHES_CUH__
