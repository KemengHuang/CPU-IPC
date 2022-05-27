#ifndef  __DEVICE_PCG_CUH__
#define __DEVICE_PCG_CUH__
#include <cuda_runtime.h>

class PCG_Data {
public:
	double* squeue;
	double3* b;
	double3* P;
	double3* r;
	double3* c;
	double3* q;
	double3* s;
	double3* dx;
	double3* tempDx;
public:
	void Malloc_DEVICE_MEM(const int& vertex_num, const int& tetradedra_num);
	void FREE_DEVICE_MEM();
};





#endif //__DEVICE_PCG_CUH__