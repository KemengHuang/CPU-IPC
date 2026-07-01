#pragma once
#ifndef FEM_SIMULATOR_H
#define FEM_SIMULATOR_H

//#include "setting.h"
//#include "mesh/Mesh.h"
#include "collision/BroadPhase.h"
//#include "FEMMeshes.cuh"
#include "core/TimeIntegrator.h"

class FEMSimulator {
public:
	FEMSimulator() { dt = 0.025; }
	~FEMSimulator() {}
	bool buildModels(unsigned int buildType, unsigned int sceneType);
	
	
	model_obj& getTriangleMeshes() { return triangle_meshes; }
	model_tet& getTetrahedraMeshes() { return tetrahedra_meshes; }
	//fiber_obj& getFiberObj() { return fiberObj; }
	int simulateStick(int& stepId);
private:
	void buildIntegrator(const int integratorType, unsigned int sceneType);
	void buildCollisionSets();
	int vertexNum;
	int tetrahedraNum;
	model_obj triangle_meshes;
	model_tet tetrahedra_meshes;
	//fiber_obj fiberObj;
	BroadPhase sh;
	Ground gd;
	std::unique_ptr<FEMIntegrator>  integrator;
	double dt;
};


#endif //FEM_SIMULATOR_H