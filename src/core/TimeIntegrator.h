#pragma once

#include "mesh/Mesh.h"
#include <memory>
#include "collision/BroadPhase.h"
#include "core/IPCSolver.h"
class FEMIntegrator {
public:
	FEMIntegrator() {};
	FEMIntegrator(int vertexNUm, int tetrahedraNum);
	~FEMIntegrator() {};
	virtual int integrate(int& stepId, int& total_cg_iterations, int& total_newton_iterations, BroadPhase& sh, Ground& gd) = 0;
public:
	int vertex_Num;
	int tetrahedra_Num;
	int modified_Num;
	int sceneType;
	int loopCount;
};



class ImplicitFEMIntegrator : public FEMIntegrator {
	friend class FEMSimulator;
public:
	ImplicitFEMIntegrator(model_tet* tetra_mesh3d, unsigned int sceneType);
	~ImplicitFEMIntegrator() {}
	virtual int integrate(int& stepId, int& total_cg_iterations, int& total_newton_iterations, BroadPhase& sh, Ground& gd);
private:

	//void Projected_Newton3D(mesh3D& mesh, double& mfsum, int& total_cg_iterations, int& total_newton_iterations);
	model_tet* meshTetes;
	double mfsum;
	int total_cg_iterations;
	int total_newton_iterations;


};
