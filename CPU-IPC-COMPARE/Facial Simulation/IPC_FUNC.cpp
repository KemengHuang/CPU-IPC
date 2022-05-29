#include "IPC_FUNC.h"
#include "IPCdistanceFuncs.h"
#include "fem_parameters.h"
#include "IPCtimeStepFuns.h"
#include <iostream>
#include <fstream>
#include <iomanip>
#include <tbb/parallel_for.h>
#include <tbb/parallel_reduce.h>
#include <tbb/spin_mutex.h>
#include "fem_timer.h"
#include <chrono>
#include "Solver.h"
using namespace FEM;
#define EXPORT_OBJ

#ifdef EXPORT_OBJ
//index count to export
int step_index = 0;
#endif
double time_total = 0;
int total_iter = 0;
void buildConstraintStartIndsWithMM(const vector<int>& activeSet,
    const std::vector<MMCVID>& MMActiveSet,
    std::vector<int>& constraintStartInds)
{
    constraintStartInds.resize(1);
    constraintStartInds[0] = 0;
    constraintStartInds.emplace_back(constraintStartInds.back() + activeSet.size());
    constraintStartInds.emplace_back(constraintStartInds.back() + MMActiveSet.size());
}


void computeXTilta(mesh3D& mesh)
{
    mesh.xTilta.resize(mesh.vertexes.size());
    Vector3d gravityDtSq = Vector3d(0, -9.8, 0) * IPC_dt * IPC_dt;
#ifdef USE_TBB
        tbb::parallel_for(0, (int)mesh.vertexes.size(), 1, [&](int vI)
#else
        for (int vI = 0; vI < mesh.vertexes.size(); vI++)
#endif
        {         
            mesh.xTilta[vI] = (mesh.V_prev[vI] + (mesh.velocities[vI] * IPC_dt + gravityDtSq));
        }
#ifdef USE_TBB
        );
#endif
   
}

void computeEGradient(const mesh3D& mesh, vector<Vector3d>& gradient) {
    //Mass part
    Matrix3d massM;
    massM.setIdentity();
#ifdef USE_TBB
    tbb::parallel_for(0, mesh.vertexNum, 1, [&](int vI)
#else
    for (int vI = 0; vI < mesh.vertexNum; vI++)
#endif
    {
        gradient[vI] += (mesh.masses[vI] * massM * (mesh.vertexes[vI] - mesh.xTilta[vI]));
    }
#ifdef USE_TBB
    );
#endif

    //internal part
    //double tolerance = 1e-6;
    double fsum = 0;
#ifdef USE_TBB
    vector<tbb::spin_mutex> countMutex(mesh.vertexNum);
    tbb::parallel_for(0, mesh.tetrahedraNum, 1, [&](int ii)
#else
    for (int ii = 0; ii < mesh.tetrahedraNum; ii++)
#endif
    {
        MatrixXd PFPX = computePFPX3D_double(mesh.DM_triangle_inverse[ii]);
        MatrixXd F = calculateDms3D_double(mesh.vertexes, mesh.tetrahedras[ii], 0) * mesh.DM_triangle_inverse[ii];

        Matrix3d PEPF = computePEPF_StableNHK3D_2_double(F, lengthRate, volumeRate);

        MatrixXd pepf = vec_double(PEPF);
        MatrixXd f = mesh.volum[ii] * PFPX.transpose() * pepf;



        for (int i = 0; i < 12; i++) {
#ifdef USE_TBB
            countMutex[mesh.tetrahedras[ii][i / 3]].lock();
            gradient[mesh.tetrahedras[ii][i / 3]][i % 3] += IPC_dt * IPC_dt * f(i, 0);
            countMutex[mesh.tetrahedras[ii][i / 3]].unlock();
#else
            gradient[mesh.tetrahedras[ii][i / 3]][i % 3] += IPC_dt * IPC_dt * f(i, 0);
#endif
        }
    }
#ifdef USE_TBB
    );
#endif
}

void computeGradientAndHessian(mesh3D& mesh, vector<Vector3d>& gradient, BHessian& BH, const Ground& grd) {
    //Mass part
    Matrix3d massM;
    massM.setIdentity();
#ifdef USE_TBB
    tbb::parallel_for(0, mesh.vertexNum, 1, [&](int vI)
#else
    for (int vI = 0; vI < mesh.vertexNum; vI++)
#endif
    {
        gradient[vI] = (mesh.masses[vI] * massM * (mesh.vertexes[vI] - mesh.xTilta[vI]));
    }
#ifdef USE_TBB
    );
#endif

    BH.H12x12.resize(mesh.tetrahedraNum);
#ifdef USE_TBB
    vector<tbb::spin_mutex> countMutex(mesh.vertexNum);
    tbb::parallel_for(0, mesh.tetrahedraNum, 1, [&](int ii)
#else
    for (int ii = 0; ii < mesh.tetrahedraNum; ii++)
#endif
    {
        MatrixXd PFPX = computePFPX3D_double(mesh.DM_triangle_inverse[ii]);
        MatrixXd F = calculateDms3D_double(mesh.vertexes, mesh.tetrahedras[ii], 0) * mesh.DM_triangle_inverse[ii];

        Matrix3d PEPF = computePEPF_StableNHK3D_2_double(F, lengthRate, volumeRate);

        MatrixXd pepf = vec_double(PEPF);
        MatrixXd f = mesh.volum[ii] * PFPX.transpose() * pepf;

        for (int i = 0; i < 4; i++) {
#ifdef USE_TBB
            countMutex[mesh.tetrahedras[ii][i]].lock();
            gradient[mesh.tetrahedras[ii][i]] += IPC_dt * IPC_dt * f.block<3, 1>(3 * i, 0);//f.template segment<3>(i * 3);//f(i, 0);
            countMutex[mesh.tetrahedras[ii][i]].unlock();
#else
            gradient[mesh.tetrahedras[ii][i / 3]][i % 3] += IPC_dt * IPC_dt * f(i, 0);
#endif
        }

        //std::cout << F << std::endl;
        MatrixXd Hq = project_StabbleNHK_2_H_3D(F, lengthRate, volumeRate);

        MatrixXd HE = mesh.volum[ii] * IPC_dt * IPC_dt * PFPX.transpose() * Hq * PFPX;
        BH.H12x12[ii] = (HE);

    }
#ifdef USE_TBB
    );
#endif
    BH.D4Index = mesh.tetrahedras;

    VectorXd constraintVals;
    int offset = 0;
    Evaluate_GroundConstraintVals(grd, mesh, constraintVals, offset);
    Matrix3d nnT = grd.normal * grd.normal.transpose();


    tbb::spin_mutex groundMutex;//, countMutex3;
    tbb::parallel_for(0, (int)constraintVals.size(), 1, [&](int cI)
        {
            double g_b, H_b;
            compute_g_b(constraintVals[cI], mesh.Hhat, g_b);
            compute_H_b(constraintVals[cI], mesh.Hhat, H_b);
            int vI = mesh.Environment_ActiveSet[cI];
            double dist = grd.normal.dot(mesh.vertexes[vI]) - grd.D;
            gradient[vI] += mesh.Kappa * g_b * 2.0 * dist * grd.normal;
            double param = 4.0 * H_b * constraintVals[cI] + 2.0 * g_b;
            if (param > 0) {
                Matrix3d Hpg = mesh.Kappa * param * nnT;
                groundMutex.lock();
                BH.H3x3.push_back(Hpg);
                BH.D1Index.push_back(vI);
                groundMutex.unlock();
            }
        }
    );

    offset = constraintVals.size();
    //cout << "self c num:";
    cout << constraintVals.size() << endl;
    Evaluate_SelfPTConstraintVals(mesh, constraintVals, offset);

    tbb::parallel_for(offset, (int)constraintVals.size(), 1, [&](int cI)
        //for (int cI = offset; cI < constraintVals.size(); ++cI) 
        {
            compute_g_b(constraintVals[cI], mesh.Hhat, constraintVals[cI]);
        }
    );
    //cout << "mesh.Self_ActiveSet.size()" << endl;
    //cout << mesh.Self_ActiveSet.size() << endl;
    compute_g_dpt(mesh, mesh.Self_ActiveSet, constraintVals, gradient, offset, mesh.Kappa);
    //compute_g_dpt_new(mesh, mesh.Self_ActiveSet, constraintVals, gradient, offset, mesh.Kappa, mesh.Hhat);

    compute_g_dee(mesh, gradient, mesh.Hhat, mesh.Kappa);

    //double totalTime = 0;
    //if (mesh.Self_ActiveSet.size() > 0) {
    //    vector<MMCVID> test_ActiveSet;
    //    for (int i = 0; i < 1000000; i++) {
    //        test_ActiveSet.emplace_back(-6, 7, 8, -1);
    //    }
    //    mesh.Self_ActiveSet = test_ActiveSet;

    //    for (int i = 0; i < 1; i++) {

    //        HighResolutionTimerForWin timer;
    //        timer.set_start();
    //        compute_H_dpt(mesh, BH, mesh.Hhat, mesh.Kappa);
    //        timer.set_end();
    //        float time0 = timer.get_millisecond();
    //        totalTime += time0;
    //    }
    //    cout << "projection time consuming:  " << totalTime / 1 << endl;
    //    system("pause");
    //}
    //


    compute_H_dpt(mesh, BH, mesh.Hhat, mesh.Kappa);
    //compute_H_dpt_new(mesh, BH, mesh.Hhat, mesh.Kappa, mesh.Hhat);
    compute_H_dee(mesh, BH, mesh.Hhat, mesh.Kappa);
}

void PCG_Precondition(const mesh3D& mesh, const BHessian& BH, const vector<Vector3d>& gradient, vector<Vector3d>& P) {
    for (int i = 0;i < BH.D4Index.size();i++) {
        for (int j = 0;j < 12;j++) {
            P[BH.D4Index[i][j / 3]][j % 3] += BH.H12x12[i](j, j);
        }
    }
    for (int i = 0;i < BH.D3Index.size();i++) {
        for (int j = 0;j < 9;j++) {
            P[BH.D3Index[i][j / 3]][j % 3] += BH.H9x9[i](j, j);
        }
    }
    for (int i = 0;i < BH.D2Index.size();i++) {
        for (int j = 0;j < 6;j++) {
            P[BH.D2Index[i][j / 3]][j % 3] += BH.H6x6[i](j, j);
        }
    }
    for (int i = 0;i < BH.D1Index.size();i++) {
        for (int j = 0;j < 3;j++) {
            P[BH.D1Index[i]][j] += BH.H3x3[i](j, j);
        }
    }
    for (int i = 0;i < mesh.vertexNum;i++) {
        P[i] += Vector3d(mesh.masses[i], mesh.masses[i], mesh.masses[i]);
        //P[i] = Vector3d(1, 1, 1);
    }
}


vector<Vector3d> cholmod_solver(BHessian& BH, std::vector<Eigen::Vector3d> gradient, const mesh3D& mesh) {

    vector<Vector3d> direction(mesh.vertexNum, Vector3d(0, 0, 0));
#ifdef NDEBUG
    CholmodSolver solver;
    std::vector<Triplet<double>> triplets = BH.toTriplets(mesh.boundaryTypes);
    int offset = triplets.size();
    //std::chrono::steady_clock::time_point begin = std::chrono::steady_clock::now();
    triplets.resize(offset + 3 * mesh.vertexNum);

    tbb::parallel_for(0, 3 * mesh.vertexNum, 1, [&](int i)
        {
            triplets[offset + i] = Triplet<double>(i, i, mesh.masses[i / 3]);
        }
    );

    Eigen::SparseMatrix<double> sparseFEM(3 * mesh.vertexNum, 3 * mesh.vertexNum);
    sparseFEM.setFromTriplets(triplets.begin(), triplets.end());
    auto vec_gradient = Eigen::VectorXd(gradient.size() * 3);

    tbb::parallel_for(0, (int)(gradient.size()), 1, [&](int i)
        {
            if (mesh.boundaryTypes[i] == 0)
                vec_gradient.block<3, 1>(3 * i, 0) = gradient[i];
            else
                vec_gradient.block<3, 1>(3 * i, 0) = Vector3d(0, 0, 0);
        }
    );
    auto result = Eigen::VectorXd(gradient.size() * 3);
    
    solver.set_pattern(sparseFEM);
    solver.solve(vec_gradient, result);
    tbb::parallel_for(0, (int)(gradient.size()), 1, [&](int i)
        {
            direction[i] = result.block<3, 1>(3 * i, 0);
        }
    );

    //std::chrono::steady_clock::time_point end = std::chrono::steady_clock::now();
    //std::cout << "Time for computing cg only = "<< std::chrono::duration_cast<std::chrono::milliseconds>(end - begin).count() << "[ms]" << std::endl;

#endif
    return direction;
}

vector<Vector3d> PCG_Solver(const mesh3D& mesh, const BHessian& BH, const vector<Vector3d>& gradient) {
    std::vector<Vector3d> P(mesh.vertexNum, Vector3d(0, 0, 0));
    PCG_Precondition(mesh, BH, gradient, P);

    vector<Vector3d> tempDeltaX(mesh.vertexNum, Vector3d(0, 0, 0));

    double deltaN = 0;
    double localOptimal = DBL_MAX;
    bool getLocalOpt = false;
    vector<Vector3d> dX(mesh.vertexNum, Vector3d(0, 0, 0));

    double delta0 = 0;
    double deltaO = 0;

    vector<Vector3d> r(mesh.vertexNum, Vector3d(0, 0, 0));
    vector<Vector3d> c(mesh.vertexNum, Vector3d(0, 0, 0));

#ifdef USE_TBB
    vector<tbb::spin_mutex> countMutex(mesh.vertexNum);
    delta0 = parallel_reduce(
        tbb::blocked_range<int>(0, mesh.vertexNum), 0.0, [&](const tbb::blocked_range<int>& r, double temp_delta) {
            for (int i = r.begin(); i != r.end(); i++) {
                double vx = 1 / (P[i][0]);// abs(P(i * 3)) > 1e-5 ? 1 / (P(i * 3)) : 1;
                double vy = 1 / (P[i][1]);// abs(P(i * 3) + 1) > 1e-5 ? 1 / (P(i * 3) + 1) : 1;
                double vz = 1 / (P[i][2]);// abs(P(i * 3) + 2) > 1e-5 ? 1 / (P(i * 3) + 2) : 1;
                Vector3d filter_b = mesh.Constraints[i] * gradient[i];
                temp_delta += filter_b[0] * filter_b[0] * vx;
                temp_delta += filter_b[1] * filter_b[1] * vy;
                temp_delta += filter_b[2] * filter_b[2] * vz;
            }
            return temp_delta;
        },
        [&](double left, double right) {
            return left + right;
        }
        );
#else
    for (int i = 0; i < mesh.vertexNum; i++) {
        double vx = 1 / (P[i][0]);// abs(P(i * 3)) > 1e-5 ? 1 / (P(i * 3)) : 1;
        double vy = 1 / (P[i][1]);// abs(P(i * 3) + 1) > 1e-5 ? 1 / (P(i * 3) + 1) : 1;
        double vz = 1 / (P[i][2]);// abs(P(i * 3) + 2) > 1e-5 ? 1 / (P(i * 3) + 2) : 1;
        Vector3d filter_b = mesh.Constraints[i] * gradient[i];
        delta0 += filter_b[0] * filter_b[0] * vx;
        delta0 += filter_b[1] * filter_b[1] * vy;
        delta0 += filter_b[2] * filter_b[2] * vz;
    }
#endif

#ifdef USE_TBB
    deltaN = parallel_reduce(
        tbb::blocked_range<int>(0, mesh.vertexNum), 0.0, [&](const tbb::blocked_range<int>& rg, double temp_deltaN) {
            for (int i = rg.begin(); i != rg.end(); i++) {
                r[i] = gradient[i];
                r[i] = mesh.Constraints[i] * r[i];
                c[i][0] = r[i][0] / P[i][0];
                c[i][1] = r[i][1] / P[i][1];
                c[i][2] = r[i][2] / P[i][2];
                c[i] = mesh.Constraints[i] * c[i];

                temp_deltaN += c[i][0] * r[i][0];
                temp_deltaN += c[i][1] * r[i][1];
                temp_deltaN += c[i][2] * r[i][2];
            }
            return temp_deltaN;
        },
        [&](double left, double right) {
            return left + right;
        }
        );
#else
    for (int i = 0; i < mesh.vertexNum; i++)
    {
        r[i] = gradient[i];
        r[i] = mesh.Constraints[i] * r[i];
        c[i][0] = r[i][0] * P[i][0];
        c[i][1] = r[i][1] * P[i][1];
        c[i][2] = r[i][2] * P[i][2];
        c[i] = mesh.Constraints[i] * c[i];

        deltaN += c[i][0] * r[i][0];
        deltaN += c[i][1] * r[i][1];
        deltaN += c[i][2] * r[i][2];
    }
#endif

    double errorRate = 1e-6;
    //std::cout << "cpu  delta0:   " << delta0 << "      deltaN:   " << deltaN << endl;
    //system("pause");
    //PCG main loop
    int cgCounts = 0;
    while (cgCounts<30000 && deltaN > errorRate * delta0) {
        cgCounts++;
        //std::cout << "delta0:   " << delta0 << "      deltaN:   " << deltaN  << endl;
        vector<Vector3d> q(mesh.vertexNum, Vector3d(0, 0, 0));
#ifdef USE_TBB
        tbb::parallel_for(0, (int)BH.D4Index.size(), 1, [&](int ii)
#else
        for (int ii = 0; ii < BH.D4Index.size(); ii++)
#endif
        {
            MatrixXd H = BH.H12x12[ii];
            VectorXd tempC(12);
            for (int i = 0; i < 12; i++) {
                tempC(i) = c[BH.D4Index[ii][i / 3]][i % 3];
            }
            VectorXd tempQ = H * tempC;
            for (int i = 0; i < 12; i++) {
#ifdef USE_TBB
                countMutex[BH.D4Index[ii][i / 3]].lock();
                q[BH.D4Index[ii][i / 3]][i % 3] += tempQ(i);
                countMutex[BH.D4Index[ii][i / 3]].unlock();
#else
                q[BH.D4Index[ii][i / 3]][i % 3] += tempQ(i);
#endif
            }
        }
#ifdef USE_TBB
        );
#endif

#ifdef USE_TBB
        tbb::parallel_for(0, (int)BH.D3Index.size(), 1, [&](int ii)
#else
        for (int ii = 0; ii < BH.D3Index.size(); ii++)
#endif
        {
            MatrixXd H = BH.H9x9[ii];

            VectorXd tempC(9);

            for (int i = 0; i < 9; i++) {
                tempC(i) = c[BH.D3Index[ii][i / 3]][i % 3];
            }
            VectorXd tempQ = H * tempC;

            for (int i = 0; i < 9; i++) {
#ifdef USE_TBB
                countMutex[BH.D3Index[ii][i / 3]].lock();
                q[BH.D3Index[ii][i / 3]][i % 3] += tempQ(i);
                countMutex[BH.D3Index[ii][i / 3]].unlock();
#else
                q[BH.D3Index[ii][i / 3]][i % 3] += tempQ(i);
#endif
            }
        }
#ifdef USE_TBB
        );
#endif

#ifdef USE_TBB
        tbb::parallel_for(0, (int)BH.D2Index.size(), 1, [&](int ii)
#else
        for (int ii = 0; ii < BH.D2Index.size(); ii++)
#endif
        {
            MatrixXd H = BH.H6x6[ii];

            VectorXd tempC(6);

            for (int i = 0; i < 6; i++) {
                tempC(i) = c[BH.D2Index[ii][i / 3]][i % 3];
            }
            VectorXd tempQ = H * tempC;

            for (int i = 0; i < 6; i++) {
#ifdef USE_TBB
                countMutex[BH.D2Index[ii][i / 3]].lock();
                q[BH.D2Index[ii][i / 3]][i % 3] += tempQ(i);
                countMutex[BH.D2Index[ii][i / 3]].unlock();
#else
                q[BH.D2Index[ii][i / 3]][i % 3] += tempQ(i);
#endif
            }
        }
#ifdef USE_TBB
        );
#endif

#ifdef USE_TBB
        tbb::parallel_for(0, mesh.vertexNum, 1, [&](int ii)
#else
        for (int ii = 0; ii < mesh.vertexNum; ii++)
#endif
        {
            q[ii] += mesh.masses[ii] * c[ii];
        }
#ifdef USE_TBB
        );
#endif


#ifdef USE_TBB
        tbb::parallel_for(0, (int)BH.D1Index.size(), 1, [&](int ii)
#else
        for (int ii = 0; ii < BH.D1Index.size(); ii++)
#endif
        {
            MatrixXd H = BH.H3x3[ii];

            VectorXd tempC(3);

            for (int i = 0; i < 3; i++) {
                tempC(i) = c[BH.D1Index[ii]][i];
            }
            VectorXd tempQ = H * tempC;

            for (int i = 0; i < 3; i++) {
#ifdef USE_TBB
                countMutex[BH.D1Index[ii]].lock();
                q[BH.D1Index[ii]][i] += tempQ(i);
                countMutex[BH.D1Index[ii]].unlock();
#else
                q[BH.D1Index[ii]][i] += tempQ(i);
#endif
            }
        }
#ifdef USE_TBB
        );
#endif

        double tempSum = 0;

#ifdef USE_TBB
        tempSum = parallel_reduce(
            tbb::blocked_range<int>(0, mesh.vertexNum), 0.0, [&](const tbb::blocked_range<int>& rg, double temp_sum) {
                for (int i = rg.begin(); i != rg.end(); i++) {
                    q[i] = mesh.Constraints[i] * q[i];
                    temp_sum += (c[i][0] * q[i][0] + c[i][1] * q[i][1] + c[i][2] * q[i][2]);
                }
                return temp_sum;
            },
            [&](double left, double right) {
                return left + right;
            }
            );
#else
        for (int i = 0; i < mesh.vertexNum; i++) {
            q[i] = mesh.Constraints[i] * q[i];
            tempSum += (c[i][0] * q[i][0] + c[i][1] * q[i][1] + c[i][2] * q[i][2]);
        }
#endif

        double alpha = deltaN / tempSum;
        //cout << "tempSum:------------------"<<tempSum << endl;
        //if(tempSum)
        deltaO = deltaN;
        deltaN = 0;
        vector<Vector3d> s(mesh.vertexNum, Vector3d(0, 0, 0));

#ifdef USE_TBB
        deltaN = parallel_reduce(
            tbb::blocked_range<int>(0, mesh.vertexNum), 0.0, [&](const tbb::blocked_range<int>& rg, double temp_deltaN) {
                for (int i = rg.begin(); i != rg.end(); i++) {
                    dX[i][0] = dX[i][0] + alpha * c[i][0];
                    dX[i][1] = dX[i][1] + alpha * c[i][1];
                    dX[i][2] = dX[i][2] + alpha * c[i][2];

                    r[i][0] = r[i][0] - alpha * q[i][0];
                    r[i][1] = r[i][1] - alpha * q[i][1];
                    r[i][2] = r[i][2] - alpha * q[i][2];

                    s[i][0] = r[i][0] / P[i][0];
                    s[i][1] = r[i][1] / P[i][1];
                    s[i][2] = r[i][2] / P[i][2];

                    temp_deltaN += (r[i][0] * s[i][0] + r[i][1] * s[i][1] + r[i][2] * s[i][2]);
                }
                return temp_deltaN;
            },
            [&](double left, double right) {
                return left + right;
            }
            );
#else
        for (int i = 0; i < mesh.vertexNum; i++) {
            dX[i][0] = dX[i][0] + alpha * c[i][0];
            dX[i][1] = dX[i][1] + alpha * c[i][1];
            dX[i][2] = dX[i][2] + alpha * c[i][2];

            r[i][0] = r[i][0] - alpha * q[i][0];
            r[i][1] = r[i][1] - alpha * q[i][1];
            r[i][2] = r[i][2] - alpha * q[i][2];

            s[i][0] = r[i][0] * P[i][0];
            s[i][1] = r[i][1] * P[i][1];
            s[i][2] = r[i][2] * P[i][2];

            deltaN += (r[i][0] * s[i][0] + r[i][1] * s[i][1] + r[i][2] * s[i][2]);
        }
#endif
        
        if (deltaN < localOptimal) {
            localOptimal = deltaN;
            getLocalOpt = true;
#ifdef USE_TBB
            tbb::parallel_for(0, mesh.vertexNum, 1, [&](int j)
#else
            for (int j = 0; j < mesh.vertexNum; j++)
#endif
            {
                tempDeltaX[j] = Vector3d(dX[j]);
            }
#ifdef USE_TBB
            );
#endif
        }

#ifdef USE_TBB
        tbb::parallel_for(0, mesh.vertexNum, 1, [&](int i)
#else
        for (int i = 0; i < mesh.vertexNum; i++)
#endif
        {
            c[i][0] = s[i][0] + (deltaN / deltaO) * c[i][0];
            c[i][1] = s[i][1] + (deltaN / deltaO) * c[i][1];
            c[i][2] = s[i][2] + (deltaN / deltaO) * c[i][2];
            c[i] = mesh.Constraints[i] * c[i];
        }
#ifdef USE_TBB
        );
#endif
    }
    printf("cg counts: %d\n", cgCounts);
    if (localOptimal)
        return tempDeltaX;
    return dX;
}

void calculateMovingDirection(const mesh3D& mesh, BHessian& BH, vector<Vector3d> gradient, vector<Vector3d>& direction) {
#ifdef NDEBUG
    direction = cholmod_solver(BH, gradient, mesh);
    //direction = PCG_Solver(mesh, BH, gradient);
#else
    direction = PCG_Solver(mesh, BH, gradient);
#endif
}

void computeEnergyVal(const mesh3D& mesh, double& energyVal, const Ground& gd, double Kappa)
{
    energyVal = 0;
    energyVal = getObjEnergy_StableNHK2_3D(mesh.vertexes, mesh, lengthRate, volumeRate);
    //energyVal += getObjEnergy_AniostroI5_3D(mesh.vertexes, mesh, lengthRate, contract_ratio);
    energyVal *= IPC_dt * IPC_dt;
    double deltaE = 0;
#ifdef USE_TBB
    deltaE = parallel_reduce(
        tbb::blocked_range<int>(0, mesh.vertexNum), 0.0, [&](const tbb::blocked_range<int>& rg, double temp_deltaE) {
            for (int i = rg.begin(); i != rg.end(); i++) {
                temp_deltaE += ((mesh.vertexes[i] - mesh.xTilta[i]).squaredNorm() * mesh.masses[i] / 2.0);
            }
            return temp_deltaE;
        },
        [&](double left, double right) {
            return left + right;
        }
        );
#else
    //double tepsum = 0;
    for (int i = 0; i < mesh.vertexNum; i++) {
        double tep = ((mesh.vertexes[i] - mesh.xTilta[i]).squaredNorm() * mesh.masses[i] / 2.0);
        deltaE += tep;
        //energyVal += ((mesh.vertexes[i] - mesh.xTilta[i]).squaredNorm() * mesh.masses[i] / 2.0);
    }
#endif
    energyVal += deltaE;
    Eigen::VectorXd constraintVals, bVals;
    int startCI = constraintVals.size();
    Evaluate_GroundConstraintVals(gd, mesh, constraintVals, startCI);
    bVals.conservativeResize(constraintVals.size());
    

    tbb::parallel_for(startCI, (int)constraintVals.size(), 1, [&](int cI)
        //for (int cI = startCI; cI < constraintVals.size(); ++cI) 
        {

            compute_b(constraintVals[cI], mesh.Hhat, bVals[cI]);

        }
    );


    startCI = constraintVals.size();
    Evaluate_SelfPTConstraintVals(mesh, constraintVals, startCI);
    bVals.conservativeResize(constraintVals.size());
    //TODO: parallelize
    tbb::parallel_for(startCI, (int)constraintVals.size(), 1, [&](int cI)
        //for (int cI = startCI; cI < constraintVals.size(); ++cI) 
        {

            compute_b(constraintVals[cI], mesh.Hhat, bVals[cI]);
            int duplication = mesh.Self_ActiveSet[cI - startCI][3];
            if (duplication < -1) {
                // PP or PE, handle duplication
                bVals[cI] *= -duplication;
            }
        }

    );

    startCI = constraintVals.size();
    Evaluate_SelfEEConstraintVals(mesh, constraintVals, startCI);
    //SelfCollisionHandler<dim>::evaluateConstraints(data, paraEEMMCVIDSet.back(), constraintVals);
    bVals.conservativeResize(constraintVals.size());

    tbb::parallel_for(startCI, (int)constraintVals.size(), 1, [&](int cI)
        //for (int cI = startCI; cI < constraintVals.size(); ++cI) 
        {

            const MMCVID& MMCVIDI = mesh.Self_EE_ActiveSet[cI - startCI];
            double eps_x, e;
            if (MMCVIDI[3] >= 0) {
                // EE
                compute_eps_x(mesh, MMCVIDI[0], MMCVIDI[1], MMCVIDI[2], MMCVIDI[3], eps_x);
                compute_e(mesh.vertexes[MMCVIDI[0]], mesh.vertexes[MMCVIDI[1]], mesh.vertexes[MMCVIDI[2]], mesh.vertexes[MMCVIDI[3]], eps_x, e);
            }
            else {
                // PP or PE
                const std::pair<int, int>& eIeJ = mesh.Self_EEeIe_ActiveSet[cI - startCI];
                const std::pair<int, int>& eI = mesh.surfEdges[eIeJ.first];
                const std::pair<int, int>& eJ = mesh.surfEdges[eIeJ.second];
                compute_eps_x(mesh, eI.first, eI.second, eJ.first, eJ.second, eps_x);
                compute_e(mesh.vertexes[eI.first], mesh.vertexes[eI.second], mesh.vertexes[eJ.first], mesh.vertexes[eJ.second], eps_x, e);
            }
            compute_b(constraintVals[cI], mesh.Hhat, bVals[cI]);
            bVals[cI] *= e;
        }

    );
    //double bE
    energyVal += Kappa * bVals.sum();
}

void stepForward(const vector<Vector3d>& dataV0,
    const vector<Vector3d>& searchDir,
    mesh3D& mesh,
    double stepSize, bool boundary_update = false)
{
    //assert(dataV0.rows() == data.V.rows());
    //assert(data.V.rows() * dim == searchDir.size());
    //assert(data.V.rows() == result.V.rows());

#ifdef USE_TBB
    tbb::parallel_for(0, mesh.vertexNum, 1, [&](int vI)
#else
    for (int vI = 0; vI < mesh.vertexNum; vI++)
#endif
    {
        if(mesh.boundaryTypes[vI]==0 || boundary_update)
            mesh.vertexes[vI] = dataV0[vI] - stepSize * searchDir[vI];
    }
#ifdef USE_TBB
    );
#endif
}





bool checkEdgeTriIntersectionIfAny(const mesh3D& mesh,
    SpatialHash& sh)
{
    Eigen::ArrayXi intersected(mesh.surface.size());
    intersected.setZero();
#ifdef USE_TBB
    tbb::parallel_for(0, (int)mesh.surface.size(), 1, [&](int sfI)
#else
    for (int sfI = 0; sfI < mesh.surface.size(); ++sfI)
#endif
    {
        const Eigen::RowVector4i& sfVInd = mesh.surface[sfI].transpose();
        int coDim_sfI = 3;
#ifndef NUSE_SH_CCS
        std::unordered_set<int> sEdgeInds;
        sh.queryTriangleForEdges(mesh.vertexes[sfVInd[0]], mesh.vertexes[sfVInd[1]], mesh.vertexes[sfVInd[2]], 0.0, sEdgeInds);
        for (const auto& eI : sEdgeInds)
#else
        for (int eI = 0; eI < mesh.surfEdges.size(); ++eI)
#endif
        {
            const auto& meshEI = mesh.surfEdges[eI];
            if (meshEI.first == sfVInd[0] || meshEI.first == sfVInd[1] || meshEI.first == sfVInd[2] || meshEI.second == sfVInd[0] || meshEI.second == sfVInd[1] || meshEI.second == sfVInd[2]) {
                continue;
            }

            int coDim_eI = 3;//mesh.vICoDim(meshEI.first);


            if (IglUtils::segTriIntersect(mesh.vertexes[meshEI.first], mesh.vertexes[meshEI.second],
                mesh.vertexes[sfVInd[0]], mesh.vertexes[sfVInd[1]], mesh.vertexes[sfVInd[2]])) {
                intersected[sfI] = 1;
                break;
            }
        }
    }
#ifdef USE_TBB
    );
#endif

    if ((intersected != 0).any()) {
        return false;
    }

    return true;
}


bool isIntersected(const Ground& grd,
    SpatialHash& sh,
    const mesh3D& mesh,
    const vector<Vector3d>& V0)
{
    Eigen::VectorXd constraint_vals;
    Evaluate_GroundConstraintVals(grd, mesh, constraint_vals, 0);
    for (int vI = 0; vI < constraint_vals.size(); ++vI) {
        if (constraint_vals[vI] <= 0.0) {
            return true;
        }
    }


    if (!checkEdgeTriIntersectionIfAny(mesh, sh)) {

        return true;
    }

    return false;
}

void buildCollisionSets(mesh3D& mesh,
    SpatialHash& sh,
    const Ground& gd,
    bool rehash = true) {
    if (rehash) {
        sh.build(mesh, mesh.averageEdgeLenth);
    }
    gd.calculateActivateSet(mesh);
    sh.calculateActivateSet(mesh);
}

bool lineSearch(mesh3D& mesh,
    SpatialHash& sh,
    const Ground& gd,
    const vector<Vector3d>& searchDir,
    const vector<Vector3d>& gradient,
    double& stepSize,
    double armijoParam,
    double lowerBound,
    double Kappa)
{
    //std::stringstream msg;



    bool stopped = false;
    double lastEnergyVal;

    computeEnergyVal(mesh, lastEnergyVal, gd, Kappa);

    //msg << "E_last = " << lastEnergyVal << " stepSize: " << stepSize << " -> ";

    // const double m = searchDir.dot(gradient);
    // const double c1m = 1.0e-4 * m;
    double c1m = 0.0;
    armijoParam = 1e-4;
    if (armijoParam > 0.0) {
#ifdef USE_TBB
        c1m = parallel_reduce(
            tbb::blocked_range<int>(0, mesh.vertexNum), 0.0, [&](const tbb::blocked_range<int>& rg, double temp_deltaE) {
                for (int i = rg.begin(); i != rg.end(); i++) {
                    temp_deltaE += gradient[i].dot(searchDir[i]);
                }
                return temp_deltaE;
            },
            [&](double left, double right) {
                return left + right;
            }
            );
#else
        //double test = 0;
        for (int i = 0; i < mesh.vertexNum; i++) {
            //test = searchDir[i].dot(searchDir[i]);
            c1m += searchDir[i].dot(searchDir[i]);
        }
#endif
    }

    vector<Vector3d> resultV0 = mesh.vertexes;

    stepForward(resultV0, searchDir, mesh, stepSize);

    bool rehash = true;
    int numOfIntersect = 0;
    sh.build(mesh, mesh.averageEdgeLenth);
    while (isIntersected(gd, sh, mesh, resultV0)) {
        numOfIntersect++;
        printf("intersect\n");
        stepSize /= 2.0;
        stepForward(resultV0, searchDir, mesh, stepSize);
        sh.build(mesh, mesh.averageEdgeLenth);
    }
    //msg << stepSize << "(safeGuard_IP) -> ";
    rehash = false;

    double testingE;

    buildCollisionSets(mesh, sh, gd, false);
    computeEnergyVal(mesh, testingE, gd, Kappa);


    //cout << "test test test:  " << abs(lastEnergyVal - mesh.restSNKE)<<"        "<< c1m<< endl;

    int numOfLineSearch = 0;
    double LFStepSize = stepSize;
    while ((testingE > lastEnergyVal + stepSize * c1m) && (stepSize > 1e-3 * LFStepSize) && abs(testingE - lastEnergyVal) / abs(lastEnergyVal - mesh.restSNKE) > 1e-8 / IPC_dt / (1 << (numOfLineSearch + 1))) {
        // fprintf(out, "%.9le %.9le\n", stepSize, testingE);
        //if (stepSize == 1.0) {
        //    stepSize /= 2.0;
        //}
        //else {
        printf("testingE:  %f       lastEnergyVal:   %f\n", testingE, lastEnergyVal);
        stepSize /= 2.0;
        //}

        ++numOfLineSearch;
        if (stepSize == 0.0) {
            stopped = true;
            break;
        }

        stepForward(resultV0, searchDir, mesh, stepSize);

        //msg << stepSize << " -> ";

        buildCollisionSets(mesh, sh, gd);
        computeEnergyVal(mesh, testingE, gd, Kappa);
    }

    if (stepSize < LFStepSize) {
        bool needRecomputeCS = false;
        while (isIntersected(gd, sh, mesh, resultV0)) {
            stepSize /= 2.0;
            numOfIntersect++;
            printf("intersect\n");
            stepForward(resultV0, searchDir, mesh, stepSize);
            sh.build(mesh, mesh.averageEdgeLenth);
            needRecomputeCS = true;
        }
        if (needRecomputeCS) {
            buildCollisionSets(mesh, sh, gd, false);
        }
    }

    //msg << stepSize << "(armijo) ";
    printf("lineSearch step: %f\n", stepSize);
    lastEnergyVal = testingE;
    //if (stepSize >= 1.0 / (1 << (1 + numOfLineSearch)) && !numOfIntersect /*&& !numOfLineSearch*/ && abs(testingE - lastEnergyVal) / abs(lastEnergyVal - mesh.restSNKE) < 1e-8 / IPC_dt / (1 << (numOfLineSearch + 1))/* / vertexNum*/ /*&& (testingE > lastEnergyVal + c1m * alpha)*/) {
    //    stopped = true;
    //}
    return stopped;
}

void suggestKappa(double& kappa, const double& Hhat, const double& bboxDiagSize2, const double& meanMass)
{
    double H_b;
    //double bboxDiagSize2 = (maxConer - minConer).squaredNorm();
    compute_H_b(1.0e-16 * bboxDiagSize2, Hhat, H_b);
    kappa = 1e11 * meanMass / (4.0e-16 * bboxDiagSize2 * H_b);
}

void upperBoundKappa(double& kappa, const double& Hhat, const double& bboxDiagSize2, const double& meanMass)
{
    double H_b;
    //double bboxDiagSize2 = (maxConer - minConer).squaredNorm();
    compute_H_b(1.0e-16 * bboxDiagSize2, Hhat, H_b);
    double kappaMax = 100 * 1e11 * meanMass / (4.0e-16 * bboxDiagSize2 * H_b);
    if (kappa > kappaMax) {
        kappa = kappaMax;
    }
}

void initKappa(mesh3D& mesh, const Ground& grd, double& kappa)
{
    std::vector<int> constraintStartInds;
    buildConstraintStartIndsWithMM(mesh.Environment_ActiveSet, mesh.Self_ActiveSet, constraintStartInds);

    if (constraintStartInds.back()) {
        vector<Vector3d> g_E(mesh.vertexNum, Vector3d(0,0,0)), g_c(mesh.vertexNum, Vector3d(0, 0, 0));
        computeEGradient(mesh, g_E);
        VectorXd constraintVal;
        int startCI = constraintStartInds[0];
        Evaluate_GroundConstraintVals(grd, mesh, constraintVal, startCI);
        //animConfig.collisionObjects[coI]->evaluateConstraints(result, activeSet[coI], constraintVal);
        
        for (int cI = startCI; cI < constraintStartInds[1]; ++cI) {
            compute_g_b(constraintVal[cI], mesh.Hhat, constraintVal[cI]);
            int vI = mesh.Environment_ActiveSet[cI];
            double dist = grd.normal.dot(mesh.vertexes[vI]) - grd.D;
            g_c[vI] += constraintVal[cI] * 2.0 * dist * grd.normal;
        }

        startCI = constraintStartInds[1];
        Evaluate_SelfPTConstraintVals(mesh, constraintVal, startCI);

        for (int cI = startCI; cI < constraintVal.size(); ++cI) {
            compute_g_b(constraintVal[cI], mesh.Hhat, constraintVal[cI]);
        }
        //cout << 2 << endl;
        compute_g_dpt(mesh, mesh.Self_ActiveSet, constraintVal, g_c, startCI, 1);

        double gsum = 0, gsnorm = 0;
        for (int i = 0;i < mesh.vertexNum;i++) {
            gsum += g_c[i].dot(g_E[i]);
            gsnorm += g_c[i].squaredNorm();
        }
        // balance current gradient at constrained DOF
        double minKappa = -gsum / gsnorm;
        if (minKappa > 0.0) {
            kappa = minKappa;
        }
        suggestKappa(mesh.Kappa, mesh.Hhat, mesh.bboxDiagSize2, mesh.meanMass);
        if (kappa < minKappa) {
            kappa = minKappa;
        }
        upperBoundKappa(mesh.Kappa, mesh.Hhat, mesh.bboxDiagSize2, mesh.meanMass);
    }

}


void postLineSearch(mesh3D& mesh, const Ground& grd, double alpha, double& kappa)
{
    if (kappa == 0.0) {
        initKappa(mesh, grd, kappa);
    }
    else {
        //TODO: avoid recomputation of constraint functions
        bool updateKappa = false;
        for (int i = 0; i < mesh.closeConstraintID.size(); ++i) {

            double d = grd.calculateGapFromObj(mesh, mesh.closeConstraintID[i]);
            if (d <= mesh.closeConstraintVal[i]) {
                updateKappa = true;
                break;
            }
        }
        if (!updateKappa) {
            for (int i = 0; i < mesh.closeMConstraintID.size(); ++i) {

                double d = SelfConstraintVal(mesh, mesh.closeMConstraintID[i]);

                if (d <= mesh.closeMConstraintVal[i]) {
                    updateKappa = true;
                    break;
                }
            }
        }
        if (updateKappa) {
            kappa *= 2.0;
            upperBoundKappa(mesh.Kappa, mesh.Hhat, mesh.bboxDiagSize2, mesh.meanMass);
        }

        mesh.closeConstraintID.resize(0);
        mesh.closeMConstraintID.resize(0);
        mesh.closeConstraintVal.resize(0);
        mesh.closeMConstraintVal.resize(0);
        Eigen::VectorXd constraintVal;

        int constraintValIndStart = constraintVal.size();
        Evaluate_GroundConstraintVals(grd, mesh, constraintVal, constraintValIndStart);

        //double bboxDiagSize2 = (mesh.maxConer - mesh.minConer).squaredNorm();
        //double dTol = 1e-18 * mesh.bboxDiagSize2;
        //animConfig.collisionObjects[coI]->evaluateConstraints(result,
        //    activeSet[coI], constraintVal);
        //double bboxDiagSize2 = (mesh.maxConer - mesh.minConer).squaredNorm();
        for (int i = 0; i < mesh.Environment_ActiveSet.size(); ++i) {
            if (constraintVal[constraintValIndStart + i] < mesh.dTol) {
                mesh.closeConstraintID.emplace_back(mesh.Environment_ActiveSet[i]);
                mesh.closeConstraintVal.emplace_back(constraintVal[constraintValIndStart + i]);
            }
        }

        constraintValIndStart = constraintVal.size();
        Evaluate_SelfPTConstraintVals(mesh, constraintVal, constraintValIndStart);

        for (int i = 0; i < mesh.Self_ActiveSet.size(); ++i) {
            // std::cout << MMActiveSet.back()[i][0] << " " << MMActiveSet.back()[i][1] << " " << MMActiveSet.back()[i][2] << " " << MMActiveSet.back()[i][3] << ", d=" << constraintVal[constraintValIndStart + i] << std::endl;
            if (constraintVal[constraintValIndStart + i] < mesh.dTol) {
                mesh.closeMConstraintID.emplace_back(mesh.Self_ActiveSet[i]);
                mesh.closeMConstraintVal.emplace_back(constraintVal[constraintValIndStart + i]);
            }
        }
    }



}

void updateBoundaryMoveDir(mesh3D& mesh, vector<Vector3d>& moveDir, double ipc_dt, double alpha) {
    double angleX = PI / 2.5 * ipc_dt * alpha;
    Matrix3d rotationL, rotationR;
    rotationL << 1, 0, 0, 0, cos(angleX), sin(angleX), 0, -sin(angleX), cos(angleX);
    rotationR << 1, 0, 0, 0, cos(angleX), -sin(angleX), 0, sin(angleX), cos(angleX);

    tbb::parallel_for(0, (int)(moveDir.size()), 1, [&](int i)
        //for (int i = 0; i < moveDir.size();i++) 
        {
            if (mesh.boundaryTypes[i] == 1) {
                moveDir[i] = rotationL * mesh.vertexes[i] - mesh.vertexes[i];
            }
            if (mesh.boundaryTypes[i] == -1) {
                moveDir[i] = rotationR * mesh.vertexes[i] - mesh.vertexes[i];
            }
        }
    );
}


double calculate_distToOpt_PN(const vector<Vector3d>& moveDir) {
    double distToOpt_PN = 0;
#ifdef USE_TBB
    distToOpt_PN = parallel_reduce(
        tbb::blocked_range<int>(0, (int)moveDir.size()), 0.0, [&](const tbb::blocked_range<int>& rg, double temp_max) {
            for (int i = rg.begin(); i != rg.end(); i++) {
                for (int jj = 0; jj < 3; jj++) {
                    if (temp_max < abs(moveDir[i][jj])) {
                        temp_max = abs(moveDir[i][jj]);
                    }
                }
            }
            return temp_max;
        },
        [&](double left, double right) {
            return left > right ? left : right;
        }
        );
#else
    for (int ii = 0; ii < moveDir.size(); ii++) {
        for (int jj = 0; jj < 3; jj++) {
            if (distToOpt_PN < abs(moveDir[ii][jj])) {
                distToOpt_PN = abs(moveDir[ii][jj]);
            }
        }
    }
#endif
    return distToOpt_PN;
}

int solve_subIP(mesh3D& mesh, SpatialHash& sh, Ground& gd, double Kappa) {
    int iterCap = 10000, k = 0;

    vector<Vector3d> moveDir(mesh.vertexNum, Vector3d(0, 0, 0));
    if (false)
    {
        double new_alpha = 1;
        updateBoundaryMoveDir(mesh, moveDir, IPC_dt, new_alpha);

        sh.build(mesh, moveDir, new_alpha, mesh.averageEdgeLenth);
        Self_largestFeasibleStepSize_CCD(mesh, sh, moveDir, 0.8, new_alpha);
        new_alpha *= 0.5;

        //updateBoundary(TetMesh, alpha);
        updateBoundaryMoveDir(mesh, moveDir, IPC_dt, new_alpha);
        vector<Vector3d> resultV0 = mesh.vertexes;
        stepForward(resultV0, moveDir, mesh, new_alpha, true);

        bool rehash = true;


        sh.build(mesh, mesh.averageEdgeLenth);
        int numOfIntersect = 0;

        while (isIntersected(gd, sh, mesh, mesh.vertexes)) {
            new_alpha /= 2.0;
            updateBoundaryMoveDir(mesh, moveDir, IPC_dt, new_alpha);
            stepForward(resultV0, moveDir, mesh, new_alpha, true);
            sh.build(mesh, mesh.averageEdgeLenth);
        }
        cout << "new_alpha:";
        cout << new_alpha << endl;
        buildCollisionSets(mesh, sh, gd, false);
    }
    double totalTimeStep = 0;
    for (; k < iterCap; ++k) {

        HighResolutionTimerForWin timer0, timer1, timer2, timer3, timer4;

        vector<Vector3d> gradient(mesh.vertexNum, Vector3d(0, 0, 0));
        BHessian BH;
        timer0.set_start();
        computeGradientAndHessian(mesh, gradient, BH, gd);
        timer0.set_end();
        timer1.set_start();
        double gradSqNorm = vector_squareNorm(gradient);
        //vector<Vector3d> moveDir;
        double distToOpt_PN = 0;

        distToOpt_PN = calculate_distToOpt_PN(moveDir);
        cout << "distToOpt_PN" << endl;
        cout << distToOpt_PN << endl;

        bool gradVanish = (distToOpt_PN < sqrt(1e-4 * mesh.bboxDiagSize2 * IPC_dt * IPC_dt));
        if (k && gradVanish && totalTimeStep > 1 - 1e-3) {
            // subproblem converged
            break;
            //return false;
        }

        calculateMovingDirection(mesh, BH, gradient, moveDir);

        //std::cout << "moveDir:" << std::endl;
        //for (int i = 0; i < 3; i++) {
        //    for (auto dir: moveDir) {
        //        std::cout << dir[i] << " ";
        //    }
        //    std::cout << std::endl;
        //}
        timer1.set_end();
        timer2.set_start();

        double alpha = 1.0, slackness_a = 0.8, slackness_m = 0.8;

        filterStepSize(mesh, moveDir, alpha);

        Environment_largestFeasibleStepSize(mesh, gd, moveDir, slackness_a, alpha);

        printf("env alpha:  %f\n", alpha);
        if (alpha < 0) {
            alpha = 1;
        }
        std::vector<std::pair<int, int>> newCandidates;

        Self_largestFeasibleStepSize(mesh, sh, moveDir, slackness_m, newCandidates, alpha);
        printf("partial self alpha:  %f\n", alpha);
        double partialCCD_alpha = alpha;

        Eigen::VectorXd pMag(mesh.surfVerts.size());
        for (int i = 0; i < pMag.size(); ++i) {
            int surfId = mesh.surfVerts[i];
            pMag[i] = moveDir[surfId].norm();
        }
        double alpha_CFL = std::sqrt(mesh.Hhat) / (pMag.maxCoeff() * 2.0);
        printf("alpha_CFL:  %f\n", alpha_CFL);
        //alpha = min(alpha, alpha_CFL);


        double fullCCD_alpha = alpha;
        sh.build(mesh, moveDir, fullCCD_alpha, mesh.averageEdgeLenth);
        Self_largestFeasibleStepSize_CCD(mesh, sh, moveDir, slackness_m, fullCCD_alpha);

        alpha = min(alpha, alpha_CFL);

        if (partialCCD_alpha > 2 * alpha_CFL) {
            alpha = min(partialCCD_alpha, fullCCD_alpha * 0.5);
            alpha = max(alpha, alpha_CFL);
        }

        cout << "CCD alpha:  " << fullCCD_alpha*0.5 << endl;
        timer2.set_end();
        timer3.set_start();

        printf("alpha:  %f\n", alpha);
        //}

        double alpha_feasible = alpha;

        bool isStop = lineSearch(mesh, sh, gd, moveDir, gradient, alpha, 0, 0, Kappa);
        timer3.set_end();
        timer4.set_start();
        postLineSearch(mesh, gd, alpha, Kappa);
        timer4.set_end();
        //cudaDeviceSynchronize();
        float time0, time1, time2, time3, time4;
        time0 = timer0.get_millisecond();
        time1 = timer1.get_millisecond();
        time2 = timer2.get_millisecond();
        time3 = timer3.get_millisecond();
        time4 = timer4.get_millisecond();
        totalTimeStep += alpha;
        if (k > 10 && isStop && totalTimeStep > 1 - 1e-3) {
            break;
        }


        printf("time0 = %f,  time1 = %f,  time2 = %f,  time3 = %f,  time4 = %f\n", time0, time1, time2, time3, time4);
    }
    printf("newton iteration:  %d    and    Kappa:  %f\n", k, mesh.Kappa);
    total_iter += k;
    return k;
}
void export_obj(const mesh3D& mesh, int index) {
    std::ofstream cloth_stream("output/surface_obj_" + std::to_string(index) + ".obj");
    cloth_stream << "# Generated by hc"
        << "\n";
    cloth_stream << std::fixed << std::setprecision(6) << "cloth\n";
    auto cloth_vertices = mesh.vertexes;
    auto cloth_triangles = mesh.surface;

    for (auto vec : cloth_vertices) {
        cloth_stream << "v " << vec(0) << " " << vec(1) << " "
            << vec(2) << "\n";
    }
    cloth_stream << "s 1\n";
    for (auto tri : cloth_triangles) {
        cloth_stream << "f " << tri(0) + 1 << " " << tri(1) + 1
            << " " << tri(2) + 1 << "\n";
    }
    cloth_stream.close();
}
int IPC_Solver(model_tet* meshTetes, SpatialHash& sh, Ground& gd) {
    std::chrono::steady_clock::time_point begin = std::chrono::steady_clock::now();
    mesh3D& mesh = meshTetes->mesh3Ds[0];
    
    upperBoundKappa(mesh.Kappa, mesh.Hhat, mesh.bboxDiagSize2, mesh.meanMass);
    if (mesh.Kappa < 1e-16) {
        suggestKappa(mesh.Kappa, mesh.Hhat, mesh.bboxDiagSize2, mesh.meanMass);
    }
    initKappa(mesh, gd, mesh.Kappa);
    
    //double dTol = 1e-18 * mesh.bboxDiagSize2;
    int k = 0;
    while (true) {
        mesh.closeConstraintID.resize(0);
        mesh.closeMConstraintID.resize(0);
        mesh.closeConstraintVal.resize(0);
        mesh.closeMConstraintVal.resize(0);
        k = solve_subIP(mesh, sh, gd, mesh.Kappa);

        VectorXd constraintVals;
        int offset = 0;//constraintVals.size();
        Evaluate_GroundConstraintVals(gd, mesh, constraintVals, offset);
        offset = constraintVals.size();
        Evaluate_SelfPTConstraintVals(mesh, constraintVals, offset);
        //offset = constraintVals.size();
        //Evaluate_SelfEEConstraintVals(mesh, constraintVals, offset);

        if (constraintVals.size()) {
            double minm = constraintVals.minCoeff();
            double maxm = constraintVals.maxCoeff();
            std::cout << minm << "    " << maxm << endl;
            if (constraintVals.minCoeff() < mesh.dTol) {
                break;
            }
            else if (constraintVals.maxCoeff() < mesh.Hhat) {
                break;
            }
        }
        else {
            break;
        }
    }

#ifdef USE_TBB
    tbb::parallel_for(0, mesh.vertexNum, 1, [&](int i)
#else
    for (int i = 0;i < mesh.vertexNum;i++)
#endif
    {
        mesh.velocities[i] = (mesh.vertexes[i] - mesh.V_prev[i]) / IPC_dt;
    }
#ifdef USE_TBB
    );
#endif

    mesh.V_prev = mesh.vertexes;
    computeXTilta(mesh);

    std::cout << "                                                finished a step" << endl;

    std::chrono::steady_clock::time_point end = std::chrono::steady_clock::now();
    std::cout << "Time for a step = "
        << std::chrono::duration_cast<std::chrono::microseconds>(end - begin).count() << "[us]" << std::endl;
    time_total += std::chrono::duration_cast<std::chrono::milliseconds>(end - begin).count();
    std::cout << "total  = "
        << time_total << "[ms]" << std::endl;
    std::cout << "total newton step:"<<total_iter << std::endl;
#ifdef EXPORT_OBJ
    std::cout << "frame number  = "
        << step_index++ << std::endl;
#endif // EXPORT_OBJ
    ofstream outTime("timeCost.txt");
    outTime << "time: " << time_total / 1000.0 << endl;
    outTime << "total iter: " << total_iter << endl;
    outTime << "frames: " << step_index << endl;
    outTime.close();
    //if (step_index == 250) {
    //    exit(0);
    //}
    return k;
}