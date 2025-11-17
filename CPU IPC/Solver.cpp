//
// Created by lamws on 2021/8/3.
//

#include "Solver.h"

#include <utility>
#include "iostream"


#include "tbb/spin_mutex.h"
#include "tbb/parallel_for.h"
#include "tbb/parallel_reduce.h"

using namespace std;

//#ifdef NDEBUG

void CholmodSolver::set_pattern(const SparseMatrix<double>& mtr) {
    numRows = static_cast<int>(mtr.rows());

    ja.conservativeResize(mtr.nonZeros());
    memcpy(ja.data(), mtr.innerIndexPtr(),
        mtr.nonZeros() * sizeof(mtr.innerIndexPtr()[0]));

    ia.conservativeResize(numRows + 1);
    memcpy(ia.data(), mtr.outerIndexPtr(),
        (numRows + 1) * sizeof(mtr.outerIndexPtr()[0]));

    a.conservativeResize(mtr.nonZeros());
    memcpy(a.data(), mtr.valuePtr(),
        mtr.nonZeros() * sizeof(mtr.valuePtr()[0]));
    if (!A) {
        A = cholmod_allocate_sparse(numRows, numRows, mtr.nonZeros(),
            true, true, -1, CHOLMOD_REAL, &cm);
        Ax = A->x;
        Ap = A->p;
        Ai = A->i;
        // -1: upper right part will be ignored during computation

        A->i = ja.data();
        A->p = ia.data();
        A->x = a.data();
    }
}

void CholmodSolver::solve(VectorXd& rhs, VectorXd& result) {
    if (!b) {
        b = cholmod_allocate_dense(numRows, 1, numRows, CHOLMOD_REAL, &cm);
        bx = b->x;
    }
    b->x = rhs.data();
    L = cholmod_analyze(A, &cm);
    cholmod_factorize(A, L, &cm);
    cholmod_dense* x;
    x = cholmod_solve(CHOLMOD_A, L, b, &cm);
    result.conservativeResize(rhs.size());
    memcpy(result.data(), x->x, result.size() * sizeof(result[0]));
    cholmod_free_dense(&x, &cm);
}


void CholmodSolver::preFactorize(const SparseMatrix<double>& mtr) {
    numRows = static_cast<int>(mtr.rows());

    ja.conservativeResize(mtr.nonZeros());
    memcpy(ja.data(), mtr.innerIndexPtr(),
        mtr.nonZeros() * sizeof(mtr.innerIndexPtr()[0]));

    ia.conservativeResize(numRows + 1);
    memcpy(ia.data(), mtr.outerIndexPtr(),
        (numRows + 1) * sizeof(mtr.outerIndexPtr()[0]));

    a.conservativeResize(mtr.nonZeros());
    memcpy(a.data(), mtr.valuePtr(),
        mtr.nonZeros() * sizeof(mtr.valuePtr()[0]));
    if (!A) {
        A = cholmod_allocate_sparse(numRows, numRows, mtr.nonZeros(),
            true, true, -1, CHOLMOD_REAL, &cm);
        Ax = A->x;
        Ap = A->p;
        Ai = A->i;
        // -1: upper right part will be ignored during computation

        A->i = ja.data();
        A->p = ia.data();
        A->x = a.data();
    }

    L = cholmod_analyze(A, &cm);
    cholmod_factorize(A, L, &cm);
}

void CholmodSolver::solve_with_preFactorize(VectorXd& rhs, VectorXd& result) {
    if (!b) {
        b = cholmod_allocate_dense(numRows, 1, numRows, CHOLMOD_REAL, &cm);
        bx = b->x;
    }
    b->x = rhs.data();

    cholmod_dense* x;
    x = cholmod_solve(CHOLMOD_A, L, b, &cm);
    result.conservativeResize(rhs.size());
    memcpy(result.data(), x->x, result.size() * sizeof(result[0]));
    cholmod_free_dense(&x, &cm);
}


CholmodSolver::CholmodSolver() {
    cholmod_start(&cm);
    A = nullptr;
    L = nullptr;
    b = nullptr;
    x_cd = y_cd = nullptr;

    Ai = Ap = Ax = nullptr;
    bx = nullptr;
    solutionx = x_cdx = y_cdx = nullptr;
}

CholmodSolver::~CholmodSolver() {
    if (A) {
        A->i = Ai;
        A->p = Ap;
        A->x = Ax;
        cholmod_free_sparse(&A, &cm);
    }

    cholmod_free_factor(&L, &cm);

    if (b) {
        b->x = bx;
        cholmod_free_dense(&b, &cm);
    }

    if (x_cd) {
        x_cd->x = x_cdx;
        cholmod_free_dense(&x_cd, &cm);
    }

    if (y_cd) {
        y_cd->x = y_cdx;
        cholmod_free_dense(&y_cd, &cm);
    }

    cholmod_finish(&cm);
}
//#endif