# Loaded by the standalone CHOLMOD build after project() enables C/C++.
# Keep CHOLMOD's BLAS/LAPACK ABI aligned with the project's oneMKL LP64/TBB
# configuration.
set(MKL_LINK static)
set(MKL_THREADING tbb_thread)
set(MKL_INTERFACE lp64)

find_package(TBB CONFIG REQUIRED)
set(TBB_tbb_FOUND TRUE)
find_package(MKL CONFIG REQUIRED)

set(BLAS_LIBRARIES MKL::MKL CACHE STRING "oneMKL BLAS provider" FORCE)
set(LAPACK_LIBRARIES MKL::MKL CACHE STRING "oneMKL LAPACK provider" FORCE)
set(BLAS_FOUND TRUE CACHE BOOL "oneMKL BLAS found" FORCE)
set(LAPACK_FOUND TRUE CACHE BOOL "oneMKL LAPACK found" FORCE)
