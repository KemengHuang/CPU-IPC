# Minimal provider shim used by the pinned SuiteSparse source build. The
# surrounding project preloads BLAS_LIBRARIES with MKL::MKL.
if(NOT BLAS_FOUND OR NOT BLAS_LIBRARIES)
    find_package(BLAS REQUIRED)
endif()
set(BLA_SIZEOF_INTEGER 4)
set(SuiteSparse_BLAS_integer int32_t)

if(WIN32)
    add_compile_definitions(BLAS64__SUFFIX=_)
endif()
