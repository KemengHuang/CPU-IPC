cmake_minimum_required(VERSION 3.18)    # support relative path, cuda_std_17, archive manipulation

message("===============================================================")
if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE "Release" CACHE STRING "Choose Release or Debug")
endif()

# project meta
# ============
project(cpu_ipc
    LANGUAGES   CXX)

set(srcs 
    "CPU IPC/ACCD.cpp"
    "CPU IPC/collisionUtil.cpp"
    "CPU IPC/fem_math.cpp"
    "CPU IPC/fem3D.cpp"
    "CPU IPC/FEMTimeIntegrator.cpp"
    "CPU IPC/gl_main.cpp"
    "CPU IPC/IPC_FUNC.cpp"
    "CPU IPC/IPCdistanceFuncs.cpp"
    "CPU IPC/IPCtimeStepFuncs.cpp"
    "CPU IPC/mesh.cpp"
    "CPU IPC/mIPC.cpp"
    "CPU IPC/Simulator.cpp"
    "CPU IPC/Solver.cpp"
    "CPU IPC/SVD.cpp"
)
add_executable(gipc ${srcs})
target_compile_features(gipc PRIVATE cxx_std_17)
add_compile_definitions(gipc PRIVATE USE_SNK NDEBUG)

find_package(Eigen3 REQUIRED)
find_package(TBB REQUIRED)
target_link_directories(gipc PRIVATE /usr/local/lib /opt/intel/oneapi/mkl/2022.0.2/lib/intel64 /opt/intel/oneapi/compiler/2022.0.2/linux/compiler/lib/intel64_lin)
target_link_libraries(gipc PRIVATE Eigen3::Eigen TBB::tbb)
target_link_libraries(gipc PRIVATE mkl_intel_lp64 mkl_intel_thread mkl_core iomp5 pthread m dl mkl_lapack95_lp64 mkl_blas95_lp64 cholmod)

if (MSVC)
    message("win build")
    target_include_directories(gipc PRIVATE 
        mLBVH/GL_LIB/glew_64/include
        mLBVH/GL_LIB/freeglut_64/include
        mLBVH)
    target_link_directories(gipc PRIVATE 
        mLBVH/GL_LIB/glew_64/libs/release
        mLBVH/GL_LIB/freeglut_64/libs/release
    )
    target_link_libraries(gipc PRIVATE
        glew32.lib
        freeglut.lib
    )
elseif(UNIX)
    find_package(OpenGL REQUIRED)
    find_package(GLUT REQUIRED)
    find_package(GLEW REQUIRED)
    target_link_libraries(gipc PRIVATE OpenGL::OpenGL GLUT::GLUT GLEW::GLEW OpenGL::GLU)
endif(MSVC)
