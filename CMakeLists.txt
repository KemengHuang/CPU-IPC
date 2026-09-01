cmake_minimum_required(VERSION 3.18)
project(cpu_ipc LANGUAGES CXX)

option(CIPC_BUILD_VIEWER "Build the GLUT/OpenGL viewer" ON)
option(CIPC_ENABLE_FRICTION "Enable lagged IPC friction" ON)
option(CIPC_ENABLE_QUADRATIC_BENDING "Enable quadratic cloth bending" ON)

# Select a default build type for single-configuration generators.
if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
    set(CMAKE_BUILD_TYPE "Release" CACHE STRING "Build type (Debug or Release)" FORCE)
    set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS "Debug" "Release" "MinSizeRel" "RelWithDebInfo")
endif()

# General compiler flags
if(MSVC)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /bigobj")
elseif(UNIX AND CMAKE_BUILD_TYPE STREQUAL "Release")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3")
endif()

message(STATUS "Build type: ${CMAKE_BUILD_TYPE}")
message(STATUS "Compiler: ${CMAKE_CXX_COMPILER_ID}")
message(STATUS "===============================================================")

# Include custom find modules
list(PREPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_LIST_DIR}/cmake")

# Find dependencies
find_package(Eigen3 REQUIRED)
find_package(TBB REQUIRED)
find_package(SuiteSparse REQUIRED CHOLMOD)
if(CIPC_BUILD_VIEWER)
    find_package(OpenGL REQUIRED)
    find_package(GLUT REQUIRED)
endif()

# Configure source files
set(IPC_SOURCE_DIR "CPU IPC")

set(CIPC_CORE_SOURCES
    "${IPC_SOURCE_DIR}/AdditiveCCD.cpp"
    "${IPC_SOURCE_DIR}/CholmodSolver.cpp"
    "${IPC_SOURCE_DIR}/CollisionBroadPhase.cpp"
    "${IPC_SOURCE_DIR}/ContactMechanics.cpp"
    "${IPC_SOURCE_DIR}/Elasticity.cpp"
    "${IPC_SOURCE_DIR}/ElasticityMath.cpp"
    "${IPC_SOURCE_DIR}/FeasibleStep.cpp"
    "${IPC_SOURCE_DIR}/Friction.cpp"
    "${IPC_SOURCE_DIR}/HingeBending.cpp"
    "${IPC_SOURCE_DIR}/IPCSolver.cpp"
    "${IPC_SOURCE_DIR}/LBVH.cpp"
    "${IPC_SOURCE_DIR}/NewtonLinearSystem.cpp"
    "${IPC_SOURCE_DIR}/RotationAwareSVD.cpp"
    "${IPC_SOURCE_DIR}/RuntimePaths.cpp"
    "${IPC_SOURCE_DIR}/SimulationMesh.cpp"
    "${IPC_SOURCE_DIR}/Simulator.cpp"
    "${IPC_SOURCE_DIR}/TetInversionGuard.cpp"
)

file(GLOB_RECURSE CIPC_HEADERS
    CONFIGURE_DEPENDS
    "${IPC_SOURCE_DIR}/*.h"
    "${IPC_SOURCE_DIR}/*.hpp"
    "${IPC_SOURCE_DIR}/*.inl"
    "${IPC_SOURCE_DIR}/*.inc"
    "${IPC_SOURCE_DIR}/*.cuh"
)

set(CIPC_VIEWER_SOURCE "${CMAKE_CURRENT_SOURCE_DIR}/${IPC_SOURCE_DIR}/ViewerMain.cpp")

add_library(cipc_core STATIC
    ${CIPC_CORE_SOURCES}
    ${CIPC_HEADERS}
)

add_executable(cipc_headless apps/cipc_headless.cpp)
if(CIPC_BUILD_VIEWER)
    add_executable(cipc "${CIPC_VIEWER_SOURCE}")
endif()

# Set C++ standard
set(CIPC_PRODUCT_TARGETS cipc_core cipc_headless)
if(CIPC_BUILD_VIEWER)
    list(APPEND CIPC_PRODUCT_TARGETS cipc)
endif()
set_target_properties(${CIPC_PRODUCT_TARGETS} PROPERTIES
    CXX_STANDARD 14
    CXX_STANDARD_REQUIRED ON
    CXX_EXTENSIONS OFF
)

target_compile_definitions(cipc_core PUBLIC
    CIPC_ASSETS_DIR="${CMAKE_CURRENT_SOURCE_DIR}/Assets/"
    CIPC_OUTPUT_DIR="${CMAKE_CURRENT_SOURCE_DIR}/Output/"
)

# Compile definitions
if(CIPC_ENABLE_FRICTION)
    target_compile_definitions(cipc_core PUBLIC USE_FRICTION)
endif()
if(CIPC_ENABLE_QUADRATIC_BENDING)
    target_compile_definitions(cipc_core PUBLIC USE_QUADRATIC_BENDING)
endif()

# Include directories
target_include_directories(cipc_core PUBLIC
    "${IPC_SOURCE_DIR}"
    "${CMAKE_CURRENT_SOURCE_DIR}"
)

# Link libraries - consolidated into one call
target_link_libraries(cipc_core PUBLIC
    Eigen3::Eigen
    TBB::tbb
    
    # SuiteSparse components
    SuiteSparse::CHOLMOD
    #SuiteSparse::AMD
    #SuiteSparse::CAMD
    #SuiteSparse::CCOLAMD
    #SuiteSparse::COLAMD
)

if(CIPC_BUILD_VIEWER)
    target_link_libraries(cipc PRIVATE
        cipc_core
        OpenGL::GLU
        GLUT::GLUT
    )
endif()

target_link_libraries(cipc_headless PRIVATE cipc_core)

# Organize source files in IDEs
source_group(TREE "${CMAKE_CURRENT_SOURCE_DIR}"
    PREFIX "Source Files"
    FILES ${CIPC_CORE_SOURCES} ${CIPC_VIEWER_SOURCE} ${CIPC_HEADERS} apps/cipc_headless.cpp
)

# Print configuration summary
message(STATUS "Configuration Summary:")
message(STATUS "  - Core library: cipc_core")
message(STATUS "  - Viewer: ${CIPC_BUILD_VIEWER}")
message(STATUS "  - Headless runner: cipc_headless")
message(STATUS "  - Friction: ${CIPC_ENABLE_FRICTION}")
message(STATUS "  - Quadratic bending: ${CIPC_ENABLE_QUADRATIC_BENDING}")
message(STATUS "  - C++ Standard: 14")
message(STATUS "  - Eigen3: ${Eigen3_FOUND}")
message(STATUS "  - TBB: ${TBB_FOUND}")
message(STATUS "  - SuiteSparse: ${SuiteSparse_FOUND}")
if(CIPC_BUILD_VIEWER)
    message(STATUS "  - OpenGL: ${OpenGL_FOUND}")
    message(STATUS "  - GLUT: ${GLUT_FOUND}")
endif()
message(STATUS "===============================================================")
