cmake_minimum_required(VERSION 3.18)
project(cpu_ipc LANGUAGES CXX)

option(CIPC_BUILD_VIEWER "Build the GLUT/OpenGL viewer" ON)
option(CIPC_ENABLE_FRICTION "Enable lagged IPC friction" ON)
option(CIPC_ENABLE_QUADRATIC_BENDING "Enable quadratic cloth bending" ON)
option(CIPC_ENABLE_METIS_ORDERING "Enable CHOLMOD's cost-aware METIS ordering fallback" ON)
option(CIPC_ENABLE_PARDISO "Build the recommended/default oneMKL PARDISO backend when available" ON)

# Select a default build type for single-configuration generators.
if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
    set(CMAKE_BUILD_TYPE "Release" CACHE STRING "Build type (Debug or Release)" FORCE)
    set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS "Debug" "Release" "MinSizeRel" "RelWithDebInfo")
endif()

if(CMAKE_CONFIGURATION_TYPES)
    message(STATUS "Build configurations: ${CMAKE_CONFIGURATION_TYPES}")
else()
    message(STATUS "Build type: ${CMAKE_BUILD_TYPE}")
endif()
message(STATUS "Compiler: ${CMAKE_CXX_COMPILER_ID}")
message(STATUS "===============================================================")

# Include custom find modules
list(PREPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_LIST_DIR}/cmake")

# Find dependencies
find_package(Eigen3 REQUIRED)
find_package(TBB REQUIRED)
set(CIPC_HAS_PARDISO OFF)
if(CIPC_ENABLE_PARDISO)
    set(MKL_LINK static)
    set(MKL_THREADING tbb_thread)
    set(MKL_INTERFACE lp64)
    # oneMKL's config checks this component variable, while recent vcpkg TBB
    # configs provide TBB::tbb without defining it.
    if(TARGET TBB::tbb)
        set(TBB_tbb_FOUND TRUE)
    endif()
    find_package(MKL CONFIG QUIET)
    if(MKL_FOUND AND TARGET MKL::MKL)
        set(CIPC_HAS_PARDISO ON)
    else()
        message(STATUS "oneMKL was not found; PARDISO backend will be unavailable")
    endif()
endif()
set(CIPC_SUITESPARSE_COMPONENTS CHOLMOD LDL)
if(CIPC_ENABLE_METIS_ORDERING)
    list(APPEND CIPC_SUITESPARSE_COMPONENTS Partition)
endif()
find_package(SuiteSparse REQUIRED COMPONENTS ${CIPC_SUITESPARSE_COMPONENTS})
if(CIPC_ENABLE_METIS_ORDERING AND NOT TARGET SuiteSparse::Partition)
    message(FATAL_ERROR
        "CIPC_ENABLE_METIS_ORDERING=ON requires a CHOLMOD build with its "
        "Partition/METIS module and a linkable cholmod_metis implementation.")
endif()
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
    "${IPC_SOURCE_DIR}/SuiteSparseLDLSolver.cpp"
    "${IPC_SOURCE_DIR}/TetInversionGuard.cpp"
)
if(CIPC_HAS_PARDISO)
    if(MSVC)
        # The Windows vcpkg oneMKL package ships a Release CRT static runtime.
        # Keep ordinary /MDd Debug builds usable by excluding this optional
        # backend there; Release/RelWithDebInfo remain enabled.
        list(APPEND CIPC_CORE_SOURCES
            "$<$<NOT:$<CONFIG:Debug>>:${IPC_SOURCE_DIR}/PardisoSolver.cpp>")
    else()
        list(APPEND CIPC_CORE_SOURCES "${IPC_SOURCE_DIR}/PardisoSolver.cpp")
    endif()
endif()

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

# ContactMechanics.cpp contains large generated fixed-size kernels on MSVC.
target_compile_options(cipc_core PRIVATE
    $<$<CXX_COMPILER_ID:MSVC>:/bigobj>
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
if(CIPC_ENABLE_METIS_ORDERING)
    target_compile_definitions(cipc_core PRIVATE CIPC_CHOLMOD_USE_METIS)
endif()
if(CIPC_HAS_PARDISO)
    if(MSVC)
        target_compile_definitions(cipc_core PUBLIC
            "$<$<NOT:$<CONFIG:Debug>>:CIPC_HAS_PARDISO>")
    else()
        target_compile_definitions(cipc_core PUBLIC CIPC_HAS_PARDISO)
    endif()
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
    SuiteSparse::AMD
    SuiteSparse::CHOLMOD
    SuiteSparse::LDL
)
if(CIPC_ENABLE_METIS_ORDERING)
    target_link_libraries(cipc_core PUBLIC SuiteSparse::Partition)
endif()
if(CIPC_HAS_PARDISO)
    if(MSVC)
        target_link_libraries(cipc_core PRIVATE
            "$<$<NOT:$<CONFIG:Debug>>:MKL::MKL>")
    else()
        target_link_libraries(cipc_core PRIVATE MKL::MKL)
    endif()
endif()

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
if(CIPC_HAS_PARDISO AND MSVC)
    message(STATUS "  - Default linear solver: PARDISO (non-Debug); SuiteSparse LDL (Debug)")
    message(STATUS "  - oneMKL PARDISO backend: ON (non-Debug configurations)")
elseif(CIPC_HAS_PARDISO)
    message(STATUS "  - Default linear solver: PARDISO")
    message(STATUS "  - oneMKL PARDISO backend: ON")
else()
    message(STATUS "  - Default linear solver: SuiteSparse LDL")
    message(STATUS "  - oneMKL PARDISO backend: OFF")
endif()
message(STATUS "  - CHOLMOD METIS ordering: ${CIPC_ENABLE_METIS_ORDERING}")
message(STATUS "  - C++ Standard: 14")
message(STATUS "  - Eigen3: ${Eigen3_FOUND}")
message(STATUS "  - TBB: ${TBB_FOUND}")
message(STATUS "  - SuiteSparse: ${SuiteSparse_FOUND}")
if(CIPC_BUILD_VIEWER)
    message(STATUS "  - OpenGL: ${OpenGL_FOUND}")
    message(STATUS "  - GLUT: ${GLUT_FOUND}")
endif()
message(STATUS "===============================================================")
