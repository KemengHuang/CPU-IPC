cmake_minimum_required(VERSION 3.18)

# ------------------------------------------------------------------------------
# Project
# ------------------------------------------------------------------------------
project(cpu_ipc LANGUAGES CXX)

# ------------------------------------------------------------------------------
# Default build type for single-configuration generators
# ------------------------------------------------------------------------------
if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
    set(CMAKE_BUILD_TYPE "Release" CACHE STRING "Build type" FORCE)
    set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS "Debug" "Release" "MinSizeRel" "RelWithDebInfo")
endif()

# ------------------------------------------------------------------------------
# vcpkg toolchain hint
# ------------------------------------------------------------------------------
if(NOT CMAKE_TOOLCHAIN_FILE)
    message(WARNING
        "CMAKE_TOOLCHAIN_FILE is not set. If you use vcpkg, pass "
        "-DCMAKE_TOOLCHAIN_FILE=<vcpkg-root>/scripts/buildsystems/vcpkg.cmake "
        "when configuring for the first time.")
endif()

if(CMAKE_CONFIGURATION_TYPES)
    message(STATUS "Build configurations: ${CMAKE_CONFIGURATION_TYPES}")
else()
    message(STATUS "Build type: ${CMAKE_BUILD_TYPE}")
endif()
message(STATUS "Compiler: ${CMAKE_CXX_COMPILER_ID}")
message(STATUS "===============================================================")

# ------------------------------------------------------------------------------
# User-configurable options
# ------------------------------------------------------------------------------
option(CIPC_USE_SNK "Use Stable Neo-Hookean (SNK) elasticity model" ON)
option(CIPC_USE_FRICTION "Enable friction support" ON)
option(CIPC_USE_QUADRATIC_BENDING "Enable quadratic bending energy" ON)
option(CIPC_USE_NEWBARRIER "Use the new log-barrier formulation (NEWB)" OFF)
option(CIPC_BUILD_TESTS "Build optional unit/verification tests" OFF)

# ------------------------------------------------------------------------------
# Custom find modules
# ------------------------------------------------------------------------------
list(PREPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_LIST_DIR}/cmake")

# ------------------------------------------------------------------------------
# Dependencies
# ------------------------------------------------------------------------------
find_package(Eigen3 REQUIRED)
find_package(TBB REQUIRED)
find_package(SuiteSparse REQUIRED COMPONENTS CHOLMOD)
find_package(OpenGL REQUIRED)
find_package(GLUT REQUIRED)
find_package(GLEW REQUIRED)

# ------------------------------------------------------------------------------
# Executable target
# ------------------------------------------------------------------------------
add_executable(cipc)

set(CIPC_SOURCE_DIR "${CMAKE_CURRENT_SOURCE_DIR}/src")

file(GLOB_RECURSE CIPC_SOURCES
    CONFIGURE_DEPENDS
    "${CIPC_SOURCE_DIR}/*.cpp"
)
list(FILTER CIPC_SOURCES EXCLUDE REGEX "_test\\.cpp$")

file(GLOB_RECURSE CIPC_HEADERS
    CONFIGURE_DEPENDS
    "${CIPC_SOURCE_DIR}/*.h"
    "${CIPC_SOURCE_DIR}/*.hpp"
    "${CIPC_SOURCE_DIR}/*.inl"
    "${CIPC_SOURCE_DIR}/*.inc"
    "${CIPC_SOURCE_DIR}/*.cuh"
)

target_sources(cipc PRIVATE
    ${CIPC_SOURCES}
    ${CIPC_HEADERS}
)

# ------------------------------------------------------------------------------
# C++ standard
# ------------------------------------------------------------------------------
set_target_properties(cipc PROPERTIES
    CXX_STANDARD 14
    CXX_STANDARD_REQUIRED ON
    CXX_EXTENSIONS OFF
)

# ------------------------------------------------------------------------------
# Compile definitions
# ------------------------------------------------------------------------------
target_compile_definitions(cipc PRIVATE
    CIPC_ASSETS_DIR="${CMAKE_CURRENT_SOURCE_DIR}/Assets/"
    CIPC_OUTPUT_DIR="${CMAKE_CURRENT_SOURCE_DIR}/Output/"
)

if(CIPC_USE_SNK)
    target_compile_definitions(cipc PRIVATE USE_SNK)
endif()

if(CIPC_USE_FRICTION)
    target_compile_definitions(cipc PRIVATE USE_FRICTION)
endif()

if(CIPC_USE_QUADRATIC_BENDING)
    target_compile_definitions(cipc PRIVATE USE_QUADRATIC_BENDING)
endif()

if(CIPC_USE_NEWBARRIER)
    target_compile_definitions(cipc PRIVATE NEWB)
endif()

# ------------------------------------------------------------------------------
# Include directories
# ------------------------------------------------------------------------------
target_include_directories(cipc PRIVATE
    "${CIPC_SOURCE_DIR}"
    "${CMAKE_CURRENT_SOURCE_DIR}"
)

# ------------------------------------------------------------------------------
# Compiler flags
# ------------------------------------------------------------------------------
if(MSVC)
    target_compile_options(cipc PRIVATE /bigobj)
endif()

if(UNIX)
    target_compile_options(cipc PRIVATE $<$<CONFIG:Release>:-O3>)
endif()

# ------------------------------------------------------------------------------
# Linked libraries
# ------------------------------------------------------------------------------
target_link_libraries(cipc PRIVATE
    Eigen3::Eigen
    TBB::tbb
    SuiteSparse::CHOLMOD
    OpenGL::GLU
    GLUT::GLUT
    GLEW::GLEW
)

# ------------------------------------------------------------------------------
# IDE organization
# ------------------------------------------------------------------------------
source_group(TREE "${CMAKE_CURRENT_SOURCE_DIR}"
    PREFIX "Source Files"
    FILES ${CIPC_SOURCES} ${CIPC_HEADERS}
)

# ------------------------------------------------------------------------------
# Configuration summary
# ------------------------------------------------------------------------------
message(STATUS "Configuration Summary:")
message(STATUS "  - Sources: ${CIPC_SOURCE_DIR}")
message(STATUS "  - C++ Standard: 14")
message(STATUS "  - Eigen3: ${Eigen3_FOUND}")
message(STATUS "  - TBB: ${TBB_FOUND}")
message(STATUS "  - SuiteSparse: ${SuiteSparse_FOUND}")
message(STATUS "  - OpenGL: ${OpenGL_FOUND}")
message(STATUS "  - GLUT: ${GLUT_FOUND}")
message(STATUS "  - GLEW: ${GLEW_FOUND}")
message(STATUS "  - USE_SNK: ${CIPC_USE_SNK}")
message(STATUS "  - USE_FRICTION: ${CIPC_USE_FRICTION}")
message(STATUS "  - USE_QUADRATIC_BENDING: ${CIPC_USE_QUADRATIC_BENDING}")
message(STATUS "  - USE_NEWBARRIER: ${CIPC_USE_NEWBARRIER}")
message(STATUS "  - BUILD_TESTS: ${CIPC_BUILD_TESTS}")
message(STATUS "===============================================================")

# ------------------------------------------------------------------------------
# Optional verification tests
# ------------------------------------------------------------------------------
if(CIPC_BUILD_TESTS)
    add_executable(bvh_test src/collision/BVH_test.cpp)
    set_target_properties(bvh_test PROPERTIES
        CXX_STANDARD 14
        CXX_STANDARD_REQUIRED ON
        CXX_EXTENSIONS OFF
    )
    target_include_directories(bvh_test PRIVATE
        "${CIPC_SOURCE_DIR}"
        "${CMAKE_CURRENT_SOURCE_DIR}"
    )
    target_link_libraries(bvh_test PRIVATE Eigen3::Eigen)
endif()
