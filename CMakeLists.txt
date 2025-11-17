cmake_minimum_required(VERSION 3.18)
project(cpu_ipc LANGUAGES CXX)

# Set build type before project() for multi-configuration generators
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
find_package(OpenGL REQUIRED)
find_package(GLUT REQUIRED)
find_package(GLEW REQUIRED)

# Create executable
add_executable(cipc)

# Configure source files
set(IPC_SOURCE_DIR "CPU IPC")

file(GLOB_RECURSE CIPC_SOURCES
    CONFIGURE_DEPENDS
    "${IPC_SOURCE_DIR}/*.cpp"
)

file(GLOB_RECURSE CIPC_HEADERS
    CONFIGURE_DEPENDS
    "${IPC_SOURCE_DIR}/*.h"
    "${IPC_SOURCE_DIR}/*.hpp"
    "${IPC_SOURCE_DIR}/*.inl"
    "${IPC_SOURCE_DIR}/*.inc"
    "${IPC_SOURCE_DIR}/*.cuh"
)

# Add sources to target
target_sources(cipc PRIVATE 
    ${CIPC_SOURCES}
    ${CIPC_HEADERS}  # Include headers for IDE organization
)

# Set C++ standard
set_target_properties(cipc PROPERTIES
    CXX_STANDARD 14
    CXX_STANDARD_REQUIRED ON
    CXX_EXTENSIONS OFF
)

target_compile_definitions(cipc PRIVATE CIPC_ASSETS_DIR="${CMAKE_CURRENT_SOURCE_DIR}/Assets/")
target_compile_definitions(cipc PRIVATE CIPC_OUTPUT_DIR="${CMAKE_CURRENT_SOURCE_DIR}/Output/")

# Compile definitions
target_compile_definitions(cipc PRIVATE 
    USE_SNK
    USE_FRICTION
    USE_QUADRATIC_BENDING
)

# Include directories
target_include_directories(cipc PRIVATE 
    "${IPC_SOURCE_DIR}"
    "${CMAKE_CURRENT_SOURCE_DIR}"
)

# Link libraries - consolidated into one call
target_link_libraries(cipc PRIVATE
    Eigen3::Eigen
    TBB::tbb
    
    # SuiteSparse components
    SuiteSparse::CHOLMOD
    #SuiteSparse::AMD
    #SuiteSparse::CAMD
    #SuiteSparse::CCOLAMD
    #SuiteSparse::COLAMD

    OpenGL::GLU
    GLUT::GLUT
    GLEW::GLEW
)

# Organize source files in IDEs
source_group(TREE "${CMAKE_CURRENT_SOURCE_DIR}" 
    PREFIX "Source Files"
    FILES ${CIPC_SOURCES} ${CIPC_HEADERS}
)

# Print configuration summary
message(STATUS "Configuration Summary:")
message(STATUS "  - Sources: ${IPC_SOURCE_DIR}")
message(STATUS "  - C++ Standard: 14")
message(STATUS "  - Eigen3: ${Eigen3_FOUND}")
message(STATUS "  - TBB: ${TBB_FOUND}")
message(STATUS "  - SuiteSparse: ${SuiteSparse_FOUND}")
message(STATUS "  - OpenGL: ${OpenGL_FOUND}")
message(STATUS "  - GLUT: ${GLUT_FOUND}")
message(STATUS "  - GLEW: ${GLEW_FOUND}")
message(STATUS "===============================================================")