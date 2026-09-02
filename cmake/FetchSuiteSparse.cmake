cmake_minimum_required(VERSION 3.23)

if(NOT CPU_IPC_SUITESPARSE_SOURCE_PARENT)
    message(FATAL_ERROR "CPU_IPC_SUITESPARSE_SOURCE_PARENT is required")
endif()
if(NOT CPU_IPC_DOWNLOAD_DIR)
    message(FATAL_ERROR "CPU_IPC_DOWNLOAD_DIR is required")
endif()

file(READ "${CMAKE_CURRENT_LIST_DIR}/../scripts/suitesparse-version.txt"
    CPU_IPC_SUITESPARSE_VERSION)
string(STRIP "${CPU_IPC_SUITESPARSE_VERSION}" CPU_IPC_SUITESPARSE_VERSION)

set(_CPU_IPC_SUITESPARSE_SOURCE
    "${CPU_IPC_SUITESPARSE_SOURCE_PARENT}/SuiteSparse-${CPU_IPC_SUITESPARSE_VERSION}")
if(EXISTS "${_CPU_IPC_SUITESPARSE_SOURCE}/CMakeLists.txt")
    message(STATUS "Reusing SuiteSparse source: ${_CPU_IPC_SUITESPARSE_SOURCE}")
else()
    file(MAKE_DIRECTORY
        "${CPU_IPC_SUITESPARSE_SOURCE_PARENT}"
        "${CPU_IPC_DOWNLOAD_DIR}")
    set(_CPU_IPC_SUITESPARSE_ARCHIVE
        "${CPU_IPC_DOWNLOAD_DIR}/DrTimothyAldenDavis-SuiteSparse-v${CPU_IPC_SUITESPARSE_VERSION}.tar.gz")
    set(_CPU_IPC_SUITESPARSE_URL
        "https://github.com/DrTimothyAldenDavis/SuiteSparse/archive/refs/tags/v${CPU_IPC_SUITESPARSE_VERSION}.tar.gz")

    message(STATUS "Downloading SuiteSparse ${CPU_IPC_SUITESPARSE_VERSION}")
    file(DOWNLOAD
        "${_CPU_IPC_SUITESPARSE_URL}"
        "${_CPU_IPC_SUITESPARSE_ARCHIVE}"
        EXPECTED_HASH
            "SHA512=4cd00b0625ef8081703139cf15c45c32a90e101bcd9c6ba38b87b988c2e76f51a73da0c312581b87e0892098c717d14782edfd6f19c92b0b0db196044cfe8337"
        STATUS _CPU_IPC_DOWNLOAD_STATUS
        SHOW_PROGRESS
        TLS_VERIFY ON)
    list(GET _CPU_IPC_DOWNLOAD_STATUS 0 _CPU_IPC_DOWNLOAD_CODE)
    list(GET _CPU_IPC_DOWNLOAD_STATUS 1 _CPU_IPC_DOWNLOAD_MESSAGE)
    if(NOT _CPU_IPC_DOWNLOAD_CODE EQUAL 0)
        message(FATAL_ERROR
            "SuiteSparse download failed: ${_CPU_IPC_DOWNLOAD_MESSAGE}")
    endif()

    message(STATUS "Extracting SuiteSparse ${CPU_IPC_SUITESPARSE_VERSION}")
    file(ARCHIVE_EXTRACT
        INPUT "${_CPU_IPC_SUITESPARSE_ARCHIVE}"
        DESTINATION "${CPU_IPC_SUITESPARSE_SOURCE_PARENT}")
    if(NOT EXISTS "${_CPU_IPC_SUITESPARSE_SOURCE}/CMakeLists.txt")
        message(FATAL_ERROR
            "SuiteSparse archive did not create ${_CPU_IPC_SUITESPARSE_SOURCE}")
    endif()
endif()

# SuiteSparse's upstream BLAS discovery expects vendor metadata not exported by
# oneMKL's CMake target. Use the same minimal provider shim as the vcpkg port.
configure_file(
    "${CMAKE_CURRENT_LIST_DIR}/SuiteSparseBLAS.cmake"
    "${_CPU_IPC_SUITESPARSE_SOURCE}/SuiteSparse_config/cmake_modules/SuiteSparseBLAS.cmake"
    COPYONLY)
