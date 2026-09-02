#!/usr/bin/env bash
set -euo pipefail

repo_root="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
build_root="${repo_root}/build-wsl"
configuration="Release"
build_viewer="OFF"
requested_vcpkg_root="${VCPKG_ROOT:-}"
cache_root="${XDG_CACHE_HOME:-${HOME}/.cache}/cpu-ipc"

while [[ $# -gt 0 ]]; do
    case "$1" in
        --vcpkg-root)
            if [[ $# -lt 2 ]]; then
                echo "--vcpkg-root requires a path" >&2
                exit 2
            fi
            requested_vcpkg_root="$2"
            shift 2
            ;;
        --build-dir)
            if [[ $# -lt 2 ]]; then
                echo "--build-dir requires a path" >&2
                exit 2
            fi
            build_root="$2"
            shift 2
            ;;
        --config)
            if [[ $# -lt 2 ]]; then
                echo "--config requires Release, RelWithDebInfo, or Debug" >&2
                exit 2
            fi
            configuration="$2"
            shift 2
            ;;
        --viewer)
            build_viewer="ON"
            shift
            ;;
        --headless-only)
            build_viewer="OFF"
            shift
            ;;
        -h|--help)
            cat <<'EOF'
Usage: ./build.sh [options]
  --vcpkg-root PATH  Use an existing Linux vcpkg checkout.
  --build-dir PATH   Build root (default: ./build-wsl).
  --config TYPE      Release, RelWithDebInfo, or Debug.
  --viewer           Build the GLUT viewer (headless is the default).
  --headless-only    Explicitly disable the viewer.
EOF
            exit 0
            ;;
        *)
            echo "Unknown option: $1" >&2
            exit 2
            ;;
    esac
done

case "$configuration" in
    Release|RelWithDebInfo|Debug) ;;
    *)
        echo "Unsupported configuration: $configuration" >&2
        exit 2
        ;;
esac

if [[ "$(uname -m)" != "x86_64" ]]; then
    echo "This oneMKL build currently supports x86_64 Ubuntu/WSL only." >&2
    exit 1
fi

required_commands=(git cmake ninja g++ make tar curl)
missing_commands=()
for command_name in "${required_commands[@]}"; do
    if ! command -v "$command_name" >/dev/null 2>&1; then
        missing_commands+=("$command_name")
    fi
done
if [[ ${#missing_commands[@]} -ne 0 ]]; then
    echo "Missing Ubuntu build tools: ${missing_commands[*]}" >&2
    echo "Install them once with:" >&2
    echo "  sudo apt-get update && sudo apt-get install -y build-essential cmake ninja-build git curl tar" >&2
    exit 1
fi

# vcpkg's bootstrap requires zip/unzip. On the tested minimal WSL image these
# were absent, so obtain just those tools in the user cache without sudo.
user_tools="${cache_root}/tools"
if ! command -v zip >/dev/null 2>&1 || ! command -v unzip >/dev/null 2>&1; then
    if ! command -v apt-get >/dev/null 2>&1 || ! command -v dpkg-deb >/dev/null 2>&1; then
        echo "zip/unzip are missing and apt-get/dpkg-deb are unavailable." >&2
        exit 1
    fi
    if [[ ! -x "${user_tools}/usr/bin/zip" || ! -x "${user_tools}/usr/bin/unzip" ]]; then
        package_cache="${cache_root}/apt-packages"
        mkdir -p "$package_cache" "$user_tools"
        (
            cd "$package_cache"
            apt-get download zip unzip
        )
        for package in "${package_cache}"/*.deb; do
            dpkg-deb -x "$package" "$user_tools"
        done
    fi
    export PATH="${user_tools}/usr/bin:${PATH}"
fi

if [[ -n "$requested_vcpkg_root" ]]; then
    if [[ -x "$requested_vcpkg_root" && ! -d "$requested_vcpkg_root" ]]; then
        vcpkg_executable="$(readlink -f "$requested_vcpkg_root")"
        vcpkg_root="$(dirname "$vcpkg_executable")"
    elif [[ -d "$requested_vcpkg_root" ]]; then
        vcpkg_root="$(cd "$requested_vcpkg_root" && pwd)"
    else
        echo "The requested vcpkg path does not exist: $requested_vcpkg_root" >&2
        exit 1
    fi
elif command -v vcpkg >/dev/null 2>&1; then
    vcpkg_executable="$(readlink -f "$(command -v vcpkg)")"
    vcpkg_root="$(dirname "$vcpkg_executable")"
else
    vcpkg_root="${cache_root}/vcpkg"
    if [[ ! -x "${vcpkg_root}/vcpkg" ]]; then
        mkdir -p "$cache_root"
        if [[ ! -d "${vcpkg_root}/.git" ]]; then
            git clone --depth 1 https://github.com/microsoft/vcpkg.git "$vcpkg_root"
        fi
        "${vcpkg_root}/bootstrap-vcpkg.sh" -disableMetrics
    fi
fi

mkdir -p "$build_root"
build_root="$(cd "$build_root" && pwd)"

vcpkg_executable="${vcpkg_root}/vcpkg"
vcpkg_toolchain="${vcpkg_root}/scripts/buildsystems/vcpkg.cmake"
if [[ ! -x "$vcpkg_executable" || ! -f "$vcpkg_toolchain" ]]; then
    echo "Invalid Linux vcpkg installation: $vcpkg_root" >&2
    exit 1
fi

triplet="x64-linux"
dependencies=(
    "eigen3:${triplet}"
    "tbb:${triplet}"
    "metis:${triplet}"
    "intel-mkl:${triplet}"
    "suitesparse-cholmod[partition]:${triplet}"
)
if [[ "$build_viewer" == "ON" ]]; then
    dependencies+=("freeglut:${triplet}")
fi
"$vcpkg_executable" install "${dependencies[@]}" --recurse

source_parent="${vcpkg_root}/buildtrees/suitesparse-cholmod/src"
cholmod_source=""
for candidate in "${source_parent}"/*.clean; do
    if [[ ! -f "${candidate}/CHOLMOD/CMakeLists.txt" ]]; then
        continue
    fi
    if [[ -z "$cholmod_source" || "$candidate" -nt "${cholmod_source%/CHOLMOD}" ]]; then
        cholmod_source="${candidate}/CHOLMOD"
    fi
done
if [[ ! -f "${cholmod_source}/CMakeLists.txt" ]]; then
    echo "Extracted CHOLMOD source was not found under $source_parent" >&2
    exit 1
fi

installed_prefix="${vcpkg_root}/installed/${triplet}"
cholmod_build="${build_root}/cholmod-mkl-lib"
cholmod_install="${build_root}/cholmod-mkl-install"
provider_file="${repo_root}/cmake/CholmodMKLProvider.cmake"

cmake -S "$cholmod_source" -B "$cholmod_build" -G Ninja \
    -DCMAKE_BUILD_TYPE=Release \
    -DCMAKE_PREFIX_PATH="$installed_prefix" \
    -DCMAKE_PROJECT_INCLUDE="$provider_file" \
    -DCMAKE_FIND_PACKAGE_PREFER_CONFIG=ON \
    -DSuiteSparse_config_DIR="${installed_prefix}/share/SuiteSparse_config" \
    -DAMD_DIR="${installed_prefix}/share/AMD" \
    -DCOLAMD_DIR="${installed_prefix}/share/COLAMD" \
    -DCAMD_DIR="${installed_prefix}/share/CAMD" \
    -DCCOLAMD_DIR="${installed_prefix}/share/CCOLAMD" \
    -DCMAKE_INSTALL_PREFIX="$cholmod_install" \
    -DCMAKE_SHARED_LINKER_FLAGS="-Wl,--exclude-libs,ALL" \
    -DBUILD_SHARED_LIBS=ON \
    -DBUILD_STATIC_LIBS=OFF \
    -DBUILD_TESTING=OFF \
    -DCHOLMOD_GPL=ON \
    -DCHOLMOD_SUPERNODAL=ON \
    -DCHOLMOD_PARTITION=ON \
    -DCHOLMOD_MATRIXOPS=OFF \
    -DCHOLMOD_MODIFY=OFF \
    -DCHOLMOD_USE_OPENMP=OFF \
    -DCHOLMOD_USE_CUDA=OFF \
    -DSUITESPARSE_USE_CUDA=OFF \
    -DSUITESPARSE_DEMOS=OFF
cmake --build "$cholmod_build" --parallel
cmake --install "$cholmod_build"

project_build="${build_root}/cpu-ipc"
cmake -S "$repo_root" -B "$project_build" -G Ninja \
    -DCMAKE_BUILD_TYPE="$configuration" \
    -DCMAKE_TOOLCHAIN_FILE="$vcpkg_toolchain" \
    -DVCPKG_TARGET_TRIPLET="$triplet" \
    -DCIPC_BUILD_VIEWER="$build_viewer" \
    -DCIPC_CHOLMOD_ROOT="$cholmod_install" \
    -DCIPC_REQUIRE_OPTIMIZED_CHOLMOD=ON
cmake --build "$project_build" --parallel

echo "CPU-IPC ${configuration} build completed: ${project_build}"
echo "Headless executable: ${project_build}/cipc_headless"
