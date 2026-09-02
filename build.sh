#!/usr/bin/env bash
set -euo pipefail

repo_root="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
build_root="${repo_root}/build-wsl"
configuration="Release"
build_viewer="OFF"
requested_vcpkg_root="${VCPKG_ROOT:-}"
cache_root="${XDG_CACHE_HOME:-${HOME}/.cache}/cpu-ipc"
install_system_packages="ON"
dependencies_only="OFF"
vcpkg_revision="$(tr -d '[:space:]' < "${repo_root}/scripts/vcpkg-revision.txt")"
suitesparse_version="$(tr -d '[:space:]' < "${repo_root}/scripts/suitesparse-version.txt")"
if [[ ! "$vcpkg_revision" =~ ^[0-9a-fA-F]{40}$ ]]; then
    echo "scripts/vcpkg-revision.txt must contain one 40-character Git revision." >&2
    exit 1
fi
if [[ ! "$suitesparse_version" =~ ^[0-9]+\.[0-9]+\.[0-9]+$ ]]; then
    echo "scripts/suitesparse-version.txt must contain a semantic version." >&2
    exit 1
fi

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
        --dependencies-only)
            dependencies_only="ON"
            shift
            ;;
        --no-system-packages)
            install_system_packages="OFF"
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
  --dependencies-only
                     Install/build dependencies without building CPU-IPC.
  --no-system-packages
                     Never invoke sudo/apt; report missing prerequisites.
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

. "${repo_root}/scripts/ensure_ubuntu_prerequisites.sh"
echo "[1/4] Checking the Ubuntu build environment"
ensure_ubuntu_prerequisites "$cache_root" "$build_viewer" "$install_system_packages"

echo "[2/4] Preparing vcpkg and numerical dependencies"
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
else
    vcpkg_root="${cache_root}/vcpkg"
    mkdir -p "$cache_root"
    needs_bootstrap="OFF"
    if [[ ! -d "${vcpkg_root}/.git" ]]; then
        mkdir -p "$vcpkg_root"
        git -C "$vcpkg_root" init
        git -C "$vcpkg_root" remote add origin https://github.com/microsoft/vcpkg.git
        needs_bootstrap="ON"
    fi
    current_revision="$(git -C "$vcpkg_root" rev-parse HEAD 2>/dev/null || true)"
    if [[ "$current_revision" != "$vcpkg_revision" ]]; then
        git -C "$vcpkg_root" fetch --depth 1 origin "$vcpkg_revision"
        git -C "$vcpkg_root" checkout --detach FETCH_HEAD
        needs_bootstrap="ON"
    fi
    if [[ "$needs_bootstrap" == "ON" || ! -x "${vcpkg_root}/vcpkg" ]]; then
        if [[ ! -x "${vcpkg_root}/bootstrap-vcpkg.sh" ]]; then
            echo "The pinned vcpkg checkout is incomplete: $vcpkg_root" >&2
            exit 1
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
    "tbb[core]:${triplet}"
    "metis:${triplet}"
    "intel-mkl:${triplet}"
)
if [[ "$build_viewer" == "ON" ]]; then
    dependencies+=("freeglut:${triplet}")
fi
"$vcpkg_executable" install "${dependencies[@]}" --recurse --no-print-usage
cmake_executable="$(resolve_cmake_executable "$vcpkg_root")"

installed_prefix="${vcpkg_root}/installed/${triplet}"
suitesparse_source_parent="${build_root}/_deps"
suitesparse_source="${suitesparse_source_parent}/SuiteSparse-${suitesparse_version}"
download_directory="${vcpkg_root}/downloads"
fetch_script="${repo_root}/cmake/FetchSuiteSparse.cmake"
cholmod_build="${build_root}/suitesparse-mkl-bundle"
cholmod_install="${build_root}/cholmod-mkl-install"
provider_file="${repo_root}/cmake/CholmodMKLProvider.cmake"
blas_shim="${repo_root}/cmake/SuiteSparseBLAS.cmake"
bundle_stamp="${cholmod_install}/.cpu-ipc-suite-sparse-bundle"
# Increment recipe when the SuiteSparse configure arguments below change.
bundle_signature="$(printf '%s\n' \
    'schema=1' \
    'recipe=1' \
    "suitesparse=${suitesparse_version}" \
    "vcpkg_revision=${vcpkg_revision}" \
    "vcpkg=${vcpkg_root}" \
    "vcpkg_executable=$(sha256sum "$vcpkg_executable" | cut -d' ' -f1)" \
    "mkl_port=$(sha256sum "${vcpkg_root}/ports/intel-mkl/vcpkg.json" | cut -d' ' -f1)" \
    "tbb_port=$(sha256sum "${vcpkg_root}/ports/tbb/vcpkg.json" | cut -d' ' -f1)" \
    "metis_port=$(sha256sum "${vcpkg_root}/ports/metis/vcpkg.json" | cut -d' ' -f1)" \
    "provider=$(sha256sum "$provider_file" | cut -d' ' -f1)" \
    "blas_shim=$(sha256sum "$blas_shim" | cut -d' ' -f1)" \
    "fetch=$(sha256sum "$fetch_script" | cut -d' ' -f1)")"

bundle_ready="ON"
required_bundle_files=(
    "lib/libsuitesparseconfig.so"
    "lib/libamd.so"
    "lib/libcamd.so"
    "lib/libccolamd.so"
    "lib/libcolamd.so"
    "lib/libcholmod.so"
)
if [[ ! -f "$bundle_stamp" || "$(cat "$bundle_stamp" 2>/dev/null || true)" != "$bundle_signature" ]]; then
    bundle_ready="OFF"
fi
for relative_path in "${required_bundle_files[@]}"; do
    [[ -f "${cholmod_install}/${relative_path}" ]] || bundle_ready="OFF"
done

if [[ "$bundle_ready" == "ON" ]]; then
    echo "[3/4] Reusing the optimized SuiteSparse/CHOLMOD bundle"
else
    "$cmake_executable" \
        "-DCPU_IPC_SUITESPARSE_SOURCE_PARENT=${suitesparse_source_parent}" \
        "-DCPU_IPC_DOWNLOAD_DIR=${download_directory}" \
        -P "$fetch_script"

    echo "[3/4] Building the optimized SuiteSparse/CHOLMOD bundle with oneMKL"
    "$cmake_executable" -S "$suitesparse_source" -B "$cholmod_build" -G Ninja \
        -DCMAKE_BUILD_TYPE=Release \
        -DCMAKE_TOOLCHAIN_FILE="$vcpkg_toolchain" \
        -DVCPKG_TARGET_TRIPLET="$triplet" \
        -DCMAKE_PREFIX_PATH="$installed_prefix" \
        -DCMAKE_PROJECT_INCLUDE="$provider_file" \
        -DCMAKE_FIND_PACKAGE_PREFER_CONFIG=ON \
        -DCMAKE_INSTALL_PREFIX="$cholmod_install" \
        -DCMAKE_SHARED_LINKER_FLAGS="-Wl,--exclude-libs,ALL" \
        -DSUITESPARSE_ENABLE_PROJECTS=cholmod \
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
        -DSUITESPARSE_USE_OPENMP=OFF \
        -DSUITESPARSE_USE_FORTRAN=OFF \
        -DSUITESPARSE_DEMOS=OFF \
        -DSUITESPARSE_USE_STRICT=ON
    "$cmake_executable" --build "$cholmod_build" --parallel
    "$cmake_executable" --install "$cholmod_build"
    printf '%s\n' "$bundle_signature" > "$bundle_stamp"
fi

if [[ "$dependencies_only" == "ON" ]]; then
    echo "CPU-IPC dependencies are ready."
    echo "vcpkg: ${vcpkg_root}"
    echo "Optimized CHOLMOD: ${cholmod_install}"
    exit 0
fi

echo "[4/4] Configuring and building CPU-IPC"
project_build="${build_root}/cpu-ipc"
"$cmake_executable" -S "$repo_root" -B "$project_build" -G Ninja \
    -DCMAKE_BUILD_TYPE="$configuration" \
    -DCMAKE_TOOLCHAIN_FILE="$vcpkg_toolchain" \
    -DVCPKG_TARGET_TRIPLET="$triplet" \
    -DCIPC_BUILD_VIEWER="$build_viewer" \
    -DCIPC_ENABLE_FRICTION=ON \
    -DCIPC_ENABLE_QUADRATIC_BENDING=ON \
    -DCIPC_ENABLE_METIS_ORDERING=ON \
    -DCIPC_ENABLE_PARDISO=ON \
    -DCIPC_CHOLMOD_ROOT="$cholmod_install" \
    -DCIPC_REQUIRE_OPTIMIZED_CHOLMOD=ON
"$cmake_executable" --build "$project_build" --parallel

echo "CPU-IPC ${configuration} build completed: ${project_build}"
echo "Headless executable: ${project_build}/cipc_headless"
