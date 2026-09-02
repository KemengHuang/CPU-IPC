# CPU-IPC

CPU-IPC is a CPU-optimized implementation of Incremental Potential Contact for tetrahedral solids and cloth. It is intended both as a concise and practical CPU simulator and as an optional reproducible benchmark/reference point when comparing other CPU IPC implementations. The recommended and default Release solver is oneMKL PARDISO with 16 threads.

The benchmark path is fully headless and records solver-stage timing together with numerical/iteration diagnostics. The implementation preserves the same contact, CCD, CFL, energy, gradient, and Hessian semantics across its alternative broad-phase and linear-solver backends, so performance comparisons are not obtained by silently changing the simulated problem.

## Quick start

No manual vcpkg, Eigen, TBB, METIS, oneMKL, SuiteSparse, CHOLMOD, or FreeGLUT setup is required.

Windows, from Command Prompt or PowerShell:

```bat
.\build.cmd
```

WSL/Ubuntu:

```bash
./build.sh
```

The first run downloads the numerical dependencies and builds the tuned CHOLMOD library, so it is substantially slower than later incremental builds. Both platforms produce the viewer and headless runner by default; use `.\build.cmd -HeadlessOnly` or `./build.sh --headless-only` on machines without a graphical environment.

### CPU-oriented optimizations

- TBB-parallel elastic/contact assembly, energy evaluation, broad-phase construction, and reductions.
- CPU Morton/Karras Linear BVH as the default broad phase, with an optimized SpatialHash backend for A/B comparisons.
- Backend-independent exact AABB filtering before CCD so LBVH and SpatialHash use equivalent effective contact pairs.
- Fixed-size Eigen element kernels, precomputed PFPX operators, and precomputed quadratic/hinge bending geometry.
- Lower-triangular sparse Hessian assembly without placeholder zeros.
- Reused Newton, energy, sparse-matrix, RHS, and solver workspaces.
- Default oneMKL PARDISO SPD solver with parallel factorization, cross-time-step symbolic reuse, adaptive METIS-permutation refresh, and per-phase metrics.
- Optional high-performance CHOLMOD comparison path: GPL supernodal factorization, embedded oneMKL LP64/TBB, measured AMD ordering, and adaptive 4/8-thread limits.
- Optional Eigen-CG backend using the same assembled system and boundary handling.
- Centralized boundary handling with validated hard Dirichlet motion and soft target penalties whose energy, gradient, and Hessian are assembled consistently.
- Current-step Newton convergence checks before CCD/line search, avoiding strict-energy backtracking on numerically vanished directions.
- Strict energy-decreasing line search, Additive CCD, and the original IPC CFL strategy.

## System requirements

- Windows: 64-bit Windows and Visual Studio/Build Tools with the **Desktop development with C++** workload. `build.cmd` detects this prerequisite and prints a ready-to-run `winget` command if it is absent. Git is needed only when the private vcpkg cache must first be downloaded or updated, with an equivalent command shown on failure. CMake is resolved from `PATH`, Visual Studio, or vcpkg's own downloaded tools.
- WSL/Ubuntu: x86-64 Ubuntu. `build.sh` detects missing compiler/build/OpenGL packages and installs them through `apt` automatically, requesting sudo only when necessary. `zip`/`unzip` can be unpacked into the user cache without sudo, and a system CMake older than 3.23 is transparently replaced by vcpkg's private tool. The viewer runs with WSLg or another working X11 `DISPLAY`; use `--headless-only` on servers. Pass `--no-system-packages` in CI or managed environments to disable automatic `apt` changes and receive the exact required command instead.

Everything else is project-managed. An explicit vcpkg option or `VCPKG_ROOT` overrides the default; otherwise the launchers bootstrap the revision pinned in `scripts/vcpkg-revision.txt`. Windows keeps that checkout in a short `%LOCALAPPDATA%/CPU-IPC` path to avoid path-length failures, while Ubuntu uses `${XDG_CACHE_HOME:-$HOME/.cache}/cpu-ipc`. Eigen, TBB core, METIS, oneMKL, and optional FreeGLUT are installed together. SuiteSparse is downloaded at a pinned version with a verified SHA-512 hash and only the required CHOLMOD/AMD family is built—OpenBLAS and TBB's unused hwloc feature are not installed.

## Build

Recommended Windows build—installs dependencies, builds and connects the tuned CHOLMOD DLL, configures CPU-IPC, and compiles Release:

```bat
.\build.cmd
```

Recommended WSL/Ubuntu 22.04 command—uses Linux packages/binaries only and writes them separately under `build-wsl`:

```bash
./build.sh
```

Use `.\build.cmd -DependenciesOnly` or `./build.sh --dependencies-only` when only preparing libraries; `.\build.cmd -Help` and `./build.sh --help` list all options. Existing vcpkg checkouts can be selected with `-VcpkgRoot <path>` or `--vcpkg-root <path>`; no repository file contains a machine-specific vcpkg path. The PowerShell entry remains available as `powershell -ExecutionPolicy Bypass -File build.ps1`. After the one-command setup, incremental compilation is simply:

```powershell
cmake --build build/cpu-ipc --config Release --parallel
```

```bash
cmake --build build-wsl/cpu-ipc --parallel
```

For best WSL compilation throughput, keeping the repository in the Linux filesystem (for example under `~/src`) is preferable, but the complete workflow was also verified directly from this repository through `/mnt/c`.

Targets:

- `cipc_core`: simulation library without OpenGL.
- `cipc`: GLUT/OpenGL viewer. Disable with `-DCIPC_BUILD_VIEWER=OFF`.
- `cipc_headless`: command-line simulation and benchmarking runner.

The build intentionally contains no test targets or CTest registration.

### Configuration options

| Option | Default | Purpose |
|---|---:|---|
| `CIPC_BUILD_VIEWER` | `ON` | Build the GLUT viewer; set `OFF` for a headless-only build. |
| `CIPC_ENABLE_FRICTION` | `ON` | Enable lagged IPC friction. |
| `CIPC_ENABLE_QUADRATIC_BENDING` | `ON` | Use quadratic isometric bending; set `OFF` for the complete dihedral-hinge model. |
| `CIPC_ENABLE_METIS_ORDERING` | `ON` | Build/link CHOLMOD Partition/METIS support for experiments; measured production ordering remains AMD. |
| `CIPC_ENABLE_PARDISO` | `ON` | Build the recommended/default oneMKL PARDISO backend; when disabled, optimized CHOLMOD becomes the default. |
| `CIPC_CHOLMOD_ROOT` | auto | Prefix of the tuned supernodal+oneMKL SuiteSparse bundle. The standard `build/cholmod-mkl-install` is detected automatically. |
| `CIPC_REQUIRE_OPTIMIZED_CHOLMOD` | `ON` | Reject accidental system/OpenBLAS CHOLMOD linkage; set `OFF` only for compatibility diagnostics. |

With METIS ordering enabled, configuration fails early unless the installed CHOLMOD exports the Partition functionality and `cholmod_metis` is linkable. The local SuiteSparse finder accepts `METIS::METIS`, `METIS::metis`, vcpkg's un-namespaced `metis` target, conventional `metis.h + libmetis` installations, and CHOLMOD builds with embedded Partition support. It also selects matching Debug/Release SuiteSparse libraries for multi-configuration builds.

PARDISO is the recommended installation and the default solver in Release-family configurations. It uses oneMKL's static LP64/TBB threading layer and defaults to 16 threads. When PARDISO is unavailable—most notably MSVC Debug—the tuned CHOLMOD backend is selected automatically. Static oneMKL increases the headless executable to about 74.0 MB (70.6 MiB) on the measured setup.

The root launchers build high-performance CHOLMOD automatically: Windows delegates to `scripts/build_cholmod_mkl.ps1`, while Ubuntu performs the equivalent build inside `build.sh`. Both build a self-contained shared SuiteSparse bundle containing only Config, AMD, CAMD, CCOLAMD, COLAMD, and GPL supernodal CHOLMOD, with static oneMKL LP64/TBB as the sole BLAS/LAPACK provider. CPU-IPC verifies `cholmod_super_numeric` and `cholmod_metis`, and deploys every required DLL on Windows. A content stamp skips this bundle entirely on unchanged incremental builds. This changes the distribution license boundary to include GPL-2.0-or-later code. A system CHOLMOD is only a compatibility diagnostic when `CIPC_REQUIRE_OPTIMIZED_CHOLMOD=OFF`.

Cross-platform non-quadratic build:

```powershell
.\build.cmd -HeadlessOnly -NonQuadraticBending -BuildDirectory build/cpu-ipc-nonquadratic
```

```bash
./build.sh --headless-only --nonquadratic-bending \
  --project-build-dir build-wsl/nonquadratic
```

## Run

Viewer:

```bash
build/cpu-ipc/Release/cipc.exe
build/cpu-ipc/Release/cipc.exe --scene twisting-mat-soft
```

Without `--scene`, the viewer runs the `twisting-mat-soft` soft-boundary example. The headless executable, `SimulationOptions`, `FEMSimulator::buildModels()`, and `scripts/benchmark.py` use the same default.

Press Space to simulate, `9` to toggle OBJ output, and `/` to toggle screenshots.

The viewer uses the working fixed-function OpenGL path only; the unused shader directory and GLEW dependency have been removed.

Headless:

```bash
build/cpu-ipc/Release/cipc_headless.exe --scene cloth-bunny --steps 20 --broad-phase lbvh --output Output/run
build/cpu-ipc/Release/cipc_headless.exe --scene cloth-bunny --steps 20 --linear-solver cholmod --no-output
build/cpu-ipc/Release/cipc_headless.exe --scene twisting-mat --steps 1 --linear-solver eigen-cg --no-output
build/cpu-ipc/Release/cipc_headless.exe --scene twisting-mat-soft --steps 20 --no-output
build/cpu-ipc/Release/cipc_headless.exe --scene bunny2 --steps 1 --linear-solver pardiso --pardiso-threads 16 --no-output
python scripts/benchmark.py --exe build/cpu-ipc/Release/cipc_headless.exe --repeats 5 --steps 20
```

WSL/Ubuntu headless executable:

```bash
./build-wsl/cpu-ipc/cipc_headless --scene bunny2 --steps 1 --no-output
```

WSL/Ubuntu viewer:

```bash
./build-wsl/cpu-ipc/cipc
./build-wsl/cpu-ipc/cipc --scene twisting-mat-soft
```

Runtime metrics are written to `metrics.csv` in the selected output directory. See [`agent_docs/README.md`](agent_docs/README.md) for architecture, algorithms, known issues, and optimization results, and [`agent_docs/12_wsl_ubuntu_build.md`](agent_docs/12_wsl_ubuntu_build.md) for the verified WSL/Ubuntu toolchain and build flow.

Scene construction is fresh by default. Checkpoint loading and writing are opt-in through `--resume` and `--write-checkpoints`.

The CPU LBVH broad phase is the default; use `--broad-phase spatial-hash` for the optimized legacy backend.

PARDISO is the project's primary and default Newton solver. It uses 16 threads unless `--pardiso-threads` overrides the limit. Tuned CHOLMOD is the automatic fallback and can be selected explicitly with `--linear-solver cholmod`; `--cholmod-threads 0` uses the measured automatic policy (4 threads below 500k matrix nonzeros, otherwise 8), while a positive value overrides it. Eigen-CG remains an additional comparison backend. All three share the same lower-triangular Hessian assembly and boundary handling.

Each Newton convergence decision uses the direction just solved from the current gradient/Hessian. A converged direction exits before CCD and line search; the strict acceptance rule remains `E_trial < E0`. This prevents a numerically vanished direction from being halved repeatedly only because parallel energy summation fluctuates in the last bit.

CPU-IPC uses a specially tuned CHOLMOD backend that is substantially faster than the generic, out-of-the-box CHOLMOD builds installed through Ubuntu's `apt install libsuitesparse-dev` or through `vcpkg install suitesparse` on Linux and Windows. It performs in the same general performance tier as PARDISO; PARDISO remains faster on the largest tested scene and is therefore the default solver. Tuned CHOLMOD serves as the automatic fallback, and both backends produce numerical results that agree within floating-point reduction error. See [11_cholmod_mkl_report.md](agent_docs/11_cholmod_mkl_report.md) for the full ordering/thread/provider A/B comparison and [10_pardiso_report.md](agent_docs/10_pardiso_report.md) for PARDISO details.

The one-command production build opts into GPL supernodal CHOLMOD; treat its binaries as GPL-covered distributions. Setting `CIPC_REQUIRE_OPTIMIZED_CHOLMOD=OFF` is only for system-package compatibility diagnostics.

## CPU IPC benchmark usage

Use a Release, headless-only build when comparing CPU implementations:

```bash
cmake -S . -B build/benchmark \
  -DCMAKE_BUILD_TYPE=Release \
  -DCIPC_BUILD_VIEWER=OFF
cmake --build build/benchmark --config Release --parallel

python scripts/benchmark.py \
  --exe build/benchmark/Release/cipc_headless.exe \
  --scene bunny2 \
  --broad-phase lbvh \
  --linear-solver pardiso \
  --pardiso-threads 16 \
  --steps 1 \
  --repeats 3
```

The benchmark launches independent processes, stores each run under `Output/benchmark/run_NN`, sums per-frame `step_ms`, and reports median/min/max. Its default `--linear-solver auto` follows the same PARDISO→optimized-CHOLMOD selection; pass a concrete backend for controlled A/B runs. `metrics.csv` separates assembly, linear solve, CCD, line search, and post-line-search time while recording Newton iterations, backtracks, accepted step sizes, contact counts, matrix nonzeros, direct-solver symbolic/numeric counts, and effective threads.

For fair comparisons:

- use the same scene, material file, timestep, bending mode, broad-phase semantics, linear-solver backend, tolerance, METIS-ordering setting, PARDISO thread count, and number of time steps;
- disable checkpoint resume and avoid viewer/rendering work;
- compare `RESULT` and physical/iteration metrics in addition to wall time;
- report the median of multiple independent runs, not only the fastest sample;
- do not compare absolute timings from different machines as if they were directly interchangeable.

## Boundary conditions

`BoundaryConditionSet` keeps hard and soft boundary data separate. Animated Dirichlet vertices are marked constrained in `boundaryTypes`; their `updateDirection(current, rest, step, alpha, dt)` callback returns the search direction used by the existing `x_new = x - p` convention, and the motion remains guarded by CCD. Static constrained vertices need no callback.

### Hard Dirichlet boundary setup

For an animated hard boundary, mark every selected vertex as Dirichlet, add it to the animated Dirichlet list, and provide the update callback:

```cpp
for (int vertex : selectedVertices) {
    mesh.boundaryTypes[vertex] =
        boundaryTypeCode(VertexBoundaryType::Dirichlet);
    mesh.boundaryConditions.dirichlet.vertexIndices.push_back(vertex);
}

mesh.boundaryConditions.dirichlet.updateDirection =
    [](const Eigen::Vector3d& current,
       const Eigen::Vector3d& rest,
       int step,
       double alpha,
       double dt) {
        // Example: prescribe an upward velocity of 0.1 units/second.
        // The callback returns p, while the solver applies x_new = x - p.
        const Eigen::Vector3d target =
            current + alpha * dt * Eigen::Vector3d(0.0, 0.1, 0.0);
        return current - target;
    };
```

`alpha` is the fractional hard-boundary step selected by CCD. Use `rest` to classify vertices without that classification changing as the mesh deforms, and use `step` for staged or time-dependent motion. For a static hard boundary, only set `boundaryTypes[vertex]` to `Dirichlet`; do not add it to `dirichlet.vertexIndices` and do not install an update callback.

### Soft target boundary setup

Soft-boundary vertices must remain free Newton unknowns. Add them only to the soft list, choose a finite positive weight, and return an absolute target position from `updateTarget`:

```cpp
for (int vertex : selectedVertices) {
    mesh.boundaryTypes[vertex] =
        boundaryTypeCode(VertexBoundaryType::Free);
    mesh.boundaryConditions.soft.vertexIndices.push_back(vertex);
}

mesh.boundaryConditions.soft.weight = 100.0;
mesh.boundaryConditions.soft.updateTarget =
    [](const Eigen::Vector3d& current,
       const Eigen::Vector3d& rest,
       int step,
       double dt) {
        // Targets are evaluated once and then frozen for the whole time step.
        return current + dt * Eigen::Vector3d(0.0, 0.1, 0.0);
    };
```

For a static soft target, omit `updateTarget` and fill `soft.targetPositions` in exactly the same order as `soft.vertexIndices`. If `targetPositions` is omitted, initialization uses the selected vertices' current positions.

The historical `update_soft_constraint_functor` path is represented by `soft.updateTarget(current, rest, step, dt)`. Its vertices remain free Newton unknowns and contribute the matched incremental-potential terms

```text
E_soft = 0.5 * weight * ||x - target||^2
g_soft = weight * (x - target)
H_soft = weight * I.
```

This is a target-position penalty (often called the project's soft or Neumann-style boundary), mathematically closer to a spring/Robin condition than a pure prescribed-traction Neumann condition. Targets are frozen once per time step so line search evaluates one consistent objective. Configuration validation rejects duplicate/out-of-range indices, non-finite targets and weights, and soft vertices accidentally marked as hard-constrained.

Configure the boundary sets after loading/selecting mesh vertices. The normal `FEMSimulator::buildModels()` path assigns `v_rest` and calls `BoundaryConditionOps::initialize(mesh)` automatically. Code that constructs a `mesh3D` outside that path must do so once after finalizing the rest positions:

```cpp
mesh.v_rest = mesh.vertexes;
BoundaryConditionOps::initialize(mesh);
```

Do not place the same vertex in both sets. `BoundaryConditionOps::initialize()` validates array sizes, uniqueness, index ranges, hard/free status, targets, and weights before simulation starts.

`--scene twisting-mat` is the original hard-Dirichlet example. `--scene twisting-mat-soft` uses the same mesh, selected end vertices, and per-step target endpoint, but keeps those vertices free and attracts them with `weight=100`; it is the runnable soft-boundary example for both the viewer and headless benchmark.

## Project layout

| Path | Responsibility |
|---|---|
| `CPU IPC/IPCSolver.*` | Time stepping, current-direction Newton convergence, CCD/CFL integration, and strict energy line search. |
| `CPU IPC/ContactMechanics.*` | Contact distances, barriers, derivatives, and feasible self-contact steps. |
| `CPU IPC/CollisionBroadPhase.*`, `LBVH.*` | SpatialHash and Linear BVH broad phases. |
| `CPU IPC/Elasticity.*`, `HingeBending.*` | Solid/cloth elasticity and bending models. |
| `CPU IPC/BoundaryConditions.*` | Dirichlet motion, soft-target updates, validation, and soft energy/derivatives. |
| `CPU IPC/NewtonLinearSystem.*`, `CholmodSolver.*`, `PardisoSolver.*` | Shared sparse assembly and optimized CHOLMOD/PARDISO/Eigen-CG backends. |
| `CPU IPC/SimulationMesh.*`, `Simulator.*` | Mesh state, IO, scene construction, and simulation ownership. |
| `CPU IPC/ViewerMain.cpp` | Fixed-function GLUT viewer. |
| `apps/cipc_headless.cpp` | Headless executable. |
| `build.cmd`, `build.ps1`, `build.sh` | One-command dependency setup and Windows/WSL/Ubuntu production builds. |
| `scripts/` | Benchmark and video utilities. |

## Manual smoke checks

The repository intentionally ships no test suite. After changes, validate the product executables directly:

```bash
build/cpu-ipc/Release/cipc_headless.exe --scene cloth-bunny --steps 1 --no-output --broad-phase lbvh
build/cpu-ipc/Release/cipc_headless.exe --scene cloth-bunny --steps 1 --no-output --broad-phase spatial-hash
build/cpu-ipc/Release/cipc_headless.exe --scene cloth-bunny --steps 1 --no-output --linear-solver cholmod
build/cpu-ipc/Release/cipc_headless.exe --scene twisting-mat --steps 1 --no-output --linear-solver eigen-cg
build/cpu-ipc/Release/cipc_headless.exe --scene twisting-mat-soft --steps 5 --no-output --linear-solver pardiso
build/cpu-ipc/Release/cipc_headless.exe --scene bunny2 --steps 1 --no-output --linear-solver pardiso --pardiso-threads 16
build/cpu-ipc-nonquadratic/Release/cipc_headless.exe --scene cloth-bunny --steps 1 --no-output
```
