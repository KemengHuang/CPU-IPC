# CPU-IPC

CPU-IPC is a CPU-optimized implementation of Incremental Potential Contact for tetrahedral solids and cloth. It is intended both as a concise and practical CPU simulator and as a reproducible benchmark/reference point when comparing other CPU IPC implementations. The recommended and default Release solver is oneMKL PARDISO with 16 threads.

The benchmark path is fully headless and records solver-stage timing together with numerical/iteration diagnostics. The implementation preserves the same contact, CCD, CFL, energy, gradient, and Hessian semantics across its alternative broad-phase and linear-solver backends, so performance comparisons are not obtained by silently changing the simulated problem.

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
- Current-step Newton convergence checks before CCD/line search, avoiding strict-energy backtracking on numerically vanished directions.
- Strict energy-decreasing line search, Additive CCD, and the original IPC CFL strategy.

## Dependencies

Windows:

```powershell
# Recommended: installs and connects every dependency and builds Release.
powershell -ExecutionPolicy Bypass -File build.ps1
```

WSL/Ubuntu (headless benchmark build by default):

```bash
bash build.sh --headless-only
```

Both launchers locate vcpkg from an explicit option, `VCPKG_ROOT`, or `PATH`. If none exists, they bootstrap an isolated copy automatically—`build/_deps/vcpkg` on Windows and `${XDG_CACHE_HOME:-$HOME/.cache}/cpu-ipc/vcpkg` on Ubuntu. All solver dependencies, including Linux oneMKL, are then installed by vcpkg.

A fresh Ubuntu installation needs a compiler and build frontends once: `sudo apt-get install build-essential cmake ninja-build git curl tar`. These were already present in the WSL installation used for validation. Its missing `zip`/`unzip` tools were downloaded and unpacked into the user cache by `build.sh`, without sudo or system package changes.

## Build

Recommended one-command Windows build—installs dependencies, builds and connects the tuned CHOLMOD DLL, configures CPU-IPC, and compiles Release:

```powershell
powershell -ExecutionPolicy Bypass -File build.ps1
```

Recommended WSL/Ubuntu 22.04 command—uses Linux packages/binaries only and writes them separately under `build-wsl`:

```bash
bash build.sh --headless-only
```

Use `powershell -ExecutionPolicy Bypass -File build.ps1 -HeadlessOnly` on Windows or `bash build.sh --viewer` on WSL to change the application targets. Existing vcpkg checkouts can be selected with `-VcpkgRoot <path>` or `--vcpkg-root <path>`; no repository file contains a machine-specific vcpkg path. The lower-level manual build remains available after one-time setup:

```bash
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build --config Release --parallel
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
| `CIPC_CHOLMOD_ROOT` | auto | Prefix of the tuned supernodal+oneMKL CHOLMOD build. The standard `build/cholmod-mkl-install` is detected automatically. |
| `CIPC_REQUIRE_OPTIMIZED_CHOLMOD` | `ON` | Reject accidental system/OpenBLAS CHOLMOD linkage; set `OFF` only for compatibility diagnostics. |

With METIS ordering enabled, configuration fails early unless the installed CHOLMOD exports the Partition functionality and `cholmod_metis` is linkable. The local SuiteSparse finder accepts `METIS::METIS`, `METIS::metis`, vcpkg's un-namespaced `metis` target, conventional `metis.h + libmetis` installations, and CHOLMOD builds with embedded Partition support. It also selects matching Debug/Release SuiteSparse libraries for multi-configuration builds.

PARDISO is the recommended installation and the default solver in Release-family configurations. It uses oneMKL's static LP64/TBB threading layer and defaults to 16 threads. When PARDISO is unavailable—most notably MSVC Debug—the tuned CHOLMOD backend is selected automatically. Static oneMKL increases the headless executable to about 74.0 MB (70.6 MiB) on the measured setup.

The root launchers build high-performance CHOLMOD automatically: Windows delegates to `scripts/build_cholmod_mkl.ps1`, while Ubuntu performs the equivalent Linux shared-library build inside `build.sh`. Both enable the GPL supernodal module and embed static oneMKL LP64/TBB; CPU-IPC verifies both `cholmod_super_numeric` and `cholmod_metis`. This changes the distribution license boundary to include GPL-2.0-or-later code. A system CHOLMOD is only a compatibility diagnostic when `CIPC_REQUIRE_OPTIMIZED_CHOLMOD=OFF`; it may use OpenBLAS or lack supernodal support.

Headless-only non-quadratic build:

```bash
cmake -S . -B build/nonquadratic \
  -DCMAKE_BUILD_TYPE=Release \
  -DCIPC_BUILD_VIEWER=OFF \
  -DCIPC_ENABLE_QUADRATIC_BENDING=OFF
cmake --build build/nonquadratic --config Release --parallel
```

## Run

Viewer:

```bash
build/Release/cipc
```

Press Space to simulate, `9` to toggle OBJ output, and `/` to toggle screenshots.

The viewer uses the working fixed-function OpenGL path only; the unused shader directory and GLEW dependency have been removed.

Headless:

```bash
build/Release/cipc_headless --scene cloth-bunny --steps 20 --broad-phase lbvh --output Output/run
build/Release/cipc_headless --scene cloth-bunny --steps 20 --linear-solver cholmod --no-output
build/Release/cipc_headless --scene twisting-mat --steps 1 --linear-solver eigen-cg --no-output
build/Release/cipc_headless --scene bunny2 --steps 1 --linear-solver pardiso --pardiso-threads 16 --no-output
python scripts/benchmark.py --exe build/Release/cipc_headless.exe --repeats 5 --steps 20
```

WSL/Ubuntu headless executable:

```bash
./build-wsl/cpu-ipc/cipc_headless --scene bunny2 --steps 1 --no-output
```

Runtime metrics are written to `metrics.csv` in the selected output directory. See [`agent_docs/README.md`](agent_docs/README.md) for architecture, algorithms, known issues, and optimization results, and [`agent_docs/12_wsl_ubuntu_build.md`](agent_docs/12_wsl_ubuntu_build.md) for the verified WSL/Ubuntu toolchain and build flow.

Scene construction is fresh by default. Checkpoint loading and writing are opt-in through `--resume` and `--write-checkpoints`.

The CPU LBVH broad phase is the default; use `--broad-phase spatial-hash` for the optimized legacy backend.

PARDISO is the project's primary and default Newton solver. It uses 16 threads unless `--pardiso-threads` overrides the limit. Tuned CHOLMOD is the automatic fallback and can be selected explicitly with `--linear-solver cholmod`; `--cholmod-threads 0` uses the measured automatic policy (4 threads below 500k matrix nonzeros, otherwise 8), while a positive value overrides it. Eigen-CG remains an additional comparison backend. All three share the same lower-triangular Hessian assembly and boundary handling.

Each Newton convergence decision uses the direction just solved from the current gradient/Hessian. A converged direction exits before CCD and line search; the strict acceptance rule remains `E_trial < E0` with `armijoCoefficient=0`. This prevents a numerically vanished direction from being halved repeatedly only because parallel energy summation fluctuates in the last bit. On the two-bunny2 scene, five independent five-step runs reduced the first 25 frames from 417 energy backtracks to 0, made all five final states identical, and reduced median total time from 4.429 s to 3.396 s. The final convergence-check factorization is still counted, so `numeric_factorizations` can be one greater than the number of accepted Newton updates.

The tuned CHOLMOD path is about 5.3× faster than the supernodal+OpenBLAS build on bunny2. With the current convergence semantics, five-run medians are: bunny2 one step 1.625 s for CHOLMOD versus 1.301 s for PARDISO; cloth-bunny five steps 0.608 s versus 0.586 s; twisting-mat five steps 0.413 s versus 0.434 s. The two production solvers are therefore in the same performance tier: generally within about 5% on the medium scenes, with PARDISO retaining about a 25% advantage on the large bunny2 scene. PARDISO remains the default, while tuned CHOLMOD is the automatic backup. Numerical results agree within floating-point reduction error. See [`agent_docs/11_cholmod_mkl_report.md`](agent_docs/11_cholmod_mkl_report.md) for the full ordering/thread/provider A/B and [`agent_docs/10_pardiso_report.md`](agent_docs/10_pardiso_report.md) for PARDISO details.

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

## Bending models

Both bending paths use the thin-plate rigidity

```text
D = E * thickness^3 / (12 * (1 - poisson_ratio^2)).
```

- Quadratic bending precomputes the constant cotangent/area matrix `Q`.
- Non-quadratic bending precomputes each interior hinge's rest angle and geometry weight `l0 / (h0 + h1)`, then uses consistent energy, gradient, and Hessian expressions.
- `TetInversionGuard` is retained as an optional utility but is not enabled in the default IPC step.

## Project layout

| Path | Responsibility |
|---|---|
| `CPU IPC/IPCSolver.*` | Time stepping, current-direction Newton convergence, CCD/CFL integration, and strict energy line search. |
| `CPU IPC/ContactMechanics.*` | Contact distances, barriers, derivatives, and feasible self-contact steps. |
| `CPU IPC/CollisionBroadPhase.*`, `LBVH.*` | SpatialHash and Linear BVH broad phases. |
| `CPU IPC/Elasticity.*`, `HingeBending.*` | Solid/cloth elasticity and bending models. |
| `CPU IPC/NewtonLinearSystem.*`, `CholmodSolver.*`, `PardisoSolver.*` | Shared sparse assembly and optimized CHOLMOD/PARDISO/Eigen-CG backends. |
| `CPU IPC/SimulationMesh.*`, `Simulator.*` | Mesh state, IO, scene construction, and simulation ownership. |
| `CPU IPC/ViewerMain.cpp` | Fixed-function GLUT viewer. |
| `apps/cipc_headless.cpp` | Headless executable. |
| `build.ps1`, `build.sh` | Portable one-command Windows and WSL/Ubuntu production builds. |
| `scripts/` | Benchmark and video utilities. |

## Manual smoke checks

The repository intentionally ships no test suite. After changes, validate the product executables directly:

```bash
build/Release/cipc_headless --scene cloth-bunny --steps 1 --no-output --broad-phase lbvh
build/Release/cipc_headless --scene cloth-bunny --steps 1 --no-output --broad-phase spatial-hash
build/Release/cipc_headless --scene cloth-bunny --steps 1 --no-output --linear-solver cholmod
build/Release/cipc_headless --scene twisting-mat --steps 1 --no-output --linear-solver eigen-cg
build/Release/cipc_headless --scene bunny2 --steps 1 --no-output --linear-solver pardiso --pardiso-threads 16
build/nonquadratic/Release/cipc_headless --scene cloth-bunny --steps 1 --no-output
```
