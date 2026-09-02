# CPU-IPC

CPU-IPC is a CPU-optimized implementation of Incremental Potential Contact for tetrahedral solids and cloth. It is intended both as a concise and practical CPU simulator and as an optional reproducible benchmark/reference point when comparing other CPU IPC implementations.

The benchmark path is fully headless and records solver-stage timing together with numerical/iteration diagnostics. The implementation preserves the same contact, CCD, CFL, energy, gradient, and Hessian semantics across its alternative broad-phase and linear-solver backends, so performance comparisons are not obtained by silently changing the simulated problem.

### CPU-oriented optimizations

- TBB-parallel elastic/contact assembly, energy evaluation, broad-phase construction, and reductions.
- CPU Morton/Karras Linear BVH as the default broad phase, with an optimized SpatialHash backend for A/B comparisons.
- Backend-independent exact AABB filtering before CCD so LBVH and SpatialHash use equivalent effective contact pairs.
- Fixed-size Eigen element kernels, precomputed PFPX operators, and precomputed quadratic/hinge bending geometry.
- Lower-triangular sparse Hessian assembly without placeholder zeros.
- Reused Newton, energy, sparse-matrix, RHS, and solver workspaces.
- Block-aware SuiteSparse LDL with reusable symbolic data and a pre-permuted upper CSC numeric path; this is the automatic fallback when PARDISO is unavailable.
- Optional CHOLMOD with symbolic-factorization reuse and a verified, cost-aware METIS nested-dissection fallback.
- Optional oneMKL PARDISO SPD backend with parallel factorization, cross-time-step symbolic reuse, adaptive METIS-permutation refresh, and per-phase metrics.
- Optional Eigen-CG backend using the same assembled system and boundary handling.
- Current-step Newton convergence checks before CCD/line search, avoiding strict-energy backtracking on numerically vanished directions.
- Strict energy-decreasing line search, Additive CCD, and the original IPC CFL strategy.

## Dependencies
Windows:
```bash
vcpkg install eigen3 freeglut tbb openblas suitesparse metis
# Optional parallel PARDISO backend:
vcpkg install intel-mkl:x64-windows
```

Ubuntu:
```bash
sudo apt install libeigen3-dev freeglut3-dev libtbb-dev libopenblas-dev libsuitesparse-dev libmetis-dev
```

## Build

```bash
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build --config Release --parallel
```

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
| `CIPC_ENABLE_METIS_ORDERING` | `ON` | Require CHOLMOD Partition/METIS support and allow CHOLMOD's fill/work heuristic to select METIS after AMD; set `OFF` for AMD-only ordering. |
| `CIPC_ENABLE_PARDISO` | `ON` | Build the oneMKL PARDISO backend when an MKL CMake package is available; otherwise keep the other solvers available. |

With METIS ordering enabled, configuration fails early unless the installed CHOLMOD exports the Partition functionality and `cholmod_metis` is linkable. The local SuiteSparse finder accepts `METIS::METIS`, `METIS::metis`, vcpkg's un-namespaced `metis` target, conventional `metis.h + libmetis` installations, and CHOLMOD builds with embedded Partition support. It also selects matching Debug/Release SuiteSparse libraries for multi-configuration builds.

PARDISO uses oneMKL's static LP64/TBB threading layer in Release-family configurations. If oneMKL is absent, configuration succeeds but selecting `--linear-solver pardiso` reports that the backend is unavailable. The Windows oneMKL static package uses the Release CRT, so MSVC Debug builds deliberately exclude only PARDISO while retaining the rest of the project; use Release or RelWithDebInfo for PARDISO benchmarking. Static oneMKL also increases the headless executable to about 74.0 MB (70.6 MiB) on the measured setup.

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

Runtime metrics are written to `metrics.csv` in the selected output directory. See [`agent_docs/README.md`](agent_docs/README.md) for architecture, algorithms, known issues, and optimization results.

Scene construction is fresh by default. Checkpoint loading and writing are opt-in through `--resume` and `--write-checkpoints`.

The CPU LBVH broad phase is the default; use `--broad-phase spatial-hash` for the optimized legacy backend.

The default Newton solver is selected from build capabilities: PARDISO when oneMKL is available in the current configuration, otherwise SuiteSparse LDL. PARDISO defaults to 16 threads, while `--pardiso-threads 0` explicitly requests the oneMKL default. Use `--linear-solver suitesparse-ldl` to force the portable block-aware LDL backend, `--linear-solver cholmod` for CHOLMOD's cost-aware AMD/METIS policy, or `--linear-solver eigen-cg` for Eigen conjugate gradient with incomplete-Cholesky preconditioning. All four backends share the same lower-triangular Hessian assembly and boundary handling.

Each Newton convergence decision uses the direction just solved from the current gradient/Hessian. A converged direction exits before CCD and line search; the strict acceptance rule remains `E_trial < E0` with `armijoCoefficient=0`. This prevents a numerically vanished direction from being halved repeatedly only because parallel energy summation fluctuates in the last bit. On the two-bunny2 scene, five independent five-step runs reduced the first 25 frames from 417 energy backtracks to 0, made all five final states identical, and reduced median total time from 4.429 s to 3.396 s. The final convergence-check factorization is still counted, so `numeric_factorizations` can be one greater than the number of accepted Newton updates.

On the development machine, alternating paired Release runs reduced five-step wall time by about 5.9% on cloth-bunny and 8.6% on twisting-mat versus the installed non-supernodal CHOLMOD build. These numbers are machine/dependency specific; use `scripts/benchmark.py` to compare locally.

For the larger two-`bunny2` scene (38,386 vertices, 115,158 DOF), one complete time step with three Newton solves had a three-run median of 1.284 s with 16-thread PARDISO, versus 8.253 s with CHOLMOD and 12.568 s with SuiteSparse LDL on the development machine. Here one `--steps 1` means one complete simulation time step/frame, not one Newton iteration. See [`agent_docs/10_pardiso_report.md`](agent_docs/10_pardiso_report.md) for the thread sweep, phase timings, memory, 50-step contact validation, and limitations.

The project deliberately does not require CHOLMOD's GPL supernodal module. If your SuiteSparse build includes it, benchmark `cholmod` separately; its relative performance can differ from the default non-supernodal vcpkg build.

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

The benchmark launches independent processes, stores each run under `Output/benchmark/run_NN`, sums per-frame `step_ms`, and reports median/min/max. Its default `--linear-solver auto` leaves backend selection to the executable, so it follows the same PARDISO→SuiteSparse-LDL capability fallback; pass a concrete backend for controlled A/B runs. `metrics.csv` also separates assembly, linear solve, CCD, line search, and post-line-search time while recording Newton iterations, backtracks, accepted step sizes, contact counts, matrix nonzeros, and direct-solver symbolic-analysis/numeric-factorization counts. PARDISO runs additionally record phase-11 analysis, phase-22 factorization, phase-33 solve time, effective thread count, and factor nonzeros.

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
| `CPU IPC/NewtonLinearSystem.*`, `SuiteSparseLDLSolver.*`, `CholmodSolver.*`, `PardisoSolver.*` | Shared sparse assembly and SuiteSparse LDL/CHOLMOD/PARDISO/Eigen-CG backends. |
| `CPU IPC/SimulationMesh.*`, `Simulator.*` | Mesh state, IO, scene construction, and simulation ownership. |
| `CPU IPC/ViewerMain.cpp` | Fixed-function GLUT viewer. |
| `apps/cipc_headless.cpp` | Headless executable. |
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
