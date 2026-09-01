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
- CHOLMOD symbolic-factorization reuse plus a verified, cost-aware METIS nested-dissection fallback.
- Optional Eigen-CG backend using the same assembled system and boundary handling.
- Strict energy-decreasing line search, Additive CCD, and the original IPC CFL strategy.

## Dependencies
Windows:
```bash
vcpkg install eigen3 freeglut tbb openblas suitesparse metis
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

With METIS ordering enabled, configuration fails early unless the installed CHOLMOD exports the Partition functionality and `cholmod_metis` is linkable. The local SuiteSparse finder accepts `METIS::METIS`, `METIS::metis`, vcpkg's un-namespaced `metis` target, conventional `metis.h + libmetis` installations, and CHOLMOD builds with embedded Partition support. It also selects matching Debug/Release SuiteSparse libraries for multi-configuration builds.

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
build/Release/cipc_headless --scene twisting-mat --steps 1 --linear-solver eigen-cg --no-output
python scripts/benchmark.py --exe build/Release/cipc_headless.exe --repeats 5 --steps 20
```

Runtime metrics are written to `metrics.csv` in the selected output directory. See [`agent_docs/README.md`](agent_docs/README.md) for architecture, algorithms, known issues, and optimization results.

Scene construction is fresh by default. Checkpoint loading and writing are opt-in through `--resume` and `--write-checkpoints`.

The CPU LBVH broad phase is the default; use `--broad-phase spatial-hash` for the optimized legacy backend.

CHOLMOD is the default Newton linear solver. Its default build enables the cost-aware AMD/METIS policy: AMD is retained for small or favorable systems, while METIS nested dissection is available when AMD predicts excessive fill or work. Eigen conjugate gradient with incomplete-Cholesky preconditioning remains available through `--linear-solver eigen-cg`; both backends share the same lower-triangular system assembly and boundary handling.

## CPU IPC benchmark usage

Use a Release, headless-only build when comparing CPU implementations:

```bash
cmake -S . -B build/benchmark \
  -DCMAKE_BUILD_TYPE=Release \
  -DCIPC_BUILD_VIEWER=OFF
cmake --build build/benchmark --config Release --parallel

python scripts/benchmark.py \
  --exe build/benchmark/Release/cipc_headless.exe \
  --scene cloth-bunny \
  --broad-phase lbvh \
  --linear-solver cholmod \
  --steps 20 \
  --repeats 5
```

The benchmark launches independent processes, stores each run under `Output/benchmark/run_NN`, sums per-frame `step_ms`, and reports median/min/max. `metrics.csv` also separates assembly, linear solve, CCD, line search, and post-line-search time while recording Newton iterations, backtracks, accepted step sizes, contact counts, matrix nonzeros, and CHOLMOD analysis/factorization counts.

For fair comparisons:

- use the same scene, material file, timestep, bending mode, broad-phase semantics, linear tolerance, METIS-ordering setting, and number of steps;
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
| `CPU IPC/IPCSolver.*` | Time stepping, Newton loop, CCD/CFL integration, and strict energy line search. |
| `CPU IPC/ContactMechanics.*` | Contact distances, barriers, derivatives, and feasible self-contact steps. |
| `CPU IPC/CollisionBroadPhase.*`, `LBVH.*` | SpatialHash and Linear BVH broad phases. |
| `CPU IPC/Elasticity.*`, `HingeBending.*` | Solid/cloth elasticity and bending models. |
| `CPU IPC/NewtonLinearSystem.*`, `CholmodSolver.*` | Sparse assembly and CHOLMOD/Eigen-CG backends. |
| `CPU IPC/SimulationMesh.*`, `Simulator.*` | Mesh state, IO, scene construction, and simulation ownership. |
| `CPU IPC/ViewerMain.cpp` | Fixed-function GLUT viewer. |
| `apps/cipc_headless.cpp` | Headless executable. |
| `scripts/` | Benchmark and video utilities. |

## Manual smoke checks

The repository intentionally ships no test suite. After changes, validate the product executables directly:

```bash
build/Release/cipc_headless --scene cloth-bunny --steps 1 --no-output --broad-phase lbvh
build/Release/cipc_headless --scene cloth-bunny --steps 1 --no-output --broad-phase spatial-hash
build/Release/cipc_headless --scene twisting-mat --steps 1 --no-output --linear-solver eigen-cg
build/nonquadratic/Release/cipc_headless --scene cloth-bunny --steps 1 --no-output
```
