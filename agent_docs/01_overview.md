# 01 — 构建、架构与核心数据结构

本版本定位为**高度 CPU 优化的 IPC 实现与 CPU IPC 对比基准**：提供无窗口产品入口、分阶段metrics、重复进程benchmark、LBVH/SpatialHash A/B和PARDISO/优化CHOLMOD/Eigen-CG A/B。对比时必须同时核对轨迹、接触、Newton/线搜索和最小距离，不能通过改变物理问题换取表面加速。

## 构建与运行

推荐构建入口（依赖细节见 `README.md` 与 `12_wsl_ubuntu_build.md`）：

```powershell
# Windows：默认同时构建 viewer 和 headless
.\build.cmd
```

```bash
# WSL/Ubuntu：默认同时构建 viewer 和 headless
./build.sh
```

两条入口均按“显式参数 → `VCPKG_ROOT` → 项目固定 revision 的托管 checkout”解析 vcpkg。Windows 会检查 MSVC C++ workload，并从 PATH、Visual Studio 或 vcpkg tools 自动解析 CMake；Ubuntu 会检测并按需通过 apt 安装基础编译包，`--no-system-packages` 可关闭该行为，系统 CMake 低于 3.23 时自动使用 vcpkg 私有版本。Windows 的短路径 checkout 位于 `%LOCALAPPDATA%/CPU-IPC/vcpkg-<revision>`，WSL 位于 `${XDG_CACHE_HOME:-$HOME/.cache}/cpu-ipc/vcpkg`；没有任何提交到仓库的本机绝对路径。vcpkg 一次准备 Eigen、`tbb[core]`、METIS、oneMKL 与可选 FreeGLUT；SuiteSparse 7.14.0 按 SHA-512 校验下载，只构建 Config/AMD/CAMD/CCOLAMD/COLAMD/CHOLMOD，不再安装无用 OpenBLAS/hwloc。Windows 产物位于 `build/cpu-ipc/<Configuration>`，WSL 使用独立的 `build-wsl/cpu-ipc`，避免两套 ABI/缓存互相污染。完成一次配置后仍可在相应 build 目录手工增量构建：

```powershell
cmake --build build/cpu-ipc --config Release --parallel
```

```bash
cmake --build build-wsl/cpu-ipc --parallel
```

- `CMakeLists.txt` 显式列出核心源文件并拆成三个 C++14 产品 targets：`cipc_core`（无 OpenGL）、`cipc_headless`、可选 `cipc` viewer；不注册 CTest/测试 target。
- CMake 选项：
  - `CIPC_BUILD_VIEWER=ON` — viewer；关掉后配置阶段无需 OpenGL/GLUT。shader/GLEW 路径已删除。
  - `CIPC_ENABLE_FRICTION=ON` — 定义 `USE_FRICTION`。
  - `CIPC_ENABLE_QUADRATIC_BENDING=ON` — 定义 `USE_QUADRATIC_BENDING`。
  - `CIPC_ENABLE_METIS_ORDERING=ON` — 构建/链接 CHOLMOD Partition/METIS 供实验；生产配置按实测固定 AMD，因为它在三个项目矩阵上都快于强制 METIS 与 multi-method AUTO。
  - `CIPC_ENABLE_PARDISO=ON` — 推荐保持开启；找到oneMKL时以PARDISO作为默认主后端。PARDISO不可用或MSVC Debug时自动选择优化CHOLMOD。
  - `CIPC_CHOLMOD_ROOT` — 指向性能版 SuiteSparse/CHOLMOD bundle prefix；Windows 标准 `build/cholmod-mkl-install` 会自动识别，WSL 一键脚本显式传入 `build-wsl/cholmod-mkl-install`。配置会从同一 prefix 解析全部 SuiteSparse 组件，验证 supernodal 与 METIS 符号，并让 CHOLMOD 使用 oneMKL LP64/TBB。
  - `CIPC_REQUIRE_OPTIMIZED_CHOLMOD=ON` — 默认拒绝误连普通 OpenBLAS CHOLMOD；根目录 `build.cmd` / `build.sh` 可一键完成依赖、连接和编译。
  - `CIPC_ASSETS_DIR` = `<repo>/Assets/`（活，`Simulator.cpp` 读参数文件/网格用）
  - `CIPC_OUTPUT_DIR` = `<repo>/Output/`（活，由 `RuntimePaths` 统一管理日志、检查点、表面和截图输出）
- SuiteSparse 查找器会校验 `cholmod_metis`，兼容 `METIS::METIS`、`METIS::metis`、vcpkg 的 `metis` target、Linux 常规 `metis.h + libmetis` 以及内嵌 Partition 的 CHOLMOD；性能 bundle 的全部组件强制来自 `CIPC_CHOLMOD_ROOT`，Windows Debug/Release 通过 DLL 边界共用 Release bundle。配置额外校验 `cholmod_super_numeric`；只有显式兼容诊断的 system fallback 才选择 OpenBLAS/系统 BLAS-LAPACK。
- MSVC 的 `/bigobj` 只施加到 `cipc_core`（`ContactMechanics.cpp` 的生成代码需要），不再污染全局 `CMAKE_CXX_FLAGS`。

运行：`cipc`打开GLUT窗口，空格开始/暂停。`cipc_headless --steps N --broad-phase lbvh`写逐帧`metrics.csv`；默认有PARDISO时使用PARDISO，否则自动使用优化CHOLMOD。`--linear-solver cholmod|pardiso|eigen-cg`可显式覆盖。PARDISO默认16线程；CHOLMOD默认按nnz自动4/8线程。一个`--steps 1`是一个完整时间步/帧。——见`06_app_layer.md`。

## 整体架构

```
cipc viewer / cipc_headless
   └─ cipc_core::FEMSimulator                  场景搭建 + 参数加载
        └─ solveIPCStep + IPCSolverContext     每时间步求解 + 实例级指标/恢复状态
                  ├─ κ 外层循环: solveBarrierSubproblem × N
                  │    └─ Newton 迭代:
                  │         computeGradientAndHessian  (弹性+软边界+障碍+摩擦)
                  │         NewtonLinearSystem::solve  (PARDISO / optimized CHOLMOD / Eigen-CG)
                  │         可行步长 (ACCD / CFL / 地面射线)
                  │         lineSearch (回溯 + 穿透防护)
                  │         postLineSearch (κ 自适应)
                  ├─ BoundaryConditionOps (软目标冻结 + Dirichlet/CCD 驱动)
                  ├─ Friction::initialize (每步开始时冻结 λ 与切空间基)
                  └─ updateVelocity / updateInertialTarget
```

- `FEMSimulator::simulateStep` 直接对 `SimulationModel::meshes.front()` 调 `solveIPCStep`；已删除仅做一次转发、没有第二实现的 integrator 层。当前仍是单 mesh 求解。
- 弹性数值内核在 `Elasticity.cpp`，旋转保持 SVD 在 `RotationAwareSVD.cpp` + 自包含 `ImplicitQR3x3SVD.h`；边界在 `BoundaryConditions.cpp`；摩擦在 `Friction.cpp`；Newton 线性系统在 `NewtonLinearSystem.cpp`；`IPCSolver.cpp` 负责时间步、Newton、CCD 与线搜索编排。
- 碰撞宽阶段默认使用 `LBVH.cpp`，`CollisionBroadPhase.cpp` 的 SpatialHash 保留为回归后端；距离、障碍与接触装配在 `ContactMechanics.cpp`，Additive CCD 在 `AdditiveCCD.cpp`。

## 核心数据结构

### `mesh3D`（`SimulationMesh.h`）——全项目中心

一个 `mesh3D` 同时容纳**四面体实体**和**布料三角形**（`ClothOverBunny` 就是布料 + bunny 合并在同一网格）。

| 字段 | 含义 |
|---|---|
| `vertexes / v_rest / velocities / V_prev` | 当前位置 / 静止位置 / 速度 / 上一步位置 |
| `tetrahedras` (Vector4i) / `triangles` (Vector3i) | 单元 |
| `surface` (Vector4i) | 表面三角形，**第 4 分量 = 所属 tet id**（布料面为 0） |
| `surfVerts / surfEdges` | 表面顶点（全局 id）、无向表面边 |
| `DM_tetrahedra_inverse` / `DM_triangle_inverse` | 静止形变矩阵逆 |
| `tetrahedraPFPX` (9×12) / `trianglePFPX` (6×9) | 初始化时预计算的静止几何导数算子，Newton 中只读 |
| `volum / areas / masses` | tet 体积；**布料 `areas` 已乘厚度=体积**；集中质量 |
| `inertialTarget` | 惯性预测位置 `x̃ = V_prev + v·dt + g·dt²`（`IPCSolver.cpp`） |
| 材料参数 | `density, YoungModulus, PoissonRate, lengthRate, volumeRate, friction, clothThicness, stretchStiffness, shearStiffness, plateRigidity=Et³/[12(1−ν²)], strainRate` |
| IPC 参数 | `Hhat`（**平方**激活距离 d̂²）、`Fhat`（摩擦 εv² 基准）、`Kappa`（障碍刚度）、`dTol`（**平方**收敛距离）、`bboxDiagSize2`、`IPC_dt`、`averageEdgeLenth`（= 平均表面边长 / 3，空间哈希格宽） |
| 活动集 | `Self_ActiveSet`（`EncodedContact`，PT/EE/PP/PE 混合）、`Self_EE_ActiveSet` + `Self_EEeIe_ActiveSet`（近平行 EE，mollified）、`Self_CCD_ActiveSet`（pair<int,int>）、`Environment_ActiveSet`（地面接触顶点 id） |
| 摩擦滞后量 | `Self_lambda_lastH / Environment_lambda_lastH`（法向力幅值）、`MMDistCoord`（最近点坐标 β/γ/η）、`MMTanBasis`（3×2 切空间基）、`*_activeSet_lastH`（上一步活动集快照） |
| 弯曲 | quadratic：`quadBendingInfo` 含 `Q⊗I₃`；hinge：`hingeBendingInfo` 含四顶点、静止角和 `l0/(h0+h1)`；二者共用 `plateRigidity` |
| 边界 | `BoundaryConditionSet` 分离 `dirichlet` 与 `soft` 集合；`boundaryTypes` 由 `VertexBoundaryType` 命名：0=自由、1=Dirichlet、≥2=外部 prescribed collider。硬回调 `updateDirection(current,rest,step,alpha,dt)` 返回 `x_new=x−p` 的方向；软回调 `updateTarget(current,rest,step,dt)` 返回绝对目标位置 |
| 开关 | `use_barrier, apply_gravity(默认true), is_quasi_static(默认false)` |
| 恢复 | `resumedFromCheckpoint`：仅当位置与 x̃ 检查点均通过顶点数校验时为 true；计时/κ 只在此时恢复 |

### `BHessian`（`ContactMechanics.h`）——分块 Hessian

弹性/障碍/摩擦的逐单元 Hessian 全部以**小块**形式追加进 4 个列表：

- `H3x3 + D1Index`：1 顶点（软边界、地面障碍、地面摩擦）
- `H6x6 + D2Index`：2 顶点（PP 接触）
- `H9x9 + D3Index`：3 顶点（PE 接触、布料三角形弹性）
- `H12x12 + D4Index`：4 顶点（PT/EE 接触、tet 弹性、弯曲）

`toTriplets(boundaryTypes, output)`用两遍并行count/prefix/fill，只输出全局下三角、自由顶点项和非零数值；不会生成占位零。输出buffer属于`NewtonLinearSystem`并跨Newton/时间步保留容量；PARDISO、CHOLMOD与Eigen-CG共用该矩阵和RHS。

### `EncodedContact` 接触编码（`EncodedContact.h`）

4 整数 id，词典序（用作 `std::map` key 计重数）。活动集条目解码规则在 `CollisionBroadPhase.cpp`：

- `[0] >= 0` → EE 对，4 个顶点 id。
- `[0] < 0` → 点基元，点 id = `-[0]-1`；然后 `[2]<0` → PP，`[3]<0` → PE，否则 PT（`[1..3]` 为三角形顶点，-1 填充）。
- PP/PE 条目经 `constraintCounter` 按 `EncodedContact` 去重，最终 `[3] = -count`：**重数 = `-contact[3]`**，在能量/梯度/λ 中作乘数。
- EE 构建时若近平行（`EECrossSqNorm < eps_x`，`eps_x = 1e-3·|restE_I|²·|restE_J|²`），`[3]` 编码为 `-eJ-2` 或 `-vertexId - surfEdges.size() - 2`，进 `Self_EE_ActiveSet` 走 mollifier 路径。

## 全局约定（改代码前必须内化）

1. **距离全部是平方值**。`Hhat = d̂²`（默认 9e-8 或参数文件 1e-3 被平方成 1e-6，再乘 `bboxDiagSize2`，`Simulator.cpp:161-163,309`）。`dTol = 1e-18·bboxDiagSize2`。比较用 `d < Hhat`。
2. **位置更新是减法**：`x ← x - α·searchDir`（`stepForward`，`IPCSolver.cpp`）。扫掠包围盒与 CCD 位移都取 `-searchDir`。
3. **向量化一律列主序**（`ElasticityMath::vectorizeColumnMajor`），与 `computePFPX*` 布局及 Eigen 默认一致。
4. `computeRotationAwareSVD` 约定：σ 按绝对值降序，仅末项可负；翻转 tet 的 `J<0` 是刻意的，SNK 全程无 J 钳制。
5. 活动集只在 `d < Hhat` 时才评估障碍函数——障碍函数对 `d ≥ dHat` 无定义保护，调用方负责预过滤。
6. TBB 并行到处可见：梯度累加仍用每顶点 `spin_mutex`；静态/CCD 主查询已改为排序 vector，减少临时 hash 与候选顺序漂移，但并行浮点归并仍不保证 bitwise 确定。
