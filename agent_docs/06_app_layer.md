# 06 — 应用层：可视化、场景、线性求解器、网格 IO

构建已分层：`cipc_core` 包含 Simulator/elasticity/IPC/contact/solver 且不依赖 OpenGL；`cipc` 仅是 viewer；`cipc_headless` 是无窗口 CLI。`CIPC_BUILD_VIEWER=OFF` 可完全跳过 OpenGL/GLUT 查找。项目已无 GLEW 依赖。

## 1. ViewerMain.cpp —— 固定管线查看器

- `main`：GLUT 初始化（1000×1000 窗口 "CPU IPC"）→ `initializeViewer` → 注册 display/reshape/keyboard/mouse/motion/idle 回调。
- shader 目录、永远为 false 的 shader 分支、VBO/VAO 全局量和 GLEW 初始化/链接已全部删除；viewer 明确只使用当前实际工作的 OpenGL 固定管线。
- 主循环：`idle` 请求重绘；`display` 渲染后若未暂停则调用一次 `simulator.simulateStep(step)`，即每显示帧一个时间步。
- `drawSimulationSurface` 只画当前 `SimulationModel::meshes.front()`：红色表面、白色表面边、灰底与线框包围盒。

**键盘**：`空格` 暂停/继续；`w/s/a/d/q/e` 平移视图；`9` 切换 OBJ；`/` 切换截图；鼠标左键拖动旋转。`1-5/k/f/m/0` 等无人消费的特效开关已删除。

**输出**：
- `saveSurfaceMesh`：默认关闭，按 `9` 切换；开启后写 `<repo>/Output/saveSurface/surf_NNNNN.obj`，直接只读引用 `SimulationModel`。
- `saveScreenshotIfDue`：默认关闭，按 `/` 切换；每 10 步写 BMP。实现会处理 RGB→BGR 和 BMP 行对齐；配合 `scripts/generate_video.py`。
- `RuntimePaths` 在启动时自动创建 `Output/tempData`、`Output/saveSurface`、`Output/saveScreen`；路径由 `CIPC_OUTPUT_DIR` 决定，与 CWD 无关。

## 2. Simulator.cpp —— 场景搭建与求解编排

场景辅助函数均位于匿名 namespace，不再污染全局符号；公开职责只保留 `FEMSimulator` 的初始化、碰撞集重建与逐步求解。

### 参数文件解析（`loadSettings`）

读 `CIPC_ASSETS_DIR + "scene/parameterSetting.txt"`。现按 key/value 解析，顺序任意；未知、重复、缺失键或非法数值会抛出明确异常。正值/非负值、`ν∈(-1,0.5)` 和 collision flag∈{0,1} 均有范围校验：

| 行 | 文件键 | 目标字段 | 现值 |
|---|---|---|---|
| 1 | volume_mesh_density | `density` | 1000 |
| 2 | volume_mesh_YoungsModulus | `YoungModulus` | 1e4 |
| 3 | poisson_rate | `PoissonRate` | 0.49 |
| 4 | friction_rate | `friction` (μ) | 0.4 |
| 5 | triangle(cloth)_mesh_thickness | `clothThicness` | 0.001 |
| 6 | triangle(cloth)_mesh_YoungsModulus | `clothYoungModulus` | 1e4 |
| 7 | triangle(cloth)_mesh_shearYoungsModulus | `shearYoungModulus` | 1e1 |
| 8 | triangle(cloth)_mesh_bendingYoungsModulus | `bendYoungModulus` | 1e4 |
| 9 | triangle(cloth)_mesh_density | `cloth_density` | 200 |
| 10 | strainRate | `strainRate` | 1e2 |
| 11 | enable_collision_handling | `use_barrier` | 1 |
| 12 | drag_coeff | `drag_coeff` | 0.0 |
| 13 | ipc_time_step | `IPC_dt` | 0.01 |
| 14 | Newton_solver_threshold | `Newton_Solver_Threshold` | 1e-2 |
| 15 | IPC_ralative_dHat | `Hhat`（随后平方） | 1e-3 |

后处理：`Hhat = relativeDHat²` → 统一 `updateMaterial` 做 Lamé/SNK/布料换算。`Fhat/Kappa/dTol = 1e-6/0/1e-18` 仍不从文件读。文件不存在时回退 `applyDefaultSettings`；默认材质也走同一换算函数，不再有 shear/bending 分母差异。

### 场景

- **`ClothOverBunny` / `buildClothOverBunnyScene`（默认）**：布料盖 bunny。当前载入 `Assets/triangleMesh/planes/plane1024.obj`（1024 顶点、1922 个三角形，scale 1），绕 X 转 π/2；按编译选项预计算 quadratic 或 hinge bending；再载入 `Assets/tetrahedraMesh/bunny.msh`（scale 0.5，offset (0,−0.5,0)）追加进同一 `mesh3D`。
- **`TwistingMat` / `buildTwistingMatScene`**：扭转垫（准静态，`is_quasi_static=true`，无重力），`ipcmesh/mat40x40.msh`，两侧顶点 `boundaryTypes=1` + `update_hard_constraint_functor` 旋转驱动。通过 `SimulationScene` 的 switch 可选；回调的 `alpha` 是 double。
- **`Bunny2` / `buildBunny2Scene`**：参考 `GPU_IPC/GPU_IPC/gl_main.cpp` 的 `initScene1`，两次追加 `Assets/tetrahedraMesh/bunny2.msh`，两只都取 scale=0.2，offset 分别为 `(0,0.65,0)` 与 `(0,0,0)`，`YoungModulus=1e5`。合并后是 38,386 顶点、159,870 tet、41,664 表面三角形；用于大稀疏系统/碰撞压力基准，不替换默认场景。

### `buildModels`（`:283-336`）全流程

1. 接收 `SimulationOptions{scene,resume,write files,checkpoint,verbose,broad phase,linear solver}`；恢复和 checkpoint 写入默认均为 false，必须显式 opt-in；`RuntimePaths::initialize` → `loadSettings` → 选择场景；重复 build 会先清理旧 model；
2. 按 vector 尺寸重数 `vertexNum/tetrahedraNum/triangleNum`（加载器会覆盖计数，必须重数）；
3. `initMesh3D`（质量、Dm 逆、面积）；`v_rest = V_prev = vertexes`；`updateInertialTarget`；
4. `bboxDiagSize2 = (maxConer−minConer).squaredNorm()`；`Hhat/Fhat/dTol` 各乘之；
5. 加入 `SimulationModel::meshes`；`calculateSurface()`；
6. `averageEdgeLenth = Σ边长 / (3·边数)`（/3 是刻意的 IPC dHat 启发式）；
7. `load_tetTempData()`：若 `Output/tempData/vertex.txt + vertexXtile.txt` 存在、两者完整且顶点数与当前网格相同则热启动，并设置 `resumedFromCheckpoint=true`；
8. 设置 broad-phase backend 并 `rebuildCollisionSets`；把线性求解选项复制进实例级 `IPCSolverContext`。

`simulateStep` 直接调用 `solveIPCStep(model_.meshes.front(), broadPhase_, ground_, solverContext_)` 并返回 Newton 迭代数。

## 3. cipc_headless 与指标

示例：

```bash
cipc_headless --scene cloth-bunny --steps 20 --broad-phase lbvh --output Output/run
cipc_headless --scene twisting-mat --steps 1 --linear-solver eigen-cg --no-output
cipc_headless --scene bunny2 --steps 1 --linear-solver pardiso --pardiso-threads 16 --no-output
```

- `--steps N` 的单位是完整仿真时间步/帧；一个 step 内可能有多次 Newton 迭代和数值分解。
- 默认不恢复、不写 checkpoint、使用 LBVH；可用 `--resume`、`--write-checkpoints` 显式开启。viewer 同样默认 fresh，不会让 `Output/tempData` 覆盖 cloth+bunny 初态。
- `--broad-phase lbvh|spatial-hash` 用于碰撞宽阶段 A/B；`--linear-solver suitesparse-ldl|cholmod|pardiso|eigen-cg` 选择线性后端，默认 SuiteSparse LDL；PARDISO 线程上限默认 16，`--pardiso-threads 0` 使用 oneMKL 默认；`--verbose` 打印 Newton/线搜索和 PARDISO phase 细节。
- 诊断选项：`--diagnose-line-search` 输出方向一阶/二阶 Taylor 对比；`--disable-barrier` 仅用于无接触隔离。临时的 `--friction-scale` 已删除，正式接口不能绕过场景材料参数改变摩擦。
- 每帧 `metrics.csv` 字段包括五阶段耗时、Newton/κ、总/能量/穿透回退、单 Newton 最大回退与 `Newton>2` 数、mean/min/max α、碰撞、最小平方距离、活动集、nnz、symbolic analyze 与 numeric factorize 次数；末尾的 `pardiso_analysis_ms/pardiso_factorization_ms/pardiso_solve_ms/linear_solver_threads/factor_nnz` 用于 PARDISO 细分，其他后端为 0。
- 最终 `RESULT` 输出位置和、平方范数和及最后一帧关键指标，供脚本/CI 解析。
- `scripts/benchmark.py` 执行多次独立运行并报告总 time-step time 的 median/min/max；支持 `--scene`（含 bunny2）、`--broad-phase`、`--linear-solver`、`--pardiso-threads`、`--steps`、`--repeats` 和 `--output`。

## 4. NewtonLinearSystem 与 Solver

`NewtonLinearSystem.cpp/.h` 是 Newton 线性系统的唯一入口：

- `assemble` 从 `BHessian` 生成下三角 triplet，加入质量/阻尼对角并按 `boundaryTypes` 将固定顶点 RHS 清零。
- SuiteSparse LDL、CHOLMOD、PARDISO 与 Eigen-CG 共用同一 `Eigen::SparseMatrix`、RHS 和解向量，再统一 scatter 回每顶点方向；这保证切换后端不会改变约束或装配语义。
- 默认 SuiteSparse LDL；可选 Eigen-CG 使用 `Eigen::ConjugateGradient<SparseMatrix, Lower, IncompleteCholesky>`，容差 `1e-6`、最大 10000 次，可由 `LinearSolverOptions` 设置。
- CHOLMOD 后端随默认构建开启 `CIPC_ENABLE_METIS_ORDERING`：配置阶段要求 SuiteSparse `Partition` 并链接验证 `cholmod_metis`；运行时保持 `nmethods=0`、指定 `default_nesdis=0` 和 weighted postorder，即先由 AMD 估计 fill/work，仅在代价较高时再尝试 METIS。关闭选项时设置 `nmethods=1 + CHOLMOD_AMD`，形成真正的 AMD-only A/B。
- PARDISO 后端由 `CIPC_ENABLE_PARDISO` 控制；oneMKL 不存在时 CMake 只关闭该后端，不影响默认 LDL。`PardisoSolver` 直接调用 phase 11/22/33，使用 SPD `mtype=2`、LP64/zero-based 索引、TBB threading 和 in-core factorization。lower CSC 与 upper CSR 的等价布局避免额外矩阵转置；RHS 使用自有副本。符号分析、METIS permutation、factor 与 phase workspace 随 `IPCSolverContext::linearSystem` 跨时间步复用，稀疏模式变化才重新分析，并用 1.2× fill 阈值自适应触发 fresh ordering。
- Windows vcpkg oneMKL 静态链接会使 Release headless 约 74.0 MB（70.6 MiB）；其静态包不兼容 `/MDd`，所以 MSVC Debug 自动不定义 PARDISO，选择该后端会报告 unavailable。Release/RelWithDebInfo 是正式 benchmark 配置。
- Eigen-CG 可用 `cipc_headless --scene twisting-mat --steps 1 --no-output --linear-solver eigen-cg` 做手工 smoke；它不要求与直接法 bitwise 相同。

`SuiteSparseLDLSolver.cpp/.h` 实现默认 CPU 直接法：

- 从共享下三角矩阵提取每顶点 3-DOF block graph，在 block graph 上做 AMD，再将 permutation 展开回标量 DOF，保留 xyz 局部性。
- 第一次或结构变化时显式生成 `PAPᵀ` 上三角 CSC、建立原 lower CSC value index 到 permuted CSC value index 的一一映射并执行 `ldl_symbolic`；同结构 Newton 只刷新数值。
- 每轮调用 `ldl_numeric + lsolve + dsolve + ltsolve`；要求所有 D 对角有限且严格为正，失败立即抛错，不把不定方向交给线搜索。
- analyze/factorize 计数与 CHOLMOD 共用 `IPCStepStats` 字段；20 步 A/B 中两直接法的逐帧整数指标完全一致。

`CholmodSolver.cpp/.h` 只实现 SuiteSparse `CholmodSolver`：

- 类已改为 RAII、不可复制：构造时检查 `cholmod_start`，析构释放当前 factor 并 `cholmod_finish`；不再持有被重绑定的 CHOLMOD-owned sparse/dense 缓冲区。
- `set_pattern(SparseMatrix)` 把压缩 Eigen CSC 的 outer/inner/value 数组复制到稳定的自有缓冲，并用非 owning `cholmod_sparse` view 传给 CHOLMOD。
- 每次都刷新数值；仅当矩阵维度或 outer/inner 索引变化时释放 factor。`solve` 在 factor 不存在时执行 `cholmod_analyze`，随后每轮执行数值 `cholmod_factorize + cholmod_solve2`；solution/Y/E dense workspace 跨 Newton 复用。
- `IPCSolverContext` 生命周期内复用同一 `NewtonLinearSystem` 及其直接求解器，因此接触拓扑不变的相邻 Newton、κ 子问题和时间步都可以跳过符号分析；接触导致稀疏模式变化时自动重分析。
- 暴露 analyze/factorize 计数给 `IPCStepStats`；当前首帧为 2 次 analyze / 9 次 factorize，说明大部分 Newton 已复用 symbolic 数据。
- `preFactorize / solve_with_preFactorize` 保留给多 RHS；所有入口都有矩阵形状、RHS 尺寸和 CHOLMOD 失败检查。
- 已移除未定义的三元组 `set_pattern`、`IJ2aI` 和无用 dense/sparse 指针成员。

## 5. SimulationMesh.cpp —— 网格 IO 与仿真状态

### 加载

- `load_tetrahedraMesh`：按行触发解析并取节点/单元行尾部 token，支持 append 多网格。文件不存在会返回 false，由场景构建抛出可捕获异常，不再 `exit(-1)`。
  - 资产：`bunny.msh`（Gmsh 2.2，1869 节点 7356 单元）；`bunny2.msh` 单份为 19,193 节点、79,935 tet；`ipcmesh/mat40x40.msh` 自称 4.1 但体是简化格式（靠“取尾部 token”侥幸解析成功）。`ipcmesh/` 下还有 Armadillo13K、rod、sphere5K、torus 等可用场景网格。
- `load_triangleMesh`（OBJ）：`v` 行 scale+offset；`f` 行扇形三角化（支持 `f a/b/c` 与 `f a b c`）；`type==2` 全顶点 `boundaryType=2`，`type==3` 面直接进 `surface`，默认进 `triangles`；末尾建 `tri_edges/_adj_points`（仅 2 三角共享边）。旧 PCG 的逐顶点 3×3 `Constraints` 已删除。资产：`CMU/plane{9,100,1024,200000}.obj`、`cloth7.obj`、`tricloth.obj`、`newtubing.obj`。
- `getSurface`（`:349-431`）：tet 4 面×6 排列哈希去重找边界面；按对顶点法线测试定向**朝外**；存 `Vector4i(v0,v1,v2,tetId)`；布料面 tetId=0 追加；`surfVerts` 按首现序唯一化；`surfEdges` 用 `std::set` 去重无向边。
- `mesh3D::InitMesh / load_test(1/2)`：程序化立方体，当前产品场景不调用。
- 已删除无人调用的 `mesh2D` 示例、`mesh_obj/model_obj/fiber_obj` 特效加载器、空 `model_tet::load_model` 和旧 EKF 灵敏度接口。程序化 unit cube 数据仅保留在 `SimulationMesh.cpp` 内部，不进入公共头文件。

### 输出 / 状态

- `output_tetTempData / load_tetTempData`：`Output/tempData/vertex.txt`（位置，恢复时也写回 `V_prev`）+ `vertexXtile.txt`（x̃）。保存使用 17 位精度；加载先读入临时数组并要求两个文件恰好包含 `vertexNum` 行，成功后才整体提交状态。
- `output_tetrahedraMesh`：节点和 tet 引用均写 1-based，并补 `$EndNodes/$EndElements`，已修复与自身加载器 round-trip 的 off-by-one。

## 6. Python 工具

- `scripts/generate_video.py`：截图目录 → 视频，支持 `--images/--output/--fps`，默认读取 `Output/saveScreen`。
- `scripts/merge_videos.py`：两组等长、同尺寸截图 → 横向比较视频，支持 `--left/--right/--output/--fps`。两脚本都通过 `__file__` 定位仓库。

## 7. 时间/日志文件

- `Output/timeCost.txt`：每帧覆盖累计 ttime0..4、总迭代、帧号、κ。
- `Output/tempData/timeCost.txt` + `vertex*.txt`：每 10 帧持久化；只有网格检查点完整且顶点数匹配时才恢复计时/κ。相同拓扑的旧检查点仍会被视为有效，重跑实验前仍应清理 `Output/tempData/`。
