# 07 — 坑、已知问题与死代码清单（改代码前必读）

## A. 当前仍存在的问题 / 风险

1. **相同拓扑的旧检查点仍可能误恢复**：检查点会校验向量数量，但尚未保存场景、拓扑 hash 或生效参数 hash。同顶点数/同拓扑的不同实验仍需使用不同 `--output` 或清理 `tempData`。
2. **性能比较仍须记录 broad phase，但物理指标必须一致**：两后端原始候选 superset 可以不同，进入 partial/full ACCD 的 pair 已由统一 AABB 条件归一化。回归要求逐帧 Newton、回退、接触数一致，α/轨迹只允许并行浮点容差；若再次分叉应视为 bug。
3. **仍非 bitwise 确定**：主 PT/EE 候选已排序，文件级全局状态也已移除，但 TBB 浮点归并与并行梯度加法顺序仍可产生约 1e−11～1e−14 的差异。回归使用容差，不比较二进制 hash。
4. **迭代上限尚未配置化**：Newton 10000、κ 64、穿透回退 64、ACCD 10000 目前是代码常量；已有失败出口，但若要针对场景调节，应进入 `SolverSettings`/配置文件。
5. **`mesh3D` 仍是 god object**：静止拓扑、材料、动态状态、接触历史和 workspace 伴随数据尚未完全分离；当前仅把 Newton workspace、solver context 和 broad phase 独立出来。
6. **PARDISO 是推荐默认Release后端**：Windows静态oneMKL与`/MDd`不兼容，因此MSVC Debug自动使用优化CHOLMOD。两条生产路径都需要oneMKL，一键脚本会统一准备。
7. **性能版 CHOLMOD 改变许可边界**：一键脚本从固定 SuiteSparse 源码构建 Config/AMD/CAMD/CCOLAMD/COLAMD/CHOLMOD 共享 bundle，启用 GPL supernodal 并把静态 oneMKL 嵌入 CHOLMOD。配置必须通过 `CIPC_CHOLMOD_ROOT` 接入完整 bundle 且验证 `cholmod_super_numeric`；不能把普通 OpenBLAS/simplicial 包的结果标成 MKL supernodal。二进制分发需遵守 GPL-2.0-or-later。
8. **Windows 与 WSL 构建目录不可复用**：两边的编译器、路径语法、库格式和 ABI 不同。Windows 固定使用 `build`，WSL 一键脚本固定使用 `build-wsl`；不要从 WSL 配置已有的 Windows build tree，也不要把 `/mnt/c/...` 或盘符路径写入仓库配置。vcpkg 只能复用同一平台的 checkout/installed tree。
9. **一键脚本只能代替库配置，不能消除平台编译器**：Windows 仍需 64-bit MSVC C++ workload；脚本会检测并给出 winget 命令。Ubuntu 默认会通过 sudo/apt 补齐缺少的编译工具和 viewer X11/OpenGL 开发包；CI 或受管环境应传 `--no-system-packages`，此时只报告精确缺项、不修改系统。自动 vcpkg 固定 revision；只有用户显式传 `VCPKG_ROOT`/路径时才允许替换，维护者升级 `vcpkg-revision.txt` 后必须重新做双平台首次安装。

## B. 保留但默认未接入的安全路径

1. `TetInversionGuard::limitStepToPreventTetInversion` 保留原 tet 防翻语义：求 `x−αp` 下有向体积降至当前值指定比例的最小正根。默认 IPC 流程仍只使用 CCD+CFL；启用时应手工验证代表 tet（20% 阈值应给出 α=0.8）并记录轨迹变化。

## C. 死代码清单（编译了但无调用方）

**求解器/算法路径：**

- `compute_g_dpt_new` / `compute_H_dpt_new`（I5 不变量障碍；当前无调用点）。
- SpatialHash 的旧 `spanSize` 钳制仍未启用。

**材料模型：**

- SNK v1：`computePEPF_StableNHK3D_double`、`project_StabbleNHK_H_3D`（含遗留 debug 输出）。
- 各向异性肌肉模型、`Iso*` 等距嵌入函数群、`__Inverse/__Inverse2x2`。

**应用层：** 当前已没有单独的休眠 shader/特效资源；viewer 只保留实际固定管线路径。

## D. 容易写错的核心约定

- **距离一律平方**；`Hhat`、`dTol`、`Fhat` 都已乘 `bboxDiagSize2`。比较写 `d < Hhat`。
- 更新是 `x − α·searchDir`；给 ACCD 的位移是 `-searchDir`。
- PP/PE 重数 = `-EncodedContact[3]`，能量、梯度、Hessian/λ 的一致缩放都要核对。
- 布料 `areas[i]` 含厚度（= 体积），能量/力直接乘它。
- `vec()` 列主序；`PFPX` 布局依赖这一点。
- 障碍函数只在 `d < Hhat` 合法；新增调用点必须先过滤活动集。
- `BHessian::toTriplets` 会按 `boundaryTypes` 过滤 Dirichlet/外部 collider 顶点，RHS 必须用同一 predicate 清零；软边界顶点必须保持 code 0，否则 `weight·I` 与梯度会被错误消元。
- 并行累加共享状态必须有明确的数据竞争策略；替换 spin mutex 时要做数值回归。
- 产品入口为 `ViewerMain.cpp` 和 `apps/cipc_headless.cpp`；核心源文件由 CMake 显式列出。项目当前不保留测试源码或 CTest target。
- CHOLMOD的排序策略会影响symbolic时间、fill-in和factorization；当前生产策略固定AMD，因为项目矩阵上的强制METIS与multi-method AUTO都更慢。Partition/METIS模块仍编入供实验。
- PARDISO 使用 `mtype=2` 且读取 upper CSR；Eigen column-major lower CSC 的 outer/inner/value 内存正好等价于同一对称矩阵的 upper CSR。不要再做一次转置，也不能把 `iparm[34]=1` 的零基索引改回一基。LP64 构建要求 `MKL_INT` 与 Eigen `StorageIndex` 都是 32-bit int。
- PARDISO 的 phase 11/22/33 状态随 `IPCSolverContext::linearSystem` 跨帧保存。模式变化时可输入上一份 METIS permutation；factor nnz 超过最近 fresh ordering 的 1.2 倍会在下一次求解强制重排。维护这一逻辑时要同时检查结果、`symbolic_analyses`、phase 时间和 fill，不能只减少 phase 11 次数。
- oneMKL 的 TBB threading layer 不由简单的 `mkl_set_num_threads` 可靠限制；当前每个 phase 用 `tbb::global_control`，运行时默认 16 线程（本机 16C/32T 最优），32 线程无收益；显式 0 才使用 oneMKL 默认。线程数属于 benchmark 配置，必须随结果报告。
- 性能脚本不安装 vcpkg CHOLMOD，而是校验下载固定 SuiteSparse 源码、显式启用 GPL supernodal 并构建 oneMKL bundle，避免无用 OpenBLAS。默认 PARDISO/system fallback 本身不强制 GPL，但一旦设置 `CIPC_CHOLMOD_ROOT` 使用性能 bundle，不能把其结果或分发许可写成普通 CHOLMOD core。

## E. 已修复记录（保留用于理解历史代码与旧文档）

- **CHOLMOD 生命周期与重复符号分析**：`CholmodSolver` 已改为 RAII；模式变化时释放 factor，同模式刷新只做数值分解；`IPCSolverContext` 持有的 `NewtonLinearSystem` 现跨 Newton/κ/时间步复用实例。
- **METIS 看似找到但未生效**：vcpkg导出target名为`metis`，旧查找器只接受`METIS::METIS`。现已归一化常见target名并链接检查`cholmod_metis`；模块保留供实验，但后续A/B确认生产ordering使用AMD更快。
- **多配置依赖混链**：SuiteSparse fallback target 现分别设置 Debug/Release imported location；vcpkg 的 `lib` 与 `debug/lib` 采用隔离查找。BLAS/LAPACK 统一 provider，避免一半来自标准 BLAS、一半来自 OpenBLAS。
- **假场景参数/空透传层**：已用 `SimulationScene::{TwistingMat,TwistingMatSoft,ClothOverBunny,Bunny2}` 接通 switch；hard/soft twisting 复用同一几何与目标运动，后者是 `weight=100` 的软边界范例。删除无第二实现、只转发到 `solveIPCStep` 的 FEMIntegrator 层。Bunny2 严格按 GPU_IPC 参考使用两份 scale=0.2 网格。
- **边界条件在重构中丢失/混用索引**：`e3df45e` 保留了 hard functor，却遗漏旧版 soft functor 的目标、能量、梯度和 Hessian。现集中为 `BoundaryConditions.cpp/.h`：hard/soft 各有独立且唯一的顶点集合；hard `updateDirection` 恢复 `rest/stepId/double alpha/dt`，明确返回 `x_new=x−p` 的方向，并在穿透回退时始终从同一个时间步起点重算 fractional motion；soft `updateTarget` 每时间步冻结目标，统一加入 `½w‖x−t‖² / w(x−t) / wI`。临时 soft twisting-mat Taylor 回归的一阶相对误差为 `3.2e-12/4.8e-11`，曲率比为 `0.999997/0.999994`；默认 hard twisting-mat 5 步终态保持不变。
- **CWD 相对输出与硬编码导出**：已由 `RuntimePaths` 统一写到 `<repo>/Output/` 并自动建目录；表面/截图均默认关闭、由按键切换；表面导出只读引用 `SimulationModel`。
- **检查点越界/半恢复风险**：加载先校验完整数量再整体提交，计时状态只在网格恢复成功时恢复，保存精度提高到 17 位。
- **旧检查点覆盖初始场景**：`SimulationOptions::resumeCheckpoint/writeCheckpoints` 已改为默认 false；viewer 和普通 API 构建默认 fresh。只有显式 `--resume` 才读取 `Output/tempData`。
- **`output_tetrahedraMesh` round-trip off-by-one**：tet 节点引用已写成 1-based，并补齐 Gmsh 结束标记。
- **Gmsh CRLF 导致 Linux tet 拓扑损坏/假性非正定**：旧 `split(" ")` 在 Linux 保留行末 `\r` 并把它当作额外 token，导致每个 tetrahedron 丢失首节点、末节点变成 `-1+offset`；Windows 文本模式会移除 `\r`，所以曾表现为仅 WSL 的 PARDISO `-4` / CHOLMOD 非正定。现按 Gmsh 2.2 的 element type/tag/node 字段使用 stream 解析，映射并校验节点 ID，且先写临时容器再整体提交。cloth-bunny 双平台均恢复为 4074 surface faces，non-quadratic 的 PARDISO/CHOLMOD 均通过。
- **零 Newton 次数除零**：平均碰撞数在 `total_iter==0` 时返回 0。
- **位置式参数解析**：已改为 key/value map，校验未知、重复、缺失、非有限与物理范围；默认与文件路径共用 `updateMaterial`。
- **Newton threshold/line-search 伪接口**：阈值已接通，并改为在线性求解后检查本轮当前方向；收敛时不再把近零方向送入 CCD/line search。`armijoCoefficient` 固定为0，严格要求 `E_trial<E0`，并有64次失败上限；旧的 `1e−3·初始步长` 截止已删除。line search 内的机器精度耗尽分支仍作为安全兜底，但不再是正常收敛路径。
- **宽阶段影响物理步长**：原 `Self_CCD_ActiveSet` 直接保存后端候选，体素假阳性会改变 partial/CFL 分支；现由统一静态 AABB+dHat 生成 partial pair，Full CCD 也在 ACCD 前统一扫掠 AABB 过滤。PT 使用入口 α，EE 使用 PT 更新后的 α；两个后端应用相同的精确区间过滤。原 CFL 合并公式保留。
- **无界循环**：ACCD、κ 外层、动画边界和穿透 safeguard 均加入上限、非有限/underflow 处理与异常诊断。
- **全局运行状态**：帧号、累计计时、碰撞和 checkpoint 状态已移入每个 Simulator 的 `IPCSolverContext`。
- **线性求解器结构**：矩阵装配、RHS、解向量和分派已移入`NewtonLinearSystem`；`LinearSolverOptions`首选PARDISO，否则使用优化CHOLMOD；Eigen-CG保留显式入口。旧自定义PCG、多RHS及`mesh3D::Constraints`已删除。
- **稀疏装配/临时分配**：只生成下三角有效 triplet；`NewtonWorkspace` 复用 gradient/BHessian/mutex，`IPCSolverContext` 持有的 `NewtonLinearSystem` 跨 Newton/κ/时间步复用 triplet/SparseMatrix/RHS/solver。
- **PARDISO 与大场景验证**：已接入 phase 指标、线程控制与 adaptive permutation，并以双 bunny2 做 1/5/50 时间步测试；50 步首次接触在 frame 43，末帧 ground/self=441/1006 且最小平方距离保持正值。完整数据见 `10_pardiso_report.md`。
- **重复 FEM/弯曲计算**：PFPX 初始化预计算，活跃 FEM 使用固定尺寸 Eigen；二次弯曲预计算 `Q⊗I₃` 并移除每轮 eigendecomposition。
- **非 quadratic hinge 不完整**：现以 `plateRigidity=Et³/[12(1−ν²)]` 为统一材料系数，预计算 `θ0` 与 `l0/(h0+h1)`；能量/梯度/Hessian 使用同一公式。CRLF tet 解析修复后已完成 Windows/WSL、PARDISO/CHOLMOD 单步对照；维护时仍必须同时构建 quadratic ON/OFF 并跑 cloth-bunny smoke。
- **摩擦与场景职责混杂**：摩擦冻结量、能量、梯度和解析 PSD Hessian 已集中到 `Friction.cpp/.h`；场景辅助函数改为用途命名并限制在 `Simulator.cpp` 内部作用域，删除未调用的 `buildSpecialPoints`。
- **SpatialHash 重建分配**：复用 voxel buckets 与 primitive occupancy，主候选使用线程局部排序 vector。
- **GPU 风格 CPU LBVH**：已按 `GPU_IPC/mlbvh` 的 Morton/Karras/face+edge/swept AABB 流程实现；现为默认 broad phase，SpatialHash 保留 A/B。
- **单体构建/GUI 绑定**：拆成 `cipc_core`、可选 viewer 和 headless；headless-only 配置无需 OpenGL，测试源码与 CTest 已按项目要求移除。
- **程序化网格索引**：`InitMesh` 已按 local corner 写 tet，`load_test(1/2)` 使用正确的 0/+8 offset，并补 boundaryTypes。
- **viewer 近裁剪面**：初始化与 resize 均统一为 0.1。
- **文件与特效遗留**：删除 shader 目录、GLEW、无效 VBO 分支、未消费的 viewer 开关、`fem_parameters.h`、2D/肌肉/肌腱加载器和旧 EKF 状态；viewer 改名 `ViewerMain.cpp`。
- **职责化命名**：`IPCSolver`、`BoundaryConditions`、`ContactMechanics`、`CollisionBroadPhase`、`AdditiveCCD`、`Elasticity`、`ElasticityMath`、`CholmodSolver`、`FeasibleStep`、`StageTimer`、`EncodedContact` 均按实际职责命名。
- **SVD 扁平化**：删除 `SVD/` 子目录和三个未调用的 JacobiSVD 包装；公开接口为 `RotationAwareSVD`，底层自包含 `ImplicitQR3x3SVD.h`。
