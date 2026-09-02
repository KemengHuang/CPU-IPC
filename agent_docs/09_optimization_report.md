# 09 — 全方位优化实施报告

本报告记录作为 CPU IPC benchmark/reference 的当前工作树优化、同机实测数据和验证边界。测试平台为本开发机 Windows/MSVC Release；不同机器只能比较趋势，不能直接复用绝对时间。

## 1. 结构与可测性

- 拆分产品 targets：`cipc_core`（无 OpenGL）、可选 `cipc` viewer、`cipc_headless`。测试源码与 CTest 后按项目要求全部移除。
- 删除只有一个实现且仅转发调用的 FEMIntegrator 层，`FEMSimulator::simulateStep` 直接进入 `solveIPCStep`。
- `CIPC_BUILD_VIEWER=OFF` 的独立 configure/build 已验证，不查找 OpenGL/GLUT；GLEW/shader 依赖已从项目删除。
- `SimulationOptions` 控制场景、恢复、运行文件、checkpoint、verbosity、broad phase 和线性求解后端。
- `IPCSolverContext` 替代旧求解文件的帧号/累计计时/碰撞全局变量，可在同一进程创建多个 Simulator。
- `IPCStepStats + metrics.csv`：逐帧记录五阶段耗时、Newton/κ、总/能量/穿透回退、mean/min/max α、碰撞、最小平方距离、活动集、nnz、analyze/factorize；PARDISO 追加 phase 11/22/33、线程数和 factor nnz。
- `scripts/benchmark.py`：独立进程重复运行，报告 median/min/max，并保留每次 metrics。
- 摩擦冻结量、能量、梯度和解析 PSD Hessian 拆到 `Friction.cpp/.h`；`IPCSolver.cpp` 只保留调用与总能量编排。
- 边界条件拆到 `BoundaryConditions.cpp/.h`：Dirichlet 与 soft-target 顶点不再共享易错索引数组；恢复重构时遗漏的 soft 目标更新及 `½w‖x−t‖² / w(x−t) / wI` 一致装配，恢复 hard callback 的 rest/step/alpha/dt 上下文，并用命名 predicate 取代 0/1/2 magic comparisons。
- Newton稀疏装配、RHS/scatter和后端分派拆到`NewtonLinearSystem.cpp/.h`。默认PARDISO，不可用时自动选择优化CHOLMOD；Eigen-CG保留显式入口，三个后端共用同一矩阵/边界条件装配。
- 删除旧自定义 PCG、多 RHS 实验入口、无用 `mesh3D::Constraints`、`FEMTimeIntegrator` 透传层和未调用的场景辅助函数；`Simulator.cpp` 辅助函数采用用途命名并限制为内部链接。
- 文件职责整理：`ViewerMain/IPCSolver/ContactMechanics/CollisionBroadPhase/AdditiveCCD/Elasticity/ElasticityMath/CholmodSolver/FeasibleStep/StageTimer` 取代含糊旧名；`EncodedContact` 取代 `MMCVID` 类型名，`SimulationModel` 取代 `model_tet`。
- 删除 shader/GLEW、无效 viewer 特效开关、`fem_parameters.h`、2D/肌肉/肌腱加载器、空 model loader、EKF 状态、未调用 ACCD broad-phase 备份与无用数学函数；Python 视频工具移入 `scripts/` 并改为参数化入口。
- SVD 扁平为 `RotationAwareSVD.cpp/.h + ImplicitQR3x3SVD.h`；数学辅助直接内联在底层头文件，三个未调用 JacobiSVD 包装删除。tet 防翻按要求独立保留为 `TetInversionGuard`，默认不接入主流程。
- 弯曲材料统一为 `plateRigidity=Et³/[12(1−ν²)]`。非 quadratic 路径拆为 `HingeBending`，预计算 `θ0` 与 `l0/(h0+h1)`；能量、精确梯度/Hessian 共用同一公式，并在实现阶段完成有限差分核对。

## 2. 稀疏求解与装配

- `CholmodSolver` RAII、模式变化检测、同模式 symbolic reuse、失败检查和 analyze/factorize 计数；`cholmod_solve2` 复用 solution/Y/E dense workspace。
- `PardisoSolver` 使用 oneMKL SPD `mtype=2`、phase 11/22/33、TBB 线程上限和自适应 METIS permutation；直接把共享 lower CSC 当作等价 upper CSR，避免额外转置。模式不变时复用全部 symbolic/factor workspace，模式变化时复用 permutation，fill 超过 fresh 参考 1.2× 则下一次强制重排。
- CHOLMOD性能路径升级为GPL supernodal+oneMKL LP64/TBB共享库；CMake通过`CIPC_CHOLMOD_ROOT`选择并验证supernodal/METIS符号，运行时采用实测AMD ordering和按nnz自适应4/8线程。普通system CHOLMOD仍作为兼容路径。
- `CIPC_ENABLE_METIS_ORDERING=ON` 保证CHOLMOD Partition/METIS模块可用于实验并链接检查`cholmod_metis`；最终生产ordering根据新A/B固定AMD，因为它快于强制METIS与multi-method AUTO。
- 修复 vcpkg 只导出无命名空间 `metis` 导致 METIS 未接入的问题，并兼容 Linux 的 `metis.h + libmetis` 及内嵌 Partition 的 CHOLMOD；SuiteSparse fallback imported targets 现在区分 Debug/Release，BLAS/LAPACK 固定为同一 provider，同时修复 `SuiteSparse_SPQR_FOUND` 拼写。
- `BHessian::toTriplets` 改为两遍并行 count/prefix/fill：仅自由 DOF、非零数值、全局下三角。
- triplet、SparseMatrix、RHS/result和三个solver workspace全部由`IPCSolverContext`中的`NewtonLinearSystem`跨Newton、κ子问题和时间步复用。

METIS接入时做了策略诊断，而不是把“调用METIS”直接当成加速。早期cloth-bunny已显示AMD不慢于METIS；新supernodal+MKL bunny2 A/B进一步得到AMD 2.464s、METIS 2.557s、multi-method AUTO约2.66s，因此生产策略固定AMD，同时保留Partition模块供后续矩阵实验。

### oneMKL PARDISO 与 bunny2

新增`SimulationScene::Bunny2`，严格参考GPU_IPC使用两份`bunny2.msh`、scale=0.2和`YoungModulus=1e5`。合并后为38,386顶点、159,870 tet、115,158 DOF；首步Hessian `nnz=2,202,090`。当前同版本中位数为PARDISO 1.301s、优化CHOLMOD 1.625s，处于同一性能量级；完整数据见`10/11`专项报告。

50 时间步接触审计中，首次 ground/self contact 出现在 frame 43；末帧为 441/1006，`min_distance2=7.3483170186983491e−7>0`。adaptive permutation 总时间 59.167 s，比每次 fresh METIS 的 73.033 s 快 19.0%，比永远固定首份 permutation 的 62.838 s 快 5.8%，并保持接触整数指标一致、终态仅有浮点舍入差。完整线程表、phase、内存、5 步中型场景和排序策略数据见 `10_pardiso_report.md`。

以下矩阵与性能数字是早期 `plane100` 场景的历史基线；当前默认场景已由用户改为 `plane1024`，不可直接横向比较：

| 指标 | 优化前 | 优化后 | 变化 |
|---|---:|---:|---:|
| `matrix_nnz` | 211,941 | 107,574 | −49.2% |
| `linear_ms`（同一代表运行） | 489.8 ms | 376.7 ms | −23.1% |

`linear_ms` 仍包含 triplet、Eigen sparse build、analyze/factorize/solve；下一步需继续细分，不能把全部时间归因于 CHOLMOD 数值分解。

## 3. FEM 与 Newton workspace

- tet/triangle PFPX 初始化预计算，活跃路径使用固定 9×12 / 6×9 类型。
- PK1 使用 Eigen Map 读取列主序 9/6 向量，避免动态 `vec_double` 临时量。
- tet/cloth Hessian 使用固定尺寸矩阵与 `noalias`。
- quadratic bending 预计算 `hessianBase=Q⊗I₃`；hinge bending 预计算静止角和无量纲几何权重；移除每 Newton 的静止几何重复计算。
- BHessian、gradient、mutex、triplet 和线性求解缓冲按 high-water capacity 复用。
- `Newton_Solver_Threshold` 已接入原有无穷范数判据，默认值保持旧行为。

代表运行的 `assembly_ms` 从约 46.1 ms 降至约 19.9 ms（−56.9%）。该项只占总时间一部分，所以不能直接等同总加速。

## 4. SpatialHash 优化

- 不再每次 `voxel.clear()` 销毁所有 bucket vector；只清空上一轮 active bucket 并复用容量，空 key 累积超过阈值时重建 map。
- surface local index、静态/扫掠 voxel ranges、edge/face occupancy、point/edge occupancy 都成为可复用成员缓冲。
- PT/EE/full-CCD scratch 使用 TBB 线程局部排序 vector，替代逐 query `unordered_set`。
- edge-triangle safeguard 也使用排序 vector。

代表运行：

| 阶段 | 优化前 | 优化后 |
|---|---:|---:|
| `ccd_ms` | 153.8 ms | 53.6 ms |
| `line_search_ms` | 252.1 ms | 71.0 ms |

线搜索下降主要来自反复 rebuild 时复用 voxel bucket/occupancy，而非改变能量或 ACCD 公式。

## 5. 参考 GPU_IPC 的 CPU LBVH

参考：`GPU_IPC/GPU_IPC/mlbvh.cu` 与 `mlbvh.cuh`。

CPU 实现 `LBVH.cpp/.h` 保留：

- face/edge 两棵独立树；
- 静态或 `x → x−αp` 扫掠 primitive AABB；
- scene normalization、30-bit Morton、`(morton<<32)|primitiveId` 唯一 key；
- Karras `determineRange/findSplit` 构造 `2N−1` 线性树；
- postorder internal AABB 与显式栈查询；
- PT/EE 共享顶点、边对顺序和全驱动过滤继续由现有 IPC 层处理。

CPU 与 GPU 的差异：CPU overlap 包含边界接触；Full CCD 采用 CPU 既有 thickness=0 的零 padding，GPU 源码查询额外使用 `sqrt(dHat)` 候选 padding。

验证：

- 历史 LBVH 审计覆盖 257 个随机盒、重复 Morton 中心、退化 z 轴、单叶与空树，并与 brute force 对比。
- 历史 headless 审计覆盖初始活动集、partial CCD pair、不同 PT/EE α 的扫掠 pair 和双后端单步轨迹。相关测试源码现已移除。
- 当前顺序 α 下 cloth-bunny 双后端 20 步：逐帧所有 Newton/回退/接触/nnz 等整数指标完全相同；mean/min α 最大差约 `4.7e−11`，最终轨迹汇总差约 `1e−12`。

性能 A/B（3 次独立运行，20 步总和中位数；这是恢复 PT→EE 顺序 α 之前的固定-sweep 历史数据）：

| 后端 | 中位总时间 | 相对 SpatialHash |
|---|---:|---:|
| 优化后的 SpatialHash | 4388.2 ms | 基准 |
| CPU LBVH | 3771.7 ms | −14.0%（约 1.16×） |

两后端现使用相同 partial/full CCD pair 与 CFL 合并结果，因此 Newton/回退/α 一致；性能差只来自宽阶段成本。LBVH 设为默认，SpatialHash 保留 `--broad-phase spatial-hash`。

## 6. 总体性能

本节的绝对性能数据采自恢复 PT→EE 顺序 α 之前的固定-sweep 版本。顺序 α 恢复属于正确性/语义调整，本轮只做数值回归，没有在隔离环境重新跑正式性能基准，因此不把并行验证运行的耗时与本节数字直接比较。

cloth-bunny、1 步、独立进程：

| 版本 | 重复数 | 中位 step time |
|---|---:|---:|
| 本轮优化前 SpatialHash | 3 | 952.0 ms |
| 最终 SpatialHash（统一 CCD pair/CFL） | 5 | 748.1 ms |
| 最终默认 LBVH（同物理步长） | 5 | 621.1 ms |

相对最初版本，最终默认 LBVH 单步中位数降低约 **34.8%**，约 **1.53×**。最终数值比早期“优化后 SpatialHash 524.7 ms”更保守：后者仍受 backend-dependent partial CCD 分支影响，不能作为最终正确版本的性能数值。

优化后 SpatialHash 20 步进程峰值 working set 实测约 70.7 MiB。该数值没有同条件“优化前”内存基线，只能作为当前运行记录，不能宣称内存下降比例。

## 7. 健壮性与配置

- 参数文件按 key 解析；未知/重复/缺失、NaN/Inf、非正材料/dt/dHat、非法 ν/collision flag 均拒绝。
- 默认设置与文件设置共用材料换算。
- ACCD 10000、κ 64、穿透/动画 safeguard 64 次上限；非有限推进、step underflow 和超限均有明确错误。
- 网格加载失败返回并由场景构建抛出异常，不再在活路径 `exit(-1)`。
- checkpoint 完整数量验证、17 位精度、恢复状态与 metrics 续写已做 10+1 步验证。

## 8. 自动验证结果

- WSL2 Ubuntu 22.04.5（GCC 11.4、CMake 4.3.3、Ninja 1.10.1）已用零参数 `./build.sh` 从系统检查、依赖安装到产品链接端到端通过；Linux oneMKL 2025.2、METIS/TBB 来自固定 `x64-linux` vcpkg，SuiteSparse 7.14/CHOLMOD 5.3.5 由校验下载的源码直接构建，`ldd` 中没有 OpenBLAS。`--dependencies-only --no-system-packages` 与 bundle stamp 增量复用分支均通过。默认 PARDISO 与显式 CHOLMOD 的 twisting-mat 单步 smoke 成功且终态仅有浮点归约误差；`./build.sh --viewer` 也已完成 OpenGL/X11/FreeGLUT 配置与链接，并通过 WSLg 3 秒启动 smoke。
- Windows `.\build.cmd -HeadlessOnly` 已从空的固定 revision vcpkg checkout 验证 MSVC/CMake 自动解析、`tbb[core]`/METIS/oneMKL 安装、无 OpenBLAS 的 SuiteSparse bundle、六个 DLL 自动部署和 CPU-IPC Release；`.\build.cmd -DependenciesOnly -HeadlessOnly` 及 bundle stamp 二次复用也通过。MSVC Debug 自动 CHOLMOD smoke 通过。
- 边界重构的临时 soft twisting-mat 验证实际执行 2 次 Newton、0 次能量回退；两次 Taylor 诊断一阶相对误差 `3.2e-12/4.8e-11`、曲率比 `0.999997/0.999994`。恢复默认 hard 配置后，Windows/WSL、PARDISO/CHOLMOD 的 5 步 `Newton=2`、`nnz=151335` 均不变，`squared_norm_sum=560.775874160993xx`，与重构前差异不超过约 `1.2e-13`。
- 安装器开发期间发现 WSL non-quadratic 曾在 PARDISO phase 22 报 `-4`、CHOLMOD 报非正定；根因不是 Hessian/solver，而是 Gmsh 行末“空格+CRLF”被旧 `split(" ")` 在 Linux 解析出额外 `\r` token，破坏了全部 tet 连接。改为结构化 Gmsh 2.2 stream parser 后，Windows/WSL 均为 4074 surface faces、8 Newton、12 energy backtracks、`nnz=161979`，PARDISO/CHOLMOD 单步终态在浮点误差内一致。
- 新 bundle 的 Windows 五次独立进程中位数：bunny2 单步 PARDISO/CHOLMOD 为 `1.310/1.686 s`，cloth-bunny 五步为 `0.575/0.621 s`，twisting-mat 五步为 `0.431/0.412 s`，终态逐次一致。同一时段旧/新 CHOLMOD bunny2 对照为 `1.670/1.686 s`（约 1%）；不把安装简化宣称为算法提速，PARDISO 仍是默认。
- MSVC Release `--clean-first` 全量构建：通过，无新增编译警告。
- headless-only（viewer OFF）独立构建：通过。
- METIS ordering ON/OFF 两套 MSVC Release 构建与 twisting-mat/cloth-bunny CHOLMOD smoke 均通过；ON 下 Eigen-CG smoke 也通过。ON 配置会实际编译/链接检查 `cholmod_metis` 后才提供 `Partition`，不会只根据包名静默假定支持。
- 多配置映射已在生成的 VS 工程中核对：Debug 链接 `debug/lib`，Release/RelWithDebInfo/MinSizeRel 链接 `lib`；Debug headless 构建和 twisting-mat CHOLMOD smoke 通过。
- PARDISO/优化CHOLMOD逐帧结果在浮点归约误差内一致；twisting-mat 5步和Eigen-CG单步smoke通过。
- Windows MSVC 与 Linux GCC 的 quadratic ON/OFF Release 配置及 cloth-bunny 代表 smoke 均通过；项目不再运行 CTest。
- oneMKL PARDISO Release、非 quadratic+PARDISO、`CIPC_ENABLE_PARDISO=OFF` 与普通 MSVC Debug 四条构建路径均通过；Windows 静态 oneMKL 因 CRT 边界在 Debug 自动排除，显式选择时按预期清晰报错。
- bunny2 PARDISO 1/5/50时间步运行通过；PARDISO/优化CHOLMOD单步终态在浮点舍入内一致。
- 默认 LBVH 与 SpatialHash、TwistingMat 三条 1 步轨迹均纳入回归。
- SpatialHash 与 LBVH 各完成 20 步运行；最小平方距离保持正值。
- viewer 去除 GLEW/shader 后完成隐藏窗口启动冒烟：成功创建并持续运行 3 秒后由验证流程主动关闭。
- `scripts/generate_video.py` 与 `scripts/merge_videos.py` 通过 AST 语法检查。
- 场景恢复审计：仓库旧 checkpoint 的 step=0 坐标和与 fresh 场景显著不同；现默认禁用恢复/写 checkpoint。当前 `plane1024` 的 cloth 顶点 0 经 X 轴 π/2 旋转后为 `(-1.414106,-0.017452,0.000108)`，追加在索引 1024 的 bunny 首顶点为 `(-0.229744,-0.5465085,-0.1556705)`。

### Line-search / CFL 修复

当前 `plane1024 + plateRigidity∝t³` 产品 smoke 中两后端轨迹一致：8 Newton、12 次能量回退、0 次穿透回退、`min_distance2=1.485724192802e−5`、`matrix_nnz=145053`。下面其余多步/scale 数字来自更早的 `plane100` 历史快照，只用于记录当时诊断结论。

bunny2 无接触初始阶段暴露了另一类“机器精度长尾”：旧顺序在求解当前方向前检查上一轮方向，导致已经达到 `gTp≈7.33e−30` 的新方向仍进入严格线搜索；TBB 能量归约仅在约 `1e−16` 末位波动，单帧会随机二分 0～52 次。5 次独立运行×前5帧共417次回退，5步中位 4429.220 ms。现改为先求本轮方向、立即按既有阈值判断，再决定是否进入 CCD/line search：同样25帧回退全部为0，5次终态完全一致，5步中位降为3396.438 ms（−23.3%）。`armijoCoefficient=0` 与严格接受条件未改变；每帧多出的最终 numeric factorization 是当前方向的收敛确认。

保持原 `beta` 辅助退出的隔离50步运行仍无穿透：总时间54.716 s、153个接受Newton步、15次能量回退；末帧ground/self=441/999、`min_distance2=7.3487834099833287e−7>0`、`nnz=2,214,492`。相比旧停止点的59.167 s/151 Newton/202回退，停止点和轨迹有小幅变化，这是按当前方向判断的预期语义变化，不应与solver误差混淆。

- `plane100` 恢复 PT→EE 顺序 α 后，首帧 LBVH/SpatialHash 均为 15 Newton、19 次能量回退、0 次穿透回退；此前固定-sweep 快照为 17 Newton/19 回退。
- 原 line search 在 α 缩到初值的 `1e−3` 后直接退出，可能接受仍不下降的试步；现按要求固定 `armijoCoefficient=0` 并严格要求 `E_trial<E0`，最多64次，否则报错，不使用允许能量上升的容差。
- 回退拆为 `energy_backtracks` 与 `intersection_backtracks`；固定-sweep 版本曾做 100 步审计（1144 Newton、954 能量回退、穿透回退 0），该数字保留为历史稳定性记录，不代表当前顺序 α 的 100 步统计。
- 已保留原 partial/CFL/full 合并公式。真正的后端差异来自 partial set 直接采用 backend 原始候选，以及旧 LBVH 未像 SpatialHash 一样按传入 α 缩小 EE 有效候选；“EE 使用 PT 更新后的 α”本身并不是错误。现 static/swept AABB 二次过滤统一 pair，并恢复顺序语义：PT 更新 α 后，两个后端都用该更小区间精确过滤 EE；ACCD 返回值仍是完整搜索方向上的全局 α。

`plane100` 顺序 α 的 20 步 LBVH/SpatialHash 都是 94 Newton、78 次能量回退（平均 0.830/Newton）、穿透回退 0；两后端逐帧整数指标完全一致，最小 `min_distance2=1.471008086921e−5 > 0`。

下面的插值实验在此前固定-sweep 快照上完成，用于比较线搜索策略本身，不应当与当前顺序 α 的 78 次回退直接横向比较：

尝试了用失败能量与方向导数做 safeguarded quadratic interpolation：

| 策略 | 20步 Newton | 能量回退 | 最大/单Newton | 总时间 |
|---|---:|---:|---:|---:|
| 二分 0.5（固定-sweep 历史基线） | 94 | 65 | 5 | 3570 ms |
| 插值 clamp [0.1,0.5] | 101 | 50 | 4 | 3920 ms |
| 插值 clamp [0.25,0.5] | 101 | 62 | 4 | 3924 ms |
| 插值 clamp [0.4,0.5] | 97 | 55 | 3 | 3694 ms |

插值能让计数更好看，但会过早选择小步、增加 Newton，所有方案总耗时都更差，因此撤回。要继续降低接触建立阶段的 3–5 次回退，需要改善 barrier/friction 的 Newton Hessian 或 activation-aware 初始步长，而不是只改回溯倍率。

方向 Taylor 审计显示总 gradient 与能量中心差分相对误差通常约 `1e−9`，接近零方向导数时约 `1e−6`，没有发现一阶装配错误。固定-sweep 诊断快照中，关闭摩擦后首帧从 17 Newton/19 回退降为 7/1，确认频繁 globalization 主要由滞后摩擦主导。

随后将摩擦 PSD 处理从通用 eigendecomposition 改为解析切向特征构造：滑动区 `(0, μλ/r)`，C1静摩擦区 `(μλf2, μλf1)`，以 lifted parallel/perpendicular rank-one PSD 和装配。在当时的固定-sweep 快照中，scale=1 的轨迹、17 Newton/19回退与旧实现完全一致；解析式也曾与旧闭式及 PSD 性质做数值核对。

曾在固定-sweep 快照中将摩擦 Hessian 单独乘 scale（能量/gradient 不变）作为 modified-Newton 诊断：scale=2.25 在20步由94 Newton/65回退降到69/7，50步为429 Newton/121回退；scale=1 的50步回退更多。但该 scale 是场景经验阻尼，不是摩擦公式缺失系数。按项目要求，正式代码已删除 Hessian scale 接口并固定使用精确 scale=1；Barrier Hessian scale 扫描同样没有作为功能保留。

### LBVH 专项复核

- 历史 LBVH 重复审计：50/50 通过。
- 最终双后端/双场景/扫掠 pair 回归 `--repeat until-fail:5`：5/5 通过。
- LBVH cloth-bunny 100 步：正常完成；`metrics.csv` 恰有连续 frame 0..99，所有阶段/κ/距离字段均为有限值。
- CFL/pair 归一化版本曾完成100步审计：最小 `min_distance2 = 1.37719209285429e-5 > 0`、穿透回退0；后续严格能量与解析摩擦 Hessian曾由20/50步实验覆盖。
- 初始两后端的 active/mollified/partial pair 完全一致；合成非平凡位移的 swept PT/EE 有效 pair 也完全一致。

验证边界：通用树 query 已逐项对 brute-force AABB overlap，静态与合成扫掠有效 pair 已做双后端等价比较，并通过 100 步集成审计；但尚未对 100 步中每一帧都另做 O(P×F+E²) brute-force pair 导出。因此仍是强工程验证，不宣称形式化证明。

## 9. 未实施或尚未完成评估的高风险项

- Eigen-CG 已可选运行，但尚未在大网格上系统比较容差、预条件、收敛鲁棒性和性能；它不会被能力感知逻辑自动选择。
- CHOLMOD GPL supernodal + oneMKL已完成隔离benchmark；剩余风险是GPL二进制分发边界，而不是性能未知。
- PARDISO为Release默认；不可用或MSVC Debug时自动使用优化CHOLMOD。下一线性热点是triplet→Eigen CSC构造与symbolic阶段。
- 优化 CHOLMOD 历史上将 bunny2 从 supernodal+OpenBLAS 约 8.67s 降到约 1.63s（约 5.3×）；固定安装 bundle 的当前五次中位数为 CHOLMOD/PARDISO：bunny2 `1.686/1.310s`、cloth 五步 `0.621/0.575s`、twisting-mat 五步 `0.412/0.431s`。PARDISO 仍作为大场景默认，完整 A/B 见 `11_cholmod_mkl_report.md`。
- lagged/modified Newton、接触 Hessian 近似。
- 直接 CSC 数值装配与固定 superset sparsity（可能改变 fill-in/内存）。
- 梯度 graph coloring/TLS/gather（当前 assembly 已降为小占比）。
- 完整拆分 `mesh3D`，并把当前仍采用负数槽位协议的 `EncodedContact` 进一步升级为显式 tagged contact。

这些项目需要更大网格 profile 或更强回归后再推进，当前不把算法风险包装成“优化”。
