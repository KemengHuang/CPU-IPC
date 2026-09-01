# 08 — 重构与性能优化路线图

本文件记录后续优化的优先级、证据、风险与验收方式。除非已有基准数据，下面的收益判断均是**静态代码审查结论**，不是实测加速比。

## 实施状态（当前工作树）

- P0 已完成：`cipc_core/viewer/headless` 分层、逐帧 CSV、历史轨迹/LBVH 审计和 benchmark 脚本；独立测试源码与 CTest 后按项目要求移除。
- P1 已完成：下三角无占位零装配、固定尺寸活跃 FEM、PFPX 预计算、二次弯曲常量 Hessian、Newton/Energy workspace。
- P2 已完成安全部分：SpatialHash bucket/occupancy 复用、排序 vector scratch、GPU 风格 CPU LBVH；梯度去锁和彻底 CSR SpatialHash 尚未实施，因为当前 assembly 已不是主瓶颈。
- P3 部分完成：运行状态、构建 target、摩擦、线性系统、viewer 与核心文件职责均已分层；旧 PCG/特效/无调用文件已清理。`MMCVID` 已改名 `EncodedContact`，但负数槽位协议尚未类型化；`mesh3D` 其余职责仍待拆分。
- P4 完成低风险接口：CHOLMOD 仍是默认，METIS 以经过依赖校验的 cost-aware fallback 启用，Eigen-CG 作为共享装配后的可选后端保留并有 smoke test；未改 Hessian scale/近似或 Newton 算法，也未宣称 Eigen-CG 性能更优。

实测详见 `09_optimization_report.md`。

## 0. 当前基线

- 默认场景：`ClothOverBunny`，当前布料 `plane1024.obj` 为 1024 顶点 / 1922 三角形，bunny 为 7356 四面体。
- 当前构建只生成产品 target，不注册 CTest。维护验证采用 quadratic ON/OFF Release 构建、cloth-bunny/twisting-mat headless smoke、双 broad-phase 对照和 benchmark。
- `metrics.csv` 逐帧记录五阶段、nnz、活动集、回退与 CHOLMOD analyze/factorize 次数。
- CHOLMOD 已在每个 `solveBarrierSubproblem` 的 `NewtonLinearSystem` 内复用；同一 CSC 模式会复用符号分析。
- 已有无窗口多步入口与容差回归；并行结果仍非 bitwise 确定。

## P0 — 先建立可测、可回归的核心（最高优先级）

### 1. 拆出 headless runner

把“仿真步进”从 GLUT `display()` 解耦，新增独立命令行入口，例如：

```text
cipc_headless --scene cloth-bunny --steps 20 --broad-phase lbvh --output <dir>
```

最低要求：

- `FEMSimulator::step()` 不依赖 OpenGL/GLUT；viewer 只负责显示和输入。
- 明确控制场景、步数、是否恢复、输出目录和随机/确定模式。
- 每帧输出能量、最小距离、κ、Newton 次数、活动集数、线搜索回退数、矩阵 nnz 与 `time0..4`。

这是后续优化的前置条件：当前“每显示帧一步”会把刷新率、窗口事件和仿真性能混在一起。

### 2. 建立数值回归

为两个场景保存短轨迹基线，至少比较：

- 最终与逐帧顶点位置（绝对/相对容差）；
- 速度、总能量、最小约束距离、是否穿透；
- Newton/κ 外循环次数和失败状态；
- deterministic 模式下的可复现 hash。

浮点并行路径不应要求 bitwise 相等；先规定可接受容差，再谈优化。

### 3. 结构化 benchmark

保留现有 `time0..4`，但改为 CSV/JSON 每帧记录，并用固定 warm-up/测量帧。至少区分：

- FEM/障碍/摩擦梯度与 Hessian；
- triplet 生成、`setFromTriplets`、CHOLMOD analyze/factorize/solve；
- SpatialHash build/query、部分 CCD、全 CCD；
- 线搜索的能量计算、活动集重建和回退次数。

## P1 — 低风险、高概率收益

### 1. 只装配 CHOLMOD 使用的下三角

**状态：已完成。** `BHessian::toTriplets` 现在只输出下三角、自由 DOF 和非零数值；CHOLMOD（`stype=-1`）与 Eigen-CG 共用该装配。优化前默认 bunny 的 7356 个 tet 单独会生成：

```text
7356 × 144 = 1,059,264 triplets
```

对称 12×12 下三角只需 78 项，即 573,768 项，尚未计算布料、弯曲和接触就静态减少约 45.8% 的 tet triplet 数。实现采用两遍 count/prefix/fill，并记录最终 nnz。

直接写 Eigen 内部 CSC 与固定 superset pattern 仍未做，继续列为需单独验证的后续项。

### 2. 固定尺寸 Eigen 内核

**状态：活跃 FEM 路径已完成。** tet/triangle 循环已使用明确的 9×12、9×9、12×12、6×9、6×6 固定尺寸类型，并配合 `.noalias()` 与固定 segment；未使用的实验材料内核未做机械重写。

- `Matrix<double,9,12>`、`Matrix<double,12,12>`、`Matrix<double,6,9>` 等；
- `.noalias()` 和固定尺寸 segment/block；
- 函数签名返回/接收明确的固定尺寸类型。

后续若继续迁移休眠内核，仍应逐个做能量、梯度和 Hessian 数值核对，避免一次性机械重写。

### 3. 避免每轮重复计算静止几何算子

**状态：已完成。** `initMesh3D` 预计算 tet 9×12 与 triangle 6×9 PFPX，Newton 内核只读缓存；当前 bunny 的额外存储约 6.4 MB。该选择已纳入现有性能报告，尚未实现“只存 Dm⁻¹、现场紧凑计算”的替代 A/B。

### 4. 二次弯曲常量化

**状态：已完成。** `prepareQuadraticBending` 预计算 `hessianBase = Q⊗I₃`；Newton 只做 `plateRigidity·dt²` 缩放和固定 matvec。非 quadratic 路径也预计算 `θ0` 与 `l0/(h0+h1)`，并由有限差分测试覆盖。

### 5. 复用 Newton 工作区

**状态：已完成主要部分。** `NewtonWorkspace` 跨迭代持有 gradient、BHessian 和顶点 mutex；`NewtonLinearSystem` 持有 triplet、SparseMatrix、RHS/result 与后端缓存；`EnergyWorkspace`/`Friction::EnergyWorkspace` 复用线搜索约束向量。容器以 clear/resize 保留 high-water capacity，未改变公式。

## P2 — 碰撞与并行热点

### 1. 重构 SpatialHash 存储

静态和扫掠 `build` 有大量重复代码，并创建 `vector<vector<int>>` 后串行合并进 `unordered_map<int, vector<int>>`。建议分两步：

1. 先抽取共同的 bounds、voxel range 和 primitive emission 逻辑，保持数据结构不变；
2. 再评估扁平 CSR 风格 occupancy（count → prefix sum → fill），复用缓冲并并行填充。

第二步可能显著减少小 vector 分配和串行 merge，但必须测量不同网格密度下的占用分布。

### 2. 去掉查询中的临时 `unordered_set`（已完成主路径）

PT/EE、Full CCD 和 edge-triangle safeguard 已使用排序 vector；PT/EE scratch 为 TBB 线程局部并复用容量。generation-stamp 未采用，因为当前 vector 方案已通过 profile 且更容易保持线程安全。

- 线程局部 `vector<int>` + sort/unique；或
- generation-stamp 数组（`visited[id] == generation`）去重。

stamp 法速度潜力高，但需线程私有或分片，避免数据竞争与 O(primitive_count × threads) 过大。

### 3. 梯度累加去锁化

tet、triangle、弯曲对每个顶点使用 `spin_mutex`，高价顶点会产生争用。候选方案：

- 预计算单元图着色，同色单元并行、颜色间串行；
- `enumerable_thread_specific` 局部梯度后归并；
- 按顶点邻接表 gather，而非单元 scatter。

三者的内存/局部性差异很大，应通过默认场景和大布料场景对比，不能凭直觉替换。

### 4. 减少线搜索重复工作

每次回退都会更新位置、重建 SpatialHash/活动集并重新计算完整能量。可分离并缓存：

- 静止常数与不随 trial step 变化的量；
- 仅在需要时做昂贵的事后相交检查；
- 统计回退原因（能量上升/地面穿透/自穿透），再决定优化哪条路径。

profile 证明反复 broadphase build 是热点后，已先完成 SpatialHash bucket 复用，再实现 GPU 风格 CPU LBVH；增量 refit 尚未做。

### 5. TBB 粒度与小任务

当前大量 `parallel_for(..., grainsize=1)`。对小活动集，调度开销可能超过计算；应基于工作量设置 blocked range 粒度，并为小 N 提供串行阈值。

## P3 — 代码结构重构

### 1. 拆分 `mesh3D` god object

建议按生命周期拆为：

- `MeshTopology`：tet/triangle/surface/edge 与邻接；
- `RestGeometry`：rest position、Dm inverse、体积/面积、弯曲算子；
- `MaterialParameters` 与 `SolverSettings`；
- `SimulationState`：x、xPrev、x̃、v、边界；
- `ContactState`：活动集、摩擦滞后量；
- `NewtonWorkspace`：临时梯度/Hessian/solver buffers。

目标不是追求类数量，而是让 const/static 数据与每轮变化数据边界清晰，为缓存和测试创造条件。

### 2. 将 `EncodedContact` 升级为 tagged contact

用显式 tagged contact 类型（PP/PE/PT/EE + 顶点 ids + multiplicity + mollifier metadata）替代 `[0]<0`、`[3]=-count` 等隐式协议。先提供集中 encode/decode 适配层，再迁移调用点，避免一次性改坏能量/梯度/Hessian的一致性。

### 3. 移除全局求解状态

**状态：求解侧已完成。** 帧号、计时、碰撞与恢复状态已进入 `IPCSolverContext`；viewer 状态限制在 `ViewerMain.cpp` 匿名 namespace。若未来支持多个窗口，再提取 `ViewerState`。

### 4. 配置与场景数据化

- 参数文件改为键值解析，拒绝未知/重复/缺失键并做范围校验；
- 场景资产、变换、边界条件、地面和输出策略放入 scene config；
- CLI 覆盖配置值，并把最终生效配置写入结果目录以便复现。

### 5. 构建目标分层

已将原单一 glob executable 拆成：

```text
cipc_core      FEM + IPC + collision + solver（无 OpenGL）
cipc           GLUT/OpenGL viewer
cipc_headless  benchmark/regression CLI
```

核心源文件已显式列入 `cipc_core`；viewer/headless 为独立产品 target，不再生成测试 target。后续仍需删除或移入 `experimental/` 的大段实验代码。

## P4 — 高风险算法级优化（最后考虑）

- Modified Newton / lagged Hessian：若多轮 Hessian 变化小，可少做装配或分解，但会改变收敛行为。
- Eigen-CG + Incomplete Cholesky 已作为可选后端接通；大网格上是否胜过 CHOLMOD、接触刚度下是否足够鲁棒仍需专门 benchmark，默认不变。
- 固定 superset sparsity pattern：可长期复用符号分解，但可能显著增加 fill-in 和内存。
- 接触 Hessian 的解析 PSD/Gauss-Newton 近似：可减少逐接触特征分解，但属于算法变更，必须核对 IPC 收敛与无穿透性质。
- LBVH 已完成并默认启用；更高风险的增量 refit 或 SAH rebuild 仍需更大场景 profile。

## 推荐实施顺序

1. 已完成：headless/回归/benchmark、下三角、workspace、固定尺寸 FEM、PFPX/弯曲静态化。
2. 已完成：SpatialHash 复用与 GPU 风格 CPU LBVH；继续按场景 A/B，默认 LBVH。
3. 下一步：细分 `linear_ms`（triplet / setFromTriplets / analyze / factorize / solve）并评估直接 CSC 数值装配。
4. 下一步：拆 `mesh3D`、把 `EncodedContact` 改为 tagged contact、保存 checkpoint 拓扑/参数 hash，逐步清理仍嵌在活跃大文件中的实验函数。
5. profile 证明 assembly 锁争用成为瓶颈时，再比较 graph coloring / TLS / gather。
6. 只有更大场景证明 CHOLMOD 是扩展瓶颈时，系统比较现有 Eigen-CG 的容差/预条件与 CHOLMOD；lagged Hessian 仍作为独立高风险实验。

## 每项优化的完成标准

- 文档说明行为/API/数据结构变化；`07_gotchas.md` 状态同步。
- quadratic ON/OFF Release 构建无新增警告，代表场景 headless smoke 正常。
- 固定场景数值回归通过，最小距离不恶化为穿透。
- 报告同一机器、同一配置下的中位数与波动，不只报告单次最快值。
- 同时报告时间、峰值内存、Newton/线搜索次数；算法迭代数变化时不能只比较 wall time。
