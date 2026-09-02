# 02 — 仿真主流程（每个时间步）

入口链：viewer `display()` 或 headless loop → `FEMSimulator::simulateStep` → `solveIPCStep`。viewer 每显示一帧 = 一个时间步，headless 按 `--steps` 显式循环。因此命令行的“一步”是完整 IPC 时间步/帧，不是一个 Newton step；一次时间步通常执行多轮 Newton、CCD、线搜索和线性分解。

## 1. `solveIPCStep(stepId, mesh, broadPhase, ground)` — 时间步驱动（`IPCSolver.cpp`）

按顺序：

1. **断点恢复**：每个 `IPCSolverContext` 最多检查一次 `<output>/tempData/timeCost.txt`，且只有位置/x̃ 检查点完整且顶点数一致时才恢复累计指标、帧号和 κ；不再使用文件级全局状态。
2. **κ 初始化**（2110-2114）：`upperBoundKappa` 封顶；`Kappa < 1e-16` 时 `suggestKappa` → `initKappa`（见 §5）。
3. **摩擦设置**（`USE_FRICTION`）：`Friction::initialize(mesh, gd)` —— 用**时间步开始时**的活动集计算每个接触的 λ、最近点坐标（MMDistCoord）、切空间基（MMTanBasis），快照 `*_activeSet_lastH`。
4. **动画硬约束**（2120-2159）：若 `update_hard_constraint_functor != nullptr`，对 `boundary_vertexes_indices` 算位移 `moveDir`，用 `Self_largestFeasibleStepSize_CCD(moveDir, slackness=0.8)` 求 `new_alpha`，施加 `stepForward(..., boundary_update=true)`，再 `new_alpha /= 2` 直到无穿透；重建碰撞集。
5. **κ 外层循环**：反复 `solveBarrierSubproblem`；每轮后评估全部地面 + 自碰撞 PT 约束值，当 `minCoeff < mesh.dTol` **或** `maxCoeff < mesh.Hhat`（或无约束）时跳出。外层现有 64 轮硬上限，超限抛出带语义的异常。
6. **速度更新**：`v = (x − V_prev)/dt`（`is_quasi_static` 则 v=0）；`V_prev = vertexes`；`updateInertialTarget`。
7. **日志/指标**：`IPCStepStats` 记录逐帧 stage ms、Newton/κ、总/能量/穿透回退、单 Newton 最大回退、超过 2 次的 Newton 数、mean/min/max α、碰撞、活动集、nnz、analyze/factorize；PARDISO 还记录 phase 11/22/33 累计时间、有效线程数与 factor nnz；写 `<output>/metrics.csv`。

## 2. `solveBarrierSubproblem` — Newton 循环

`iterationLimit = 10000`。每个 barrier 子问题创建一次 `NewtonWorkspace`，复用 gradient、BHessian 和 mutex；它引用 `IPCSolverContext` 持有的共享 `NewtonLinearSystem`。后者持有 triplet、SparseMatrix、RHS/result 与四个线性后端的缓存，并跨 κ 子问题和时间步存活；网格顶点数变化时才重建。每次迭代严格按序：

| 顺序 | 步骤 | 位置 | 计时器 |
|---|---|---|---|
| 1 | `computeGradientAndHessian(mesh, gradient, BH, gd)`，返回碰撞数 | `:1971` | time0 |
| 2 | `NewtonLinearSystem::solve`：统一装配后按配置调用 SuiteSparse LDL、CHOLMOD、PARDISO 或 Eigen-CG；三个直接法在稀疏结构不变时复用符号分析 | `NewtonLinearSystem.cpp` | time1 |
| 3 | 用**本轮刚求出的方向**做收敛检查：`directionInfinityNorm(moveDir)`；阈值 = `Newton_Solver_Threshold·sqrt(bboxDiagSize2)·dt`；低于阈值则记录本轮 assembly/linear 时间并在 CCD/line search 前 break | `solveBarrierSubproblem` | — |
| 4 | 可行步长上限（见 §3） | `:1990-2024` | time2 |
| 5 | `lineSearch(...)` 回溯 | `:2031` | time3 |
| 6 | `postLineSearch(mesh, gd, Kappa)` κ 自适应 + 近约束簿记 | `IPCSolver.cpp` | time4 |

`collisionNum` 来自 `compute_g_dpt`：每个 EE/PT 约束 +1，PP/PE 加重数（`ContactMechanics.cpp`），纯诊断用。

当前方向收敛时仍完成了一次 numeric factorization，因此 `numeric_factorizations` 通常比实际接受的 `newton` 步数多 1；这是判断当前解是否收敛所需的确认求解，不是漏计的 Newton 更新。

## 3. 可行步长上限（`:1990-2024`，`alpha` 初值 1，slackness=0.8）

1. `limitStepByGround(mesh, gd, moveDir, 0.8, alpha)`（`FeasibleStep.cpp`）：地面射线-平面解析解；`alpha<=0` 重置为 1。
2. `use_barrier` 时：
   - `Self_largestFeasibleStepSize`（部分 CCD）遍历经过统一静态 AABB+dHat 判定的 `Self_CCD_ActiveSet` → `partialCCD_alpha`；该集合与宽阶段后端无关；
   - CFL 风格上限 `alpha_CFL = sqrt(Hhat) / (2 · max_{surfVerts}|moveDir|)`（`:2012`）；
   - 宽阶段按入口 α 沿轨迹重建扫掠结构；PT ACCD 先更新安全 α，随后两个后端都用该更准确的 PT α 对 EE 做统一扫掠 AABB 精确过滤；ACCD 始终用完整方向返回全局参数，最终得到后端无关的 `fullCCD_alpha`；
   - 合并：`alpha = min(alpha, alpha_CFL)`；若 `partialCCD_alpha > 2·alpha_CFL`，`alpha = clamp(min(partial, full), ≥ alpha_CFL)`。
3. `filterStepSize`（tet 防翻转，体积收缩为 0.2 倍的三次方程根）不在当前主流程中调用。

## 4. `lineSearch`（`IPCSolver.cpp`）

1. `computeEnergy` 得当前能量 `lastEnergyVal`（公式见 §6）。
2. 按项目要求 `armijoCoefficient=0`，接受判据为严格单调 `E_trial < E0`，不使用允许能量上升的 roundoff tolerance；方向导数仍用于 Taylor 诊断。
3. 穿透防护：当 `isIntersected`（地面 gap ≤ 0，或边-三角相交），`stepSize /= 2` 并重建宽阶段，直到无穿透；现有 64 次上限和 step underflow 异常，不再可能静默死循环。
4. 能量回溯最多 64 次，每次 α 减半、重建碰撞集并重估能量；已删除旧的 `α≤1e−3·初始α` 就退出并可能接受能量上升步的逻辑。失败或 underflow 会抛出明确异常。`stepForward`：`x = V0 − α·searchDir`。
5. 步长被缩减后再做一次穿透防护（`:1753-1768`）。

## 5. 障碍刚度 κ 的生命周期

- `suggestKappa`（`:1775`）：`κ = 1e13·meanMass / (4e-16·bboxDiagSize2·H_b(1e-16·bboxDiagSize2, Hhat))` —— 标定准则：零距离处障碍 Hessian 与惯性项平衡。
- `upperBoundKappa`（`:1781`）：封顶 100×建议值。
- `initKappa`（`:1792`）：算弹性+惯性梯度 `g_E` 与单位 κ 接触梯度 `g_c`（地面 + 自碰撞 PT），`κ = −⟨g_c,g_E⟩/‖g_c‖²`（若为正），再钳到 `[suggest, upperBound]`。
- `postLineSearch`（`:1846`）：κ==0 → `initKappa`。否则检查记录的“接近约束”（`closeConstraintID/closeMConstraintID`）：若有任何距离较记录值**变小**，`κ *= 2`（再封顶）。随后把当前所有值 `< dTol` 的约束重新记录进接近列表（`:1877-1907`）。

## 6. 增量势能 E(x)（`computeEnergy`，`IPCSolver.cpp`）

```
E = dt²·(E_SNH(四面体) + E_BaraffWitkin(布料) + E_bend)
  + Σ_v ½·m_v·‖x_v − x̃_v‖²                (惯性, x̃ = V_prev + v·dt + g·dt², g=-9.8)
  + Σ_v ½·drag_coeff·m_v·‖x_v − V_prev‖²    (阻尼, 默认 drag=0)
  + κ·( Σ_ground b(d) + Σ_PT b(d)·重数 + Σ_EE e(x)·b(d) )
  + μ·Σ_c λ_c·f(‖u_t‖)                      (USE_FRICTION, f = 光滑钳制, ε² = Fhat·dt²)
```

`computeGradientAndHessian`（`:566-831`）按同一公式装配梯度/Hessian，细节：

- 惯性梯度 `m(x−x̃)`；质量对角 `m(1+drag)` 由 `NewtonLinearSystem::assemble` 加入 Hessian。
- tet：TBB parallel_for，`P = computePEPF_StableNHK3D_2_double(F, μ, λ)`，梯度 `f = volum·PFPXᵀ·vec(P)`（每顶点 spin_mutex 累加），Hessian `volum·dt²·PFPXᵀ·Hq·PFPX`（Hq 已 PD 投影）→ `BH.H12x12`。
- 布料：Baraff-Witkin，9×9 块 → `BH.H9x9`。
- 弯曲：两条路径共用 `plateRigidity=Et³/[12(1−ν²)]`。quadratic 使用预计算 `Q⊗I₃`；非 quadratic 使用预计算 `θ0`、`l0/(h0+h1)` 和二面角精确导数，Hessian 在装配前做 PSD 投影 → `H12x12`。
- 地面障碍：每激活顶点 `d=(n·x−D)²`，`g += κ·g_b(d)·2√d·n`；Hessian 仅当 `4H_b·d + 2g_b > 0` 才写入 `κ·param·nnᵀ`（逐项 PD 判定，不走 makePD）→ `H3x3`。
- 自碰撞 PT：`compute_g_b` 就地作用约束值后 `compute_g_dpt(...,Kappa)`；EE（mollified）：`compute_g_dee` / `compute_H_dee`；PT Hessian：`compute_H_dpt`。
- 摩擦（`USE_FRICTION`）：`Friction::addGradient/addHessian(mesh, grd, ..., eps2 = Fhat·dt², coef = mesh.friction)`；实现集中在 `Friction.cpp`。
- `computeEGradient`（`:528`）是仅弹性+惯性的变体，只被 `initKappa` 用。

## 7. 线性求解（`NewtonLinearSystem.cpp/.h`）

- 四个后端共用完全相同的 `BHessian` 下三角 triplet、质量对角和固定顶点 RHS 清零逻辑，避免后端间装配语义分叉。
- 回退 `LinearSolverBackend::SuiteSparseLDL`：当当前配置没有 PARDISO 时自动使用；将标量稀疏图折叠到每顶点 3-DOF block 做 AMD，展开 permutation 后只构建一次 `PAPᵀ` 上三角 CSC；同结构 Newton 只按缓存映射刷新 14.5 万个数值，再复用 symbolic 数据做 LDLᵀ numeric factorization。
- 可选 `LinearSolverBackend::Cholmod`：只有 CSC 结构变化才重新 analyze，保留 cost-aware AMD→METIS；使用 `cholmod_solve2` 复用 dense solution/workspace。
- 首选 `LinearSolverBackend::Pardiso`：当前配置定义 `CIPC_HAS_PARDISO` 时自动成为默认；使用 `mtype=2` 的实对称正定模式。Eigen 的 lower CSC 对称地解释为 upper CSR，避免转置/复制矩阵；phase 11 分析、22 数值分解、33 求解分别计时。同模式跨 Newton/κ/时间步跳过 phase 11；模式变化时优先复用已有 METIS permutation，若 factor nnz 超过最近 fresh ordering 的 1.2 倍则下一轮重新排序。默认限制为 16 线程；`--pardiso-threads 0` 使用 oneMKL 默认，其他正数通过 `tbb::global_control` 限制每个 phase。
- 可选 `LinearSolverBackend::EigenConjugateGradient`：Eigen CG + `IncompleteCholesky<Lower>`，容差/最大迭代来自 `LinearSolverOptions`。CLI 用 `--linear-solver eigen-cg`；可通过 twisting-mat 单步手工 smoke。它是正式可选项，不是死代码。
- MSVC + vcpkg oneMKL 的静态库使用 Release CRT，因此 PARDISO 只在非 Debug 配置定义 `CIPC_HAS_PARDISO`；普通 Debug 构建仍可使用其余后端，显式选 PARDISO 会抛出清晰错误。
- 旧的自定义 PCG、多 RHS `cholmod_solver_EPF` 与其 `mesh3D::Constraints` 投影状态已删除。

## 8. 每步结束

`updateVelocity`（v 更新）→ `V_prev = x` → `updateInertialTarget`。碰撞集不在时间步末重建；活动集在每次 Newton 迭代的线搜索中按需重建。
