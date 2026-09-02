# 03 — IPC 核心：距离、障碍函数、活动集与装配

涉及文件：`CPU IPC/ContactMechanics.cpp`（解析距离/障碍/导数）、`CollisionBroadPhase.cpp`（活动集）和 `IPCSolver.cpp`（总装配）。

## 1. 距离函数（全部是平方距离）

| 函数 | 位置 | 公式 |
|---|---|---|
| `d_PP` | `ContactMechanics.cpp` | `‖v0−v1‖²` |
| `d_PE` | `:863` | `\|(v1−v0)×(v2−v0)\|² / \|v2−v1\|²` |
| `d_PT` | `:353` | `((v0−v1)·b)² / \|b\|²`，`b = (v2−v1)×(v3−v1)` |
| `d_EE` | `:820` | `((v2−v0)·b)² / \|b\|²`，`b = (v1−v0)×(v3−v2)` |
| `computeEECrossSqNorm` | `:364` | `\|(v1−v0)×(v3−v2)\|²`（平行度，mollifier 用） |
| 地面 | `CollisionBroadPhase.cpp` | `(n·x − D)²` |

Additive CCD 内部保留自包含距离实现；未调用的备用 broad-phase 函数已删除。

## 2. 障碍函数（经典 IPC；d = 平方距离，dHat = `mesh.Hhat` = d̂²）

`ContactMechanics.cpp`：

```
b(d)  = −(d−dHat)² · ln(d/dHat)
b'(d) = −2(d−dHat)·ln(d/dHat) − (d−dHat)²/d
b''(d)= −2·ln(d/dHat) − 4(d−dHat)/d + (d−dHat)²/d²
```

仅在 `d < dHat` 时有定义——活动集构建时已按 `d < Hhat` 预过滤，调用方不做二次检查。

## 3. 最近特征分类

- `dType_PT`（`:199`）：0/1/2 = PP（对 t0/t1/t2），3/4/5 = PE（边 t0t1/t1t2/t2t0），6 = 真点-面。
- `dType_EE`（`:251`，Ericson 算法 + 近平行守卫 `:291-305`：`\|u×v\|² < 1e-20·a·c` 时回退端点情形）：
  0=PP(eI.0,eJ.0)，1=PP(eI.0,eJ.1)，2=PE(eI.0→eJ)，3=PP(eI.1,eJ.0)，4=PP(eI.1,eJ.1)，5=PE(eI.1→eJ)，6=PE(eJ.0→eI)，7=PE(eJ.1→eI)，8=真 EE。

## 4. 活动集构建（`SpatialHash::calculateActivateSet`；名称保留但宽阶段可为 LBVH）

默认查询由 CPU face/edge LBVH 提供，SpatialHash 可切回；两者都用半径 `sqrt(Hhat)` 做静态候选，精确距离判定为 `d < Hhat`。并行于 `surfVerts`（PT）和 `surfEdges`（EE，`eJ > eI` 去重）；候选输出排序，线程内 scratch 复用。

过滤规则：
- 4 个顶点全部满足 `isExternalColliderBoundary(boundaryTypes[v])`（code ≥2，均为不参与求解的外部 prescribed collider）→ 跳过；普通 Dirichlet 顶点(code 1)及软边界自由顶点不会被这条规则误删；
- 点-面：点是三角形角点 → 跳过；
- 边-边：共享顶点且 `eI > eJ` → 跳过。

产出（编码规则见 `01_overview.md` 的 `EncodedContact` 节）：
- `Self_ActiveSet`：PT + 常规 EE + 去重后的 PP/PE（`[3]=-count`）；
- `Self_EE_ActiveSet` + `Self_EEeIe_ActiveSet`（(eI,eJ) 索引对，真 EE 时为 (-1,-1)）：近平行对，`eps_x = 1e-3·|restE_I|²·|restE_J|²`（`compute_eps_x`，`:814-818`）；
- `Self_CCD_ActiveSet`：宽阶段候选再经过统一的静态 AABB overlap（padding=`sqrt(Hhat)`）、共享顶点和全驱动过滤后的 PT/EE primitive pairs，编码为 `(-svI-1,sfI)` 或 `(eI,eJ)`。该二次过滤消除了 SpatialHash 体素假阳性与 LBVH 紧 AABB 候选差异，确保 partial CCD/CFL 分支输入一致。

摩擦切空间基**不在这里**算，而在每时间步的 `Friction::initialize`（`Friction.cpp`）按 `Self_ActiveSet` 逐约束计算并冻结。

地面侧：`Ground::calculateActivateSet`（`CollisionBroadPhase.cpp`）并行遍历 `surfVerts`，gap² < Hhat 的全局顶点 id 进 `Environment_ActiveSet`。

## 5. 梯度装配（活跃路径）

- `Evaluate_GroundConstraintVals / _SelfPTConstraintVals / _SelfEEConstraintVals`（`:2561-2609`）：并行求值三组约束（平方距离），按 `offset` 追加（协议：地面 → 自 PT → 自 EE 依次堆叠）。
- `compute_g_dpt`（返回碰撞计数）：逐 `Self_ActiveSet` 条目算距离梯度，乘 `κ·b'(d)` 与 PP/PE 重数 `-contact[3]`，散到每顶点梯度。
- `compute_g_dee`（`:3090`，近平行 EE）：`∇(κ·e·b) = κ·(b·∇e + e·b'·∇d)`；先加 `κ·b·∇e`，再逐约束调 `compute_g_dpt` 传 `e·b'(d)`。

## 6. Hessian 装配

- `compute_H_dpt`（`:3522`，TBB parallel_for + 3 互斥锁）：逐条目 `H_IP = (κ·b'')·g·gᵀ + (κ·b')·H_d`，然后 `IglUtils::makePD`（6×6/9×9/12×12）→ 推入 `BH.H6x6/H9x9/H12x12` + 索引表。
- `compute_H_dee`（`:3631`）：完整乘积法则
  `H = κ·b'·(∇d∇eᵀ + ∇e∇dᵀ) + κ·b·H_e + κ·e·b''·∇d∇dᵀ + κ·e·b'·H_d`，再 `makePD<12>`。EE 集合里的 PP/PE 条目通过顶点映射 `indMap`（`:3692-3741`）把 6×6/9×9 嵌入 12×12。
- 各块最终经 `BHessian::toTriplets(boundaryTypes)` 与弹性块、质量对角一起进稀疏矩阵（机制见 `01_overview.md`）。

`IglUtils::makePD`（`ContactMechanics.h`）：SelfAdjointEigenSolver；最小特征值 ≥0 直接返回；否则把负特征值钳为 0 重建 `V·D·Vᵀ`。接触障碍块与部分弹性路径仍使用它；解析摩擦 Hessian 和二次弯曲常量 Hessian不使用通用投影。

## 7. 解析梯度/Hessian 的出处

大段代码是 **MATLAB Symbolic Math Toolbox 8.3** 生成（保留出处注释）：

- `computeEECrossSqNormGradient/Hessian`（`:374/:449`，2019-11-01）
- `g_PE`（`:871`）、`g_PT`（`:926`）、`H_PE`（`:1002`）、`H_PT`（`:1335`，2019-06-10）
- `g_EE`（`:1806`）、`H_EE`（`:1900`，约 600 行到 `:2511`，2019-06-14）

`g_PP/H_PP` 是手写的平凡形式（`:848-861`）。每个标量数组内核后有 Eigen 包装函数。

## 8. EE mollifier e(x)（`compute_e`，`:751`）

```
e = 1                                        若 EECrossSqNorm ≥ eps_x
e = (2 − c/eps_x)·(c/eps_x)  (q(c), c=crossSqNorm)  否则
```

`compute_e_g / compute_e_H`（`:768/:788`）用链式法则过生成的 crossSqNorm 梯度/Hessian：`q' = 2(1−c/eps_x)/eps_x`，`q'' = −2/eps_x²`。

## 9. 未接入的替代障碍内核

- `compute_g_dpt_new` / `compute_H_dpt_new` 仍保留在 `ContactMechanics.cpp/.h`：另一套基于“偏移最近特征构造 F + 不变量 I5 = nᵀFᵀFn”的实验障碍表述。
- 当前 `IPCSolver.cpp` 没有调用点；若未来启用，必须同时补能量、梯度、Hessian 一致性测试。
