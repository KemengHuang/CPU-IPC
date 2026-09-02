# 10 — oneMKL PARDISO 与 bunny2 专项报告

本报告记录 oneMKL PARDISO 后端的设计、可复现基准、数值核对和采用边界。测试日期为 2026-09-02；绝对时间仅代表本开发机，跨机器应比较趋势并重新测量。

## 1. 先明确“一个 step”

本文的 `--steps 1` 指一个完整 IPC **simulation time step / frame**，不是一次 Newton step。一个时间步内部包含 κ 外层、若干 Newton 迭代、Hessian 装配、线性分析/分解/求解、CCD、CFL 与严格下降线搜索。

bunny2 首时间步为 3 次 Newton，因此 PARDISO phase 22 与 phase 33 各执行 3 次；phase 11 因模式复用只执行 1 次。`step_ms` 是整帧墙钟时间，`linear_ms` 包括 triplet/CSC 构造与全部 solver 工作，三个 `pardiso_*_ms` 只统计对应 PARDISO phase 的累计墙钟时间。

## 2. 构建与接口

- `CIPC_ENABLE_PARDISO=ON` 是推荐保持开启的默认主路径。找到 `MKL CONFIG` 时定义 `CIPC_HAS_PARDISO`、编译 `PardisoSolver.cpp` 并默认选择 PARDISO；找不到时仅为兼容性继续构建 SuiteSparse LDL、CHOLMOD 与 Eigen-CG。
- Windows/vcpkg：`vcpkg install intel-mkl:x64-windows`。当前接入固定 `MKL_LINK=static`、`MKL_INTERFACE=lp64`、`MKL_THREADING=tbb_thread`。
- 当前 vcpkg oneMKL 2025.2 静态包使用 Release CRT。为避免 `/MD` 与 `/MDd` 的 `LNK2038`，MSVC Debug 配置自动只排除 PARDISO；Release、RelWithDebInfo 和其余 Debug 产品仍正常。静态链接后的 Release `cipc_headless.exe` 在本机为 74,003,968 bytes（约 74.0 MB / 70.6 MiB）。
- CLI：`--linear-solver pardiso --pardiso-threads N`。运行时默认 `N=16`，显式 `N=0` 才采用 oneMKL 默认；正数通过每个 phase 外层的 `tbb::global_control(max_allowed_parallelism)` 生效。直接调用 `mkl_set_num_threads` 对 TBB threading layer 的实测限制无效，因此未采用。
- `scripts/benchmark.py` 已支持 `--scene bunny2`、`--linear-solver auto|pardiso|...` 和 `--pardiso-threads`；solver 默认 `auto` 以继承 executable 的能力选择，PARDISO 线程默认同样为 16。

官方语义参考：oneMKL 的 [PARDISO 接口与并行直接法说明](https://www.intel.com/content/www/us/en/docs/onemkl/developer-reference-c/2026-0/onemkl-pardiso-parallel-direct-sparse-solver-iface.html)、[phase 参数](https://www.intel.com/content/www/us/en/docs/onemkl/developer-reference-c/2026-0/pardiso.html) 与 [iparm/permutation 参数](https://www.intel.com/content/www/us/en/docs/onemkl/developer-reference-c/2026-0/pardiso-iparm-parameter.html)。

## 3. 求解器实现

### 矩阵与 phase

- Newton Hessian 已由项目统一 PSD 处理并加正质量对角，PARDISO 使用 real symmetric positive-definite `mtype=2`。
- 项目共享矩阵是 Eigen column-major lower CSC。对称矩阵的 lower CSC 在内存上正好等价于同一矩阵的 upper CSR，因此可把 `value/outer/inner` 直接交给 PARDISO，避免 O(nnz) 转置和第二份矩阵。
- LP64、double、zero-based C indexing：`iparm[27]=0`、`iparm[34]=1`，CSR storage；分析/数值分解/求解分别调用 phase 11/22/33，析构或结构切换调用 phase −1。
- RHS 复制到 solver 自有 buffer，结果统一 scatter 回每顶点方向；非有限结果与非零 PARDISO error 均立即抛异常。

### 生命周期与 permutation

`NewtonLinearSystem` 现由 `IPCSolverContext` 持有，不再局限于一个 barrier 子问题，因此 sparse buffer 与三个直接法的 symbolic 状态可跨 Newton、κ 外层和时间步复用。每帧 metrics 用调用前后的累计计数/计时差值，避免跨帧重复计账。

最终两步 smoke 中 frame 0 为 `analyze=1`，相同 pattern 的 frame 1 为 `analyze=0`；第二帧 phase 22/33 约为 89/61 ms，验证了复用不是只停留在代码结构上。

PARDISO 完整比较 outer/inner pattern：

1. 初次分析让 PARDISO 用 METIS nested dissection，并返回 permutation。
2. 维度相同但接触使 pattern 变化时，先复用该 permutation，减少重新排序成本。
3. 若 phase 22 报告的 factor nnz 超过最近一次 fresh ordering 参考值的 1.2 倍，则标记下一次求解做 fresh METIS；新 fill 成为下一参考。

这个策略不是只追求最少 phase 11，而是在排序时间和 fill/factorization 之间取实测平衡。

### 指标

`metrics.csv` 新增：

```text
pardiso_analysis_ms,pardiso_factorization_ms,pardiso_solve_ms,
linear_solver_threads,factor_nnz
```

`symbolic_analyses` 与 `numeric_factorizations` 仍是跨直接后端的统一口径。非 PARDISO 后端的新字段为 0。

## 4. bunny2 场景

场景参考 `GPU_IPC/GPU_IPC/gl_main.cpp` 的 `initScene1`：

- 两次载入 `Assets/tetrahedraMesh/bunny2.msh`；
- 两只均 `scale=0.2`；
- offset 分别为 `(0,0.65,0)` 与 `(0,0,0)`；
- `YoungModulus=1e5`。

单份网格为 19,193 顶点 / 79,935 tet；合并场景为：

| 指标 | 数值 |
|---|---:|
| vertices | 38,386 |
| tetrahedra | 159,870 |
| surface triangles | 41,664 |
| Newton DOF | 115,158 |
| 首步 matrix nnz | 2,202,090 |
| 首次 PARDISO factor nnz | 25,102,596 |

## 5. 测试环境与方法

- CPU：AMD Ryzen 9 9950X3D，16 cores / 32 threads。
- 系统/编译：Windows、MSVC、Release。
- oneMKL：2025.2.0，static LP64 + TBB threading。
- broad phase：默认 LBVH；物理参数、严格能量下降、原 CFL、顺序 PT→EE CCD、Hessian scale=1 均保持不变。
- 时间表除特别注明外为 3 个独立进程的中位数；多后端采用交替运行以降低温度/调度偏差。

复现示例：

```bash
python scripts/benchmark.py \
  --exe build/Release/cipc_headless.exe \
  --scene bunny2 \
  --linear-solver pardiso \
  --pardiso-threads 16 \
  --steps 1 \
  --repeats 3
```

## 6. PARDISO 线程扩展

bunny2、一个完整时间步、每档 3 次中位数：

| threads | step_ms | linear_ms | phase 11 ms | phase 22 ms | phase 33 ms |
|---:|---:|---:|---:|---:|---:|
| 1 | 2880.008 | 2673.688 | 503.040 | 1428.984 | 200.212 |
| 2 | 2008.860 | 1793.839 | 472.109 | 701.393 | 124.283 |
| 4 | 1612.079 | 1418.935 | 454.023 | 361.720 | 108.259 |
| 8 | 1488.740 | 1287.906 | 469.836 | 215.989 | 101.901 |
| 16 | 1409.121 | 1215.252 | 460.533 | 150.185 | 102.699 |
| 32 | 1404.253 | 1229.308 | 471.068 | 151.777 | 105.164 |

结论：phase 22 从 1→16 线程加速 9.51×，但完整 time step 只加速约 2.04×，因为 phase 11、Eigen CSC 构造、装配/CCD/线搜索不随 PARDISO factorization 同比例缩放。32 线程的总数值落在噪声内，`linear_ms`、phase 22/33 反而略高；5 步测试的 32 线程中位数 5.104 s 也慢于 16 线程约 4.786 s，因此本机推荐 16。

## 7. 与现有直接法比较

### 大 bunny2：1 个时间步，3 次中位数

| solver | step_ms | linear_ms | 相对 PARDISO |
|---|---:|---:|---:|
| PARDISO, 16 threads | 1284.371 | 1167.917 | 1.00× |
| CHOLMOD（当前非 supernodal vcpkg） | 8253.445 | 8143.121 | 6.43× slower |
| SuiteSparse LDL | 12568.366 | 12456.280 | 9.79× slower |

三者均为 3 Newton；位置和、平方范数和、nnz 与接触指标在浮点舍入内一致。这里的巨大收益来自大 SPD sparse factorization，并不表示 PARDISO 在所有小场景都必然有相同比例。

### 中型场景：5 个完整时间步，3 次中位数

| scene / solver | total step_ms | total linear_ms | PARDISO 加速 |
|---|---:|---:|---:|
| cloth-bunny / PARDISO(16) | 585.018 | 450.187 | — |
| cloth-bunny / SuiteSparse LDL | 973.519 | 837.246 | 1.66× |
| cloth-bunny / CHOLMOD | 1049.260 | 910.539 | 1.79× |
| twisting-mat / PARDISO(16) | 467.444 | 338.971 | — |
| twisting-mat / SuiteSparse LDL | 628.311 | 506.454 | 1.34× |
| twisting-mat / CHOLMOD | 674.331 | 544.592 | 1.44× |

这组数据已包含 `NewtonLinearSystem` 跨 time step 复用；该生命周期优化也同时帮助 LDL/CHOLMOD，比较不是只给 PARDISO 特权。

## 8. 内存

bunny2 单时间步、50 ms 轮询的一次运行记录：

| solver | wall | process CPU | peak working set | peak private |
|---|---:|---:|---:|---:|
| PARDISO(16) | 1.70 s | 8.97 s | 1050.9 MiB | 1229.7 MiB |
| CHOLMOD | 8.62 s | 13.39 s | 1045.2 MiB | 1146.8 MiB |
| SuiteSparse LDL | 13.0 s | 17.86 s | 1085.8 MiB | 1180.2 MiB |

PARDISO 的 working set 与两者相近，private bytes 比 CHOLMOD 高约 7%；这是单次轮询样本，不应解释成精确分配剖析。

## 9. 50 步接触与 ordering 策略（旧停止点基线）

以下三组 ordering A/B 在“用上一轮方向判断收敛”的旧停止点下完成，用于比较 permutation 策略；当前方向判定改变了绝对 Newton/回退数，不能把表中总时间当作现版本基线。bunny2 前 42 帧主要检验跨帧 reuse；frame 43 首次出现 ground=694/self=3，frame 49 为 ground=441/self=1006。三个排序策略都得到相同整数接触指标和正的最小距离：

| 策略 | total | linear | phase 11 | phase 22 | phase 33 | 末帧 factor nnz |
|---|---:|---:|---:|---:|---:|---:|
| pattern 变化总是 fresh METIS | 73.033 s | 58.358 s | 18.400 s | 8.424 s | 5.595 s | 27.45 M |
| 永久复用首份 permutation | 62.838 s | 49.884 s | 8.393 s | 11.188 s | 6.519 s | 46.69 M |
| adaptive 1.2× fill（当前） | 59.167 s | 47.964 s | 8.791 s | 8.884 s | 6.341 s | 35.60 M |

adaptive 比 always-fresh 快 19.0%，比 permanently-fixed 快 5.8%。它允许少量额外 phase 11 来阻止固定 permutation 的 fill 膨胀。

旧停止点 adaptive 50 步末态：

```text
sum_x             = -4725.2909366452941
sum_y             = -24659.348104246626
sum_z             = 1568.4896282611007
squared_norm_sum  = 20347.976494341674
min_distance2     = 7.3483170186983491e-7
matrix_nnz        = 2214465
```

always-fresh 与 adaptive 的坐标和差约 `1e-10`，碰撞整数指标一致，最小距离差约 `1e-18`，属于并行归约/直接法浮点路径差异。

改为用本轮刚求出的方向判断后，在保持原 `beta` 辅助退出的隔离50步运行中，总时间54.716 s、153个接受Newton步、15次能量回退；末帧ground/self=441/999、`min_distance2=7.3487834099833287e−7>0`、`matrix_nnz=2,214,492`。与旧停止点相比轨迹和接触有小幅变化，但仍保持严格下降与正距离。

## 10. 正确性与当前方向收敛

- PARDISO、CHOLMOD、SuiteSparse LDL 使用同一 Hessian、RHS、边界处理、CCD pair、CFL 与 line search；后端只替换线性求解。
- 首步跨后端结果在浮点舍入内一致；50 步 ordering A/B 保持接触数和正最小距离。
- `armijoCoefficient=0`、严格 `E_trial<E0`、摩擦 Hessian scale=1 均未改变。
- 无接触且接近机器精度时，旧停止顺序会把本轮近零方向继续送入严格比较：实测 `gTp≈7.33e−30`、能量只在约 `1e−16` 末位摆动，单帧可随机二分至52次。按项目当前要求，正式顺序改为先完成本轮线性求解，再用该当前方向判断；收敛时不进入CCD/line search。5次独立运行×前5帧的回退由417降为0，终态完全一致，5步中位由4429.220降至3396.438 ms。该调整不允许相等/上升能量，也没有修改Hessian或scale。
- PARDISO `iparm[1]=3` 的 parallel nested dissection 在当前 MKL/TBB/场景上更慢，已恢复 METIS `iparm[1]=2`。

## 11. 采用建议与后续瓶颈

- 大型 CPU benchmark 或 10 万 DOF 级 SPD Newton 系统：优先试 `pardiso --pardiso-threads <physical cores>`；本机为 16。
- 当前配置提供 PARDISO 时，能力感知默认值直接选择 PARDISO（16线程）；无 oneMKL、显式关闭或 MSVC Debug 时自动回退 SuiteSparse LDL。仍可用 `--linear-solver suitesparse-ldl` 强制选择较小、可移植的后端。
- 下一性能目标应是直接构建/更新 CSC、减少 `setFromTriplets` 与 pattern comparison 成本，以及降低 phase 11 前后的串行工作；在本机继续从 16 增到 32 线程没有收益。
- 若评估固定 superset pattern，必须同时报告 symbolic 次数、factor nnz、峰值内存、Newton/线搜索和最终轨迹；少做分析但让 fill 长期膨胀并不是净优化。

后续GPL supernodal+oneMKL CHOLMOD已把与PARDISO的差距显著缩小；本文早期CHOLMOD绝对时间属于旧OpenBLAS路径。当前同版本对比见`11_cholmod_mkl_report.md`。
