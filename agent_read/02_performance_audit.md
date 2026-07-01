# 性能热点审计

> 只读审计，未修改任何文件。  
> 重点：识别计算密集区，并标注**数值安全**的优化方向。

## 1. 现有性能剖析

`IPC_FUNC.cpp` 已在 `solve_subIP` 中使用 `HighResolutionTimerForWin` 分段计时：

| 阶段 | 计时变量 | 主要函数 |
|------|---------|---------|
| 梯度/海森组装 | `time0` | `computeGradientAndHessian` |
| 线性求解 | `time1` | `calculateMovingDirection` / `cholmod_solver` |
| CCD 步长 | `time2` | `Self_largestFeasibleStepSize` / `_CCD` |
| 线搜索 | `time3` | `lineSearch` |
| 后处理 / κ 更新 | `time4` | `postLineSearch` |

计时器粒度较粗，未细分 FEM / 障碍 / 摩擦内部耗时。

## 2. FEM 弹性核函数（`fem3D.cpp`）

| 函数 | 主要开销 |
|------|---------|
| `getObjEnergy_StableNHK2_3D` | 每 tet SVD（`QRSVD`） |
| `computePEPF_StableNHK3D_2_double` | 每 tet SVD + 手动 cofactor 矩阵 |
| `project_StabbleNHK_2_H_3D` | 每 tet 构建 9×9 解析海森 + `IglUtils::makePD<double,9>`（完整特征分解） |
| `computePFPX3D_double` | 返回 `MatrixXd`（9×12），每次调用堆分配；大小本可在编译期确定 |
| `vec_double` / `vec_float` | 向量化时动态分配 `MatrixXd` |

### 优化方向（数值安全）

1. **预计算 `PFPX`**：`DM_tetrahedra_inverse` 恒定，可在初始化时存储 `Matrix<double,9,12>` 与 `Matrix<double,6,9>`，避免每迭代重算。
2. **固定尺寸替代 `MatrixXd`**：9×12、6×9、12×1 等全部用 `Eigen::Matrix<double, M, N>`，避免堆分配并提升向量化。
3. **降低/移除每 tet 的 PSD 投影**：Stable Neo-Hookean 解析海森在可逆构型下本身 PSD，可删除 `makePD` 或改用更便宜的 3×3 拉伸块截断。
4. **`vec_double` 改为固定尺寸 Map**：使用 `Map<Matrix<double,9,1>>` 等。

## 3. IPC 障碍与距离函数（`IPCdistanceFuncs.cpp`）

| 函数 | 主要开销 |
|------|---------|
| `compute_g_dpt` | **串行**遍历 `Self_ActiveSet` 累加障碍梯度 |
| `compute_g_dee` | **串行**遍历 `Self_EE_ActiveSet` |
| `compute_H_dpt` | 并行遍历活动接触，但每个接触调用 `makePD`，并在 `spin_mutex` 下 push_back |
| `compute_H_dee` | EE 接触 + mollifier，符号表达式 + `makePD` |
| `Evaluate_*ConstraintVals` | `conservativeResize` + TBB 并行距离计算 |
| `Self_largestFeasibleStepSize_CCD` | 全对 CCD，每个 primitive 分配 `unordered_set` |

### 优化方向（数值安全）

1. **并行化障碍梯度累加**：为每个线程分配局部 `vector<Vector3d>`，最后确定性归约。只要保持与串行相同的求和顺序，结果不变。
2. **预分配 `BHessian` 容量**：在 `compute_H_dpt`/`compute_H_dee` 前 `reserve` 各尺寸块，避免锁下 reallocation。
3. **缓存 `eps_x`**：相同边对在 EE mollifier 中重复计算，可缓存。
4. **用 sort/unique 替代 `unordered_set`**：`calculateActivateSet` 与 CCD 中每查询新建 `unordered_set`，大量小内存分配。
5. **`std::map<MMCVID,int>` 去重替换为 hash map 或 sort-merge**。

## 4. 空间哈希与 CCD（`collisionUtil.cpp`、`ACCD.cpp`）

| 函数 | 主要开销 |
|------|---------|
| `SpatialHash::build` | 构建 `unordered_map<int, vector<int>>` 体素网格 |
| `SpatialHash::calculateActivateSet` | 每个表面点/边查询邻近 primitive，分配 `unordered_set`/`map` |
| `queryPointForTriangles` / `queryEdgeForEdgesWithBBoxCheck` | 每查询分配 `unordered_set` |
| `point_triangle_ccd` / `edge_edge_ccd` | while 循环中每次重新进行距离类型分类 |
| `Self_largestFeasibleStepSize_CCD` | 全对 CCD，每 primitive 分配集合 |

### 优化方向（数值安全）

1. **线程局部预分配 flat buffer + sort/unique** 替代 `unordered_set`。
2. **CCD 距离类型一次性分类**：在 while 外先分类，循环内复用。
3. **避免线搜索内重复 build 空间哈希**。

## 5. 线性求解器（`Solver.cpp`、`IPC_FUNC.cpp`）

| 函数 | 主要开销 |
|------|---------|
| `cholmod_solver` | 每次牛顿迭代**完整符号+数值分解**（`set_pattern` + `solve`） |
| `PCG_Solver` | 手卷块 PCG，矩阵向量乘使用每顶点 `spin_mutex` |
| `PCG_Precondition` | 手卷 3×3 高斯约旦求逆（`__Inverse2`） |

### 优化方向（数值安全）

1. **复用 CHOLMOD 符号分解**：活动接触模式不变时，先 `preFactorize` 一次，后续用 `solve_with_preFactorize` 只更新数值。
2. **模式变化时仅重新 analyze**：`set_pattern` 后应只在 pattern 变化时调用 `cholmod_analyze`。
3. **PCG 锁优化**：使用 atomic 累加或网格着色消除 `spin_mutex` 竞争。
4. **3×3 逆用 Eigen 固定尺寸 `inverse()`** 替代手卷高斯约旦。

## 6. 内存分配热点

| 来源 | 分配模式 | 影响 |
|------|---------|------|
| `computePFPX3D_double` / `computePFPX32D_double` | 每元素返回 `MatrixXd` | 每梯度/能量/海森调用都堆分配 |
| `vec_double` / `vec_float` | 返回 `MatrixXd` | 每次向量化都堆分配 |
| `BHessian BH;` 每次牛顿迭代重建 | 所有元素/接触海森存储重新分配 | 大量 realloc |
| 每个 collision query 新建 `unordered_set<int>` | 大量小分配 | 分配器压力 |
| `std::map<MMCVID,int>` | 树插入 + 分配 | 慢且分配频繁 |
| `Eigen::VectorXd` / `MatrixXd` 临时变量 | PCG 与海森块 mat-vec 中 | 栈/堆临时对象 |
| `cholmod_solver` / `Eigen_CG_solver` | 每次重建 triplet buffer 与稀疏矩阵 | 额外分配 |

## 7. SVD / 稠密线性代数（`SVD.cpp`）

- `QRSVD` 每 tet 在梯度、海森、能量评估中多次调用。
- 可使用非 SVD 的不变量形式 Stable Neo-Hookean 来消除 SVD 调用，同时保持能量/梯度/海森值相同。

## 8. 优先级最高的优化建议（保持数值一致）

1. **消除 FEM 核中的每元素堆分配**（`PFPX` 预计算、固定尺寸矩阵、`vec_double` 优化）。
2. **移除或廉价化每 tet 的 `makePD`**。
3. **并行化串行的障碍梯度累加**。
4. **预分配 `BHessian` 容量**。
5. **复用 CHOLMOD 符号分解**。
6. **消除 collision 查询中的 `unordered_set` 小分配**。
7. **CCD 中一次性分类距离类型**。
8. **改善 TBB 划分粒度与锁策略**。
9. **编译器向量化**：在 CMake 中加入 `/arch:AVX2`（MSVC）或 `-march=native`（GCC/Clang），并保持 `/fp:precise`。
