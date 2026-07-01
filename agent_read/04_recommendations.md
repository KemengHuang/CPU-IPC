# 整理 / 重构 / 优化 优先级建议

## 第一阶段：代码梳理与可读性（任务 2）

目标：在不改变数值行为的前提下提升可读性。

### 高优先级

1. **删除死代码**
   - `FEMIntegrator::solveCollision` 未定义声明
   - `mesh.h` 中 `model_obj`、`fiber_obj`、`mesh2D` 等未使用类型
   - `gl_main.cpp` 中 `fiberObj`、`drawFiber`、`drawSurface` 等无效开关
   - `Simulator.cpp` 中 `case1` 或实现场景分发
   - `fem3D.cpp` 中 `__Inverse` / `__Inverse2x2` 死代码
   - `IPC_FUNC.cpp` 中 `getEKF` / `test_jacobian` 大段注释代码

2. **修复拼写与命名**
   - `averageEdgeLenth` → `averageEdgeLength`
   - `minConer` / `maxConer` → `minCorner` / `maxCorner`
   - `spectialPontsArray` → `specialPointsArray`
   - `xTilta` → `xTilde`
   - `Aniostropic` → `Anisotropic`
   - `StabbleNHK` → `StableNHK`
   - `IPCtimeStepFuns.h` → `IPCtimeStepFuncs.h`（同步 `.cpp` 包含）

3. **消除头文件污染**
   - 移除头文件中的 `using namespace std;` / `using namespace Eigen;`
   - 统一使用 `#pragma once`
   - 减少不必要包含，使用前向声明

4. **const 正确性与传参优化**
   - `Matrix2d`/`Matrix3d`/`Vector3d` 等大型 Eigen 对象改为 `const&`
   - 审计只读函数加 `const`

5. **移除 `exit(-1)` / `system("pause")` / `exit(0)`**
   - 网格 I/O 返回 bool/异常
   - 求解器返回状态码

### 中优先级

6. **拆分过长函数**
   - `solve_subIP` → 初始化、梯度海森、求解方向、CCD、线搜索、后处理
   - `computeGradientAndHessian` → FEM 弹性、壳、弯曲、接触、摩擦辅助函数
   - `computeEnergyVal` → 各能量项辅助函数
   - `SpatialHash::calculateActivateSet` → 每 case 辅助函数

7. **集中魔法数字**
   - `SimulationParameters` 结构体
   - `kDefaultTimeStep`、`kPCGTolerance`、`kMaxCGIterations`、`kCFLFactor` 等命名常量

## 第二阶段：重构（任务 3）

目标：降低耦合、提高可维护性。

### 架构层

1. **分离可视化、场景配置与求解库**
   - `gl_main.cpp` 仅负责渲染与事件循环
   - 将 `case1`/`case2` 场景配置抽离到 `Scenes/` 或 `SceneBuilder`
   - 求解器核心尽量不依赖 OpenGL

2. **抽象线性求解器接口**
   - 定义 `LinearSolver` 接口
   - `CholmodSolver` 实现该接口
   - 将 `cholmod_solver`、`Eigen_CG_solver`、`PCG_Solver` 统一封装

3. **统一稀疏/块海森表示**
   - 简化 `BHessian::toTriplets`
   - 用模板泛化 D1/D2/D3/D4 处理

4. **减少 `mesh3D` 职责**
   - I/O、FEM 数据、IPC 状态、摩擦状态可逐步分离
   - 但至少先集中参数到 `SimulationParameters`

### 数值安全层

5. **保留符号生成代码为只读模块**
   - 将 `g_*`、`H_*` 等函数集中并标注"auto-generated, do not edit"
   - 如需修改，通过脚本重新生成

6. **引入回归测试骨架**
   - 每 10 帧检查点对比
   - 单步/多步 `saveSurface/` 顶点位置对比
   - `timeCost.txt` 关键指标对比

## 第三阶段：性能优化（任务 4）

目标：在数值一致的前提下提升运行效率。

### 高影响力、低风险

1. **预计算 `PFPX`**
   - 在 `mesh3D` 初始化时计算并存储 `Matrix<double,9,12>`（tet）和 `Matrix<double,6,9>`（tri）
   - 从 `computeGradientAndHessian`、`computeEGradient`、`getObjEnergy_StableNHK2_3D` 中移除 `computePFPX*` 调用

2. **固定尺寸 Eigen 替代 `MatrixXd`**
   - `PFPX`、`vec_double` 返回类型、临时 `VectorXd(12)` 等
   - 消除每元素堆分配

3. **预分配 `BHessian` 容量**
   - 在 `solve_subIP` 中根据 `FEM_blocks + active_set_size` reserve 各向量

4. **并行化串行障碍梯度累加**
   - `compute_g_dpt`、`compute_g_dee` 改为线程局部缓冲区 + 确定性归约
   - 保持求和顺序不变

5. **复用 CHOLMOD 符号分解**
   - 检测活动接触模式是否变化
   - 模式不变：一次 `preFactorize` + 多次 `solve_with_preFactorize`
   - 模式变化：重新 analyze

### 中影响力

6. **降低/移除每 tet `makePD`**
   - Stable Neo-Hookean 解析海森在可逆区已 PSD
   - 可改为可选的早期退出或 3×3 拉伸块截断

7. **消除 collision 查询中的小分配**
   - 用线程局部 `std::vector<int>` + sort/unique 替代 `unordered_set`
   - 用 hash map 或 sort-merge 替代 `std::map<MMCVID,int>`

8. **CCD 内循环优化**
   - 在 while 外一次性分类距离类型
   - 循环内复用

9. **TBB 使用优化**
   - 使用 `blocked_range` + 合理 grain size
   - 梯度/PCG 中减少 `spin_mutex` 竞争（atomic 或着色）

### 编译器/构建层

10. **启用向量化指令**
    - MSVC: `/arch:AVX2`
    - GCC/Clang: `-march=native`
    - 保持 `/fp:precise` 而非 `/fp:fast`，确保数值一致

## 数值验证策略

每次改动后执行：

1. 编译通过（Debug 与 Release）。
2. 运行默认场景（`case2`：布料落向 bunny）至少 50 帧。
3. 对比以下指标与 baseline：
   - `saveSurface/surf_XXXXX.obj` 顶点位置（L2 误差 < 1e-6 相对或绝对）
   - `timeCost.txt` 中每步牛顿迭代数、碰撞数、总能量
   - 无负体积
   - 线搜索后无自相交
   - κ 演进一致

若出现差异，立即回退并用更小粒度（函数级）逐步验证。
