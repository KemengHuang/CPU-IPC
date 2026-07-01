# 代码质量与可读性审计

> 只读审计，未修改任何文件。

## 1. 死代码 / 未使用成员

| 问题 | 文件 | 位置 | 建议 |
|------|------|------|------|
| `FEMIntegrator::solveCollision()` 声明但未定义 | `FEMTimeIntegrator.h` | ~17 | 删除或实现 |
| `model_obj`、`fiber_obj`、`mesh2D` 等类型未使用 | `mesh.h` | ~10–271 | 删除或归档 |
| `fiber_obj fiberObj` 全局变量未使用 | `gl_main.cpp` | ~34 | 删除 |
| `FEMSimulator::triangle_meshes` 从未填充 | `Simulator.h` | ~27 | 删除或接入 |
| `case1` 从未调用；`buildType`/`sceneType` 参数被忽略 | `Simulator.cpp` | ~195–235, ~283–291 | 实现分发或移除参数 |
| `drawFiber`、`drawSurface` 等布尔开关无效 | `gl_main.cpp` | ~36–42, ~437–443 | 删除或使用 |
| `__Inverse` / `__Inverse2x2` 未使用（使用 Eigen `.inverse()`） | `fem3D.cpp` | ~427–537 | 删除 |
| `computeGradientAndHessian` 返回 0 无意义 | `IPC_FUNC.cpp` | ~566, ~1971 | 改为 `void` 或真正计数 |
| 大段注释掉的 `getEKF` / `test_jacobian` | `IPC_FUNC.cpp` | ~1001–1063 | 删除 |

## 2. 命名不一致 / 拼写错误

| 错误/不一致 | 建议 | 位置 |
|------------|------|------|
| `averageEdgeLenth` | `averageEdgeLength` | `mesh.h` 等 |
| `minConer` / `maxConer` | `minCorner` / `maxCorner` | `mesh.h` 等 |
| `spectialPontsArray` | `specialPointsArray` | `mesh.h` 等 |
| `xTilta` | `xTilde` | 多处 |
| `Aniostropic` | `Anisotropic` | `fem3D.h/cpp` |
| `StabbleNHK` | `StableNHK` | `fem3D.h/cpp` |
| `baraffwitkin` / `baraffwitkint` | `baraffWitkin` | `fem3D.h/cpp` |
| `tetrahedras` | `tetrahedra` | `mesh.h` 等 |
| `IPCtimeStepFuns.h` | `IPCtimeStepFuncs.h` | 文件名 |
| 混合 `snake_case`、`camelCase`、`PascalCase` | 统一约定 | 全项目 |

## 3. 过长函数 / 重复代码

| 函数 | 文件 | 约行数 | 建议 |
|------|------|--------|------|
| `solve_subIP` | `IPC_FUNC.cpp` | ~312 | 拆分为梯度/海森组装、方向求解、CCD/线搜索、后处理 |
| `computeGradientAndHessian` | `IPC_FUNC.cpp` | ~263 | 提取 tet 弹性、三角形弹性、弯曲、障碍、摩擦为独立辅助函数 |
| `compute_fiction_hessian` | `IPC_FUNC.cpp` | ~244 | PP/PE/PT/EE 四块高度重复，模板化或统一实现 |
| `computeEnergyVal` | `IPC_FUNC.cpp` | ~193 | 提取各能量项 |
| `SpatialHash::calculateActivateSet` | `collisionUtil.cpp` | ~293 | 提取 PT/EE 各 case 辅助函数 |
| `SpatialHash::build(..., searchDir, ...)` | `collisionUtil.cpp` | ~212 | 两个重载共享体素网格计算 |
| `BHessian::toTriplets` | `IPCdistanceFuncs.cpp` | ~181 | 用统一泛型例程替换 D1/D2/D3/D4 四段重复 |
| `project_StabbleNHK_2_H_3D` | `fem3D.cpp` | ~139 | 删除 `if (false)` 死分支 |
| `load_triangleMesh` | `mesh.cpp` | ~188 | 拆分 OBJ 解析与网格维护 |

## 4. 魔法数字 / 硬编码常量

| 文件 | 位置 | 建议 |
|------|------|------|
| `Simulator.cpp` | `DefaultSettings` / `LoadSettings` | 集中为 `SimulationParameters` 结构体 |
| `IPC_FUNC.cpp` | 容差、迭代上限、CFL 因子 | 命名常量：`kPCGTolerance`、`kMaxCGIterations` 等 |
| `IPCtimeStepFuncs.cpp` | `computeInjectiveStepSize_3d(..., 1.0e-6, 0.2, ...)` | 命名参数 |
| `IPC_FUNC.cpp` | 重力 `-9.8` | 参数化 |
| `IPC_FUNC.cpp` | `suggestKappa` / `upperBoundKappa` 中 `1e-16`、`1e13` | 命名并加注释 |
| `IPCdistanceFuncs.cpp` | `eps_x = 1e-3 * edgeLen0² * edgeLen1²` | 命名系数 |
| `ACCD.cpp`、`IPCdistanceFuncs.cpp` | 距离类型返回码 0–8 | 改为 `enum class DistanceType` |

## 5. 头文件 / 包含结构

| 问题 | 文件 | 建议 |
|------|------|------|
| `using namespace std;` / `using namespace Eigen;` 在头文件中 | `mesh.h`、`SVD.h`、`ACCD.h`、`collisionUtil.h`、`Solver.h` | 头文件中禁止使用 using namespace |
| `#pragma once` 与 `#ifndef` 同时存在 | 多数头文件 | 统一用 `#pragma once` |
| `Solver.h` 直接包含 `<suitesparse/cholmod.h>` | `Solver.h` | 尽量前向声明或 PIMPL |
| `fem3D.h` 包含 `collisionUtil.h`、`fem_math.h` | `fem3D.h` | 尽量前向声明 |
| `FEMTimeIntegrator.h` 保留 CUDA 注释包含 | `FEMTimeIntegrator.h` | 删除注释 |
| `Solver.h` 包含 `cfloat` 但未使用 | `Solver.h` | 删除 |
| `fem_math.h` 依赖 `SVD.h` 仅取 `Vector3d` | `fem_math.h` | 直接包含 `<Eigen/Core>` |

## 6. 错误处理

| 问题 | 文件 | 建议 |
|------|------|------|
| 文件 I/O 使用 `exit(-1)` | `mesh.cpp` | 返回 `false` 或抛异常 |
| `printf`/`fprintf` 报错 | `mesh.cpp`、`IPC_FUNC.cpp` | 统一用 `std::cerr` 或异常 |
| `assert` 用于可能因输入失败的场景 | `ACCD.cpp`、`IPC_FUNC.cpp` | 输入校验改为错误处理；assert 只保留内部不变量 |
| `system("pause")` | `fem3D.cpp` | 删除 |
| `PCG_Solver` 中 `exit(0)` | `IPC_FUNC.cpp` | 返回状态码 |

## 7. const 正确性与传参

| 问题 | 文件 | 建议 |
|------|------|------|
| `Matrix2d` / `Matrix3d` 按值传递 | `SVD.h/cpp` | 改为 `const&` |
| `vec_double` 接受 `MatrixXd` 按值 | `fem_math.h/cpp` | 改为 `const MatrixXd&` |
| `distance(Vector3d x, Vector3d y)` 按值 | `fem3D.cpp` | 改为 `const Vector3d&` |
| `computePEPF_Aniostropic3D_double` 中 `Vector3d direction` 按值 | `fem3D.h/cpp` | 改为 `const Vector3d&` |
| `std::vector<Vector3d>` 多处按值传递 | `fem_math.cpp`、`IPC_FUNC.cpp` | 改为 `const&` |
| 求解函数应只读却非 const | `IPC_FUNC.cpp` | 审计并加 const |

## 8. 内存管理

| 问题 | 文件 | 建议 |
|------|------|------|
| 截图缓冲区用 `malloc`/`free` | `gl_main.cpp` | 用 `std::vector<unsigned char>` 或 `std::unique_ptr<[]>` |
| shader 加载用 `FILE*` 与固定 10KB 缓冲区 | `gl_main.cpp` | 用 `std::ifstream` / `std::string` |
| `CholmodSolver` 手动管理 CHOLMOD 指针 | `Solver.h/cpp` | 用 `std::unique_ptr` 自定义删除器或至少文档化所有权 |
| `BHessian::toTriplets` 先填充 dummy triplet | `IPCdistanceFuncs.cpp` | 直接 `push_back` 或 `reserve` |
| `BHessian` 每个牛顿迭代重新构造 | `IPC_FUNC.cpp` | 复用并 clear |
| `compute_H_dpt` 在 `spin_mutex` 下 push_back | `IPCdistanceFuncs.cpp` | 预 reserve 容量，避免 reallocation |
