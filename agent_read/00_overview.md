# CPU-IPC 项目总览

> 生成时间：2026-07-01  
> 项目路径：`c:\Users\KemengHuang\Desktop\kimi-gipc\cipc\CPU-IPC`  
> 当前分支：`refactor`  
> 语言：C++14

## 1. 项目目标

CPU-IPC 是一个基于 **增量势接触（Incremental Potential Contact, IPC）** 的 CPU 有限元仿真器，支持：

- 可变形体（四面体，Stable Neo-Hookean）
- 布料/壳（三角形，Baraff-Witkin 膜 + 弯曲）
- 自接触与地面对象接触
- 摩擦力
- OpenGL 可视化

## 2. 顶层数据流

```
gl_main.cpp
    └── FEMSimulator::buildModels(...) / simulateStick(...)
            └── ImplicitFEMIntegrator::integrate(...)
                    └── IPC_Solver(...)
                            ├── computeGradientAndHessian(...)      ← 构建梯度/海森
                            ├── calculateMovingDirection(...)       ← CHOLMOD 求解
                            ├── Self/Environment_largestFeasibleStepSize(...) ← CCD
                            ├── lineSearch(...)                     ← 线搜索
                            └── postLineSearch(...)                 ← κ 自适应
```

## 3. 主要文件职责

| 文件/模块 | 职责 |
|----------|------|
| `gl_main.cpp` | GLUT + OpenGL 入口、渲染、截图、键盘鼠标交互 |
| `Simulator.h/cpp` | 场景构建、参数加载、模型管理 |
| `FEMTimeIntegrator.h/cpp` | 时间积分器接口，当前为隐式欧拉 |
| `IPC_FUNC.h/cpp` | IPC 求解器核心：牛顿迭代、线搜索、κ 自适应、摩擦、线性求解器封装 |
| `IPCdistanceFuncs.h/cpp` | 距离函数、障碍函数、符号生成的梯度/海森、CCD 步长 |
| `IPCtimeStepFuncs.h/cpp` | 体积翻转保护、地面 CCD、二次/三次方程求根 |
| `fem3D.h/cpp` | Stable Neo-Hookean、Baraff-Witkin 壳、弯曲的能量/梯度/海森 |
| `mesh.h/cpp` | `mesh3D` 数据结构、网格加载（`.msh` / `.obj`）、表面提取 |
| `collisionUtil.h/cpp` | 空间哈希（SpatialHash）、地面（Ground）、活动集计算 |
| `ACCD.h/cpp` | 仿射不变 CCD（point-triangle / edge-edge） |
| `FrictionUtils.hpp` | 摩擦力切向基、相对位移、光滑静态摩擦模型 |
| `Solver.h/cpp` | CHOLMOD 稀疏直接求解器封装 |
| `SVD.h/cpp`, `SVD/` | 3×3 SVD 封装（Implicit QR SVD） |
| `fem_math.h/cpp` | 向量/矩阵工具、多项式求根 |
| `fem_parameters.h` | 全局常量命名空间 `FEM` |
| `mIPC.h/cpp` | `MMCVID` 约束标识符 |

## 4. 外部依赖

- **Eigen3**：稠密/稀疏线性代数
- **Intel TBB**：`parallel_for`、`parallel_reduce`、`spin_mutex`
- **SuiteSparse / CHOLMOD**：稀疏直接求解器
- **OpenGL + GLUT + GLEW**：可视化
- **OpenBLAS**：CHOLMOD/Eigen 底层 BLAS（可选）

## 5. 构建系统

- `CMakeLists.txt`：单目标可执行文件 `cipc`
- 默认 `Release` 配置
- 编译定义（CMake 注入）：
  - `CIPC_ASSETS_DIR`
  - `CIPC_OUTPUT_DIR`
  - `USE_SNK`
  - `USE_FRICTION`
  - `USE_QUADRATIC_BENDING`

## 6. 代码规模

`CPU IPC/` 目录下约 13,000+ 行 C++ 代码，核心集中在：

- `IPCdistanceFuncs.cpp`：~3,900 行（障碍与距离导数）
- `IPC_FUNC.cpp`：~2,300 行（求解器主循环）
- `fem3D.cpp`：~1,800 行（弹性核函数）
- `mesh.cpp`：~1,150 行（网格数据结构）
- `collisionUtil.cpp`：~1,050 行（空间哈希与碰撞检测）
- `ACCD.cpp`：~600 行

## 7. 当前无自动化测试

- 没有单元测试、回归测试或基准测试框架
- 验证方式：OpenGL 可视化、`saveSurface/` 输出、`timeCost.txt` 计时、每 10 帧的 `tempData/` 检查点
- 建议后续引入检查点对比作为回归验证
