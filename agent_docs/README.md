# Agent Docs — CPU-IPC

面向接手开发/维护本项目的 AI agent 与开发者的实现文档。文档基于对源码的全量通读整理，所有结论均带 `文件:行号` 引用，修改代码前建议先核对行号（代码变动后可能漂移）。

## 项目是什么

**CPU-IPC**：单机 CPU 版 Incremental Potential Contact（IPC, Li et al. 2020）软体/布料仿真器。从 CUDA 代码移植而来（残留 `.cuh` 注释），使用：

- **Stable Neo-Hookean（SNK）** 四面体弹性模型（Smith et al. 2018）
- **Baraff-Witkin** 布料拉伸/剪切 + 二次弯曲（Bergou/Bridson 风格）
- 经典 IPC 对数障碍接触（PT/EE/PP/PE）+ 近平行 EE mollifier + 滞后摩擦
- Newton 迭代 + 能力感知默认后端（有 oneMKL 则 PARDISO，否则块感知 SuiteSparse LDL）/ 可选 CHOLMOD 与 Eigen-CG + 回溯线搜索 + Additive CCD
- TBB 并行；GLUT/OpenGL 固定管线可视化

## 文档索引

| 文件 | 内容 |
|---|---|
| `01_overview.md` | 构建/运行、整体架构、核心数据结构与全局约定 |
| `02_simulation_pipeline.md` | 每帧完整流程：solveIPCStep → barrier 子问题/Newton → 线搜索 → CCD → κ 更新 |
| `03_ipc_core.md` | 障碍函数、距离函数、`EncodedContact` 活动集编码、梯度/Hessian 装配 |
| `04_fem_model.md` | SNK 弹性、布料、弯曲、摩擦模型的公式与实现 |
| `05_collision_ccd.md` | SpatialHash、Ground、ACCD、可行步长 |
| `06_app_layer.md` | ViewerMain 可视化、Simulator 场景搭建、CholmodSolver、网格 IO、资产与脚本 |
| `07_gotchas.md` | **必读**：已知 bug、死代码清单、易踩的坑 |
| `08_optimization_roadmap.md` | 结构与性能优化路线、优先级、收益/风险和验证要求 |
| `09_optimization_report.md` | 已实施优化、实测结果、回归与仍保留的高风险事项 |
| `10_pardiso_report.md` | oneMKL PARDISO 实现、bunny2 大场景、线程/阶段/内存基准与 50 步正确性验证 |
| `11_cholmod_mkl_report.md` | GPL supernodal+oneMKL CHOLMOD 构建、ordering/线程A/B及与PARDISO的当前性能对比 |

## 30 秒上手

1. 装依赖（见 `01_overview.md`），`cmake -B build && cmake --build build --config Release`，得到 `cipc` viewer、`cipc_headless` 和 `cipc_core`。项目不再生成测试 target。
2. 运行 `cipc`：打开 1000×1000 GLUT 窗口，**空格键**开始/暂停仿真；每显示一帧 = 一个 IPC 时间步。资产和输出均使用编译期绝对路径，不再依赖当前工作目录。
3. `SimulationScene` 已接通三个场景：默认 `ClothOverBunny`、`TwistingMat`，以及参考 GPU_IPC 的双 `bunny2` 大场景（每只 scale=0.2）。场景参数在 `Assets/scene/parameterSetting.txt`，现按键名解析并校验未知、重复、缺失及越界值。
4. 时间步/Newton/线搜索编排在 `CPU IPC/IPCSolver.cpp`（`solveIPCStep` / `solveBarrierSubproblem`）；接触导数在 `ContactMechanics.cpp`，摩擦在 `Friction.cpp`，线性后端在 `NewtonLinearSystem.cpp`。
5. 无窗口运行：`cipc_headless --steps 20 --output Output/run`；PARDISO 大场景：`cipc_headless --scene bunny2 --steps 1 --linear-solver pardiso --pardiso-threads 16 --no-output`；批量基准：`python scripts/benchmark.py --exe <cipc_headless> --repeats 5 --steps 1`。这里一个 step 是完整仿真时间步，不是单次 Newton。

## 给 agent 的忠告

- 改动前先读 `07_gotchas.md` 和 `08_optimization_roadmap.md`，区分仍存在的问题、已修复记录和未接入的实验代码。
- 全项目**距离一律用平方值**，`Hhat` 是平方激活距离。比较/阈值写错数量级是最常见错误。
- 搜索方向是**减法**更新：`x ← x - α·searchDir`。所有 CCD 位移都用 `-searchDir`。

## 文档维护约定

- 每次修改代码时，同步更新受影响的 `agent_docs`；行为、路径、接口、默认值和已知问题不得只改代码不改文档。
- 每项性能优化必须记录基线、测试场景、指标与数值回归结果；没有测量时只能标为“静态推断”，不能宣称确定提速。
- 修复 `07_gotchas.md` 中的问题后，将其移入“已修复记录”，不要直接删除历史背景。
- 提交前至少执行 Release 构建、代表场景 headless smoke、`git diff --check`；修改弯曲编译分支时还应分别构建 quadratic ON/OFF。性能变化用固定参数的 `scripts/benchmark.py` 记录。
