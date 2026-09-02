# 11 — CHOLMOD supernodal + oneMKL 性能报告

本报告记录将CHOLMOD配置到本项目实测最佳状态的实现、许可边界与A/B结果。测试平台与`10_pardiso_report.md`相同：Windows/MSVC Release、Ryzen 9 9950X3D（16C/32T）、oneMKL 2025.2。所有结果使用当前方向Newton收敛判定、相同LBVH/CCD/CFL/严格能量语义。

## 1. 最终配置

- 固定 SuiteSparse 7.14.0 / CHOLMOD 5.3.5，启用 GPL supernodal 与 Partition/METIS；关闭 CUDA、CHOLMOD OpenMP、MatrixOps 和 Modify。
- BLAS/LAPACK：oneMKL static LP64，TBB threading；MKL 代码嵌入共享 `cholmod`，不依赖 OpenBLAS/LAPACK。SuiteSparse Config/AMD/CAMD/CCOLAMD/COLAMD 与 CHOLMOD 作为同一 bundle 构建。
- `common.supernodal=CHOLMOD_AUTO`；三个项目场景均实际选择supernodal。
- ordering固定AMD。Partition/METIS仍编译供实验，但强制METIS和multi-method AUTO在项目矩阵上都更慢。
- `--cholmod-threads 0`为默认自适应：矩阵`nnz<500,000`用4线程，否则8线程；正数显式覆盖。
- `CIPC_CHOLMOD_ROOT` 选择优化 prefix；CMake 强制全部 SuiteSparse 组件来自该 prefix，链接检查 `cholmod_super_numeric` 和 `cholmod_metis`，Windows 构建后复制六个 SuiteSparse DLL 与正确 TBB runtime。

推荐一键构建（依赖、CHOLMOD、连接、CPU-IPC Release全部完成）：

```bat
.\build.cmd
```

WSL/Ubuntu 22.04 对应入口为：

```bash
./build.sh
```

Windows 底层 `build_cholmod_mkl.ps1` 默认安装到 `build/cholmod-mkl-install`，主工程自动识别；Linux 脚本安装到隔离的 `build-wsl/cholmod-mkl-install` 并显式传入。两者读取固定 vcpkg revision 与 SuiteSparse 版本/校验和，vcpkg 只安装 Eigen、TBB core、METIS、oneMKL 和可选 FreeGLUT；SuiteSparse bundle 直接从源码构建，避免 vcpkg CHOLMOD 间接引入 OpenBLAS。内容签名允许后续构建直接复用 bundle。也可手动设置 `-DCIPC_CHOLMOD_ROOT=<prefix>`。Linux 的详细环境与复现记录见 `12_wsl_ubuntu_build.md`。

## 2. 许可边界

CHOLMOD supernodal 在当前 SuiteSparse 源码中属于 GPL-2.0-or-later 模块。性能 bundle 因此改变二进制分发的许可边界；不能把它描述为原先的非 GPL CHOLMOD core。PARDISO 与显式关闭性能要求的 system-CHOLMOD 诊断路径仍可独立构建。

## 3. Provider与算法A/B

bunny2，一个完整time step，3次独立进程中位数；所有运行得到相同终态：

| CHOLMOD配置 | step time | 相对最终配置 |
|---|---:|---:|
| supernodal + OpenBLAS 0.3.29 | 8.669 s | 5.33× slower |
| oneMKL + forced METIS | 2.557 s | 1.57× slower |
| oneMKL + AMD，未限制TBB | 2.464 s | 1.52× slower |
| oneMKL + AMD，16线程 | 2.202 s | 1.36× slower |
| oneMKL + AMD，12线程 | 2.041 s | 1.26× slower |
| oneMKL + AMD，8线程（最终大矩阵策略） | 约1.63 s | 1.00× |
| oneMKL + AMD，6线程 | 1.730 s | 1.06× slower |
| oneMKL + AMD，4线程 | 1.733 s | 1.07× slower |

结论：本项目中最大的收益来自“真正启用supernodal并让其内部直接绑定MKL”；只在应用末端补链MKL不会替换已构建CHOLMOD DLL的BLAS。METIS也不是无条件加速，当前矩阵由AMD胜出。8线程优于16线程，主要受MKL任务调度和缓存局部性影响。

## 4. 最终与PARDISO对比

当前固定依赖 bundle、5次独立进程中位数：

| 场景 | 测量步数 | CHOLMOD MKL | PARDISO 16 | 结论 |
|---|---:|---:|---:|---|
| bunny2 | 1 | 1.686 s | 1.310 s | PARDISO约28.7%更快 |
| cloth-bunny | 5 | 0.621 s | 0.575 s | PARDISO约8.0%更快 |
| twisting-mat | 5 | 0.412 s | 0.431 s | CHOLMOD约4.6%更快 |

三场景的终态、Newton与接触语义一致，只存在并行归约舍入差。PARDISO在大系统上仍有明确优势，因此保持默认；优化CHOLMOD已从“慢数倍的兼容后端”提升为中型场景上非常接近PARDISO的有效A/B后端。同一时段用旧单 DLL 与新完整 bundle 做 bunny2 五次对照，中位数为 `1.670/1.686 s`，差约 1%，说明去掉 vcpkg/OpenBLAS 安装和拆分 SuiteSparse 运行库没有实质改变 CHOLMOD 数值性能。

## 5. 验证

- 优化CHOLMOD Release与MSVC Debug均构建并运行；Debug通过DLL边界使用Release MKL/TBB runtime，不把静态MKL链接进Debug executable。
- `CIPC_CHOLMOD_ROOT=`的system/OpenBLAS兼容配置独立构建并通过smoke。
- Windows 产品目录自动包含 `suitesparseconfig/amd/camd/ccolamd/colamd/cholmod` 与 TBB DLL，不含 `openblas.dll` 或 `lapack.dll`；`cholmod.dll` 大小约 26.6 MB。
- Windows/WSL quadratic 与 non-quadratic Release、普通 Debug、PARDISO 关闭与 system-CHOLMOD fallback 均完成构建和 headless smoke；性能数字只来自 quadratic 默认路径。Linux non-quadratic 曾有的非正定现象来自 CRLF-sensitive Gmsh 拓扑解析，修复后两直接求解器均通过。
- WSL2 Ubuntu 22.04.5 使用 SuiteSparse 7.14.0/CHOLMOD 5.3.5 与 oneMKL 2025.2 完成 Release 全量构建；共享 bundle 通过 `--exclude-libs` 隐藏内嵌静态 MKL 符号，主程序链接无重复符号警告。默认 PARDISO 和显式 CHOLMOD 均完成 twisting-mat 单步 smoke，`ldd` 无 `not found` 或 OpenBLAS。
- Windows 从空的固定 vcpkg checkout 完成无 OpenBLAS/hwloc 的首次依赖安装、SHA-512 SuiteSparse 下载、bundle 构建、Release/Debug 产品编译和运行；第二次依赖运行命中内容签名，不重复配置 SuiteSparse。
