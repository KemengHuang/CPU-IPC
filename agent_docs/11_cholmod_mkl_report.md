# 11 — CHOLMOD supernodal + oneMKL 性能报告

本报告记录将CHOLMOD配置到本项目实测最佳状态的实现、许可边界与A/B结果。测试平台与`10_pardiso_report.md`相同：Windows/MSVC Release、Ryzen 9 9950X3D（16C/32T）、oneMKL 2025.2。所有结果使用当前方向Newton收敛判定、相同LBVH/CCD/CFL/严格能量语义。

## 1. 最终配置

- CHOLMOD 5.3.0，启用GPL supernodal与Partition/METIS模块；关闭CUDA、CHOLMOD OpenMP、MatrixOps和Modify。
- BLAS/LAPACK：oneMKL static LP64，TBB threading；MKL代码嵌入共享`cholmod.dll`，不依赖OpenBLAS/LAPACK DLL。
- `common.supernodal=CHOLMOD_AUTO`；三个项目场景均实际选择supernodal。
- ordering固定AMD。Partition/METIS仍编译供实验，但强制METIS和multi-method AUTO在项目矩阵上都更慢。
- `--cholmod-threads 0`为默认自适应：矩阵`nnz<500,000`用4线程，否则8线程；正数显式覆盖。
- `CIPC_CHOLMOD_ROOT`选择优化prefix；CMake链接检查`cholmod_super_numeric`和`cholmod_metis`，Windows构建后复制优化DLL与Release TBB runtime。

构建：

```powershell
powershell -ExecutionPolicy Bypass -File scripts/build_cholmod_mkl.ps1 `
  -VcpkgRoot D:/VCPKG/vcpkg
cmake -S . -B build -DCMAKE_TOOLCHAIN_FILE=D:/VCPKG/vcpkg/scripts/buildsystems/vcpkg.cmake
```

脚本默认安装到`build/cholmod-mkl-install`，主工程自动识别；也可显式传`-DCIPC_CHOLMOD_ROOT=<prefix>`。

## 2. 许可边界

CHOLMOD supernodal在当前SuiteSparse/vcpkg中属于GPL-2.0-or-later模块。性能版DLL因此改变二进制分发的许可边界；不能把它描述为原先的非GPL CHOLMOD core。PARDISO、LDL与system CHOLMOD兼容路径仍可独立构建。

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

当前代码、5次独立进程中位数：

| 场景 | 测量步数 | CHOLMOD MKL | PARDISO 16 | 结论 |
|---|---:|---:|---:|---|
| bunny2 | 1 | 1.625 s | 1.301 s | PARDISO约24.9%更快 |
| cloth-bunny | 5 | 0.608 s | 0.586 s | PARDISO约3.6%更快 |
| twisting-mat | 5 | 0.413 s | 0.434 s | CHOLMOD约5.0%更快 |

三场景的终态、Newton与接触语义一致，只存在并行归约舍入差。PARDISO在大系统上仍有明确优势，因此保持默认；优化CHOLMOD已从“慢数倍的兼容后端”提升为中型场景上非常接近PARDISO的有效A/B后端。

## 5. 验证

- 优化CHOLMOD Release与MSVC Debug均构建并运行；Debug通过DLL边界使用Release MKL/TBB runtime，不把静态MKL链接进Debug executable。
- `CIPC_CHOLMOD_ROOT=`的system/OpenBLAS兼容配置独立构建并通过smoke。
- 优化DLL依赖中含`tbb12.dll`，不含`openblas.dll`或`lapack.dll`；大小约26.6MB。
- quadratic/non-quadratic Release、普通Debug、PARDISO关闭与system-CHOLMOD fallback均完成构建和headless smoke；性能数字只来自quadratic默认路径。
