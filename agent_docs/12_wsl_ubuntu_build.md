# 12 — WSL / Ubuntu 构建与配置

本页记录 CPU-IPC 的 Linux 原生构建入口及在当前 WSL 环境中的端到端验证。目标不是从 WSL 调用 Windows 的 `.exe` 或复用 Windows 库，而是用 GCC、`x64-linux` vcpkg、Linux oneMKL 和 `libcholmod.so` 生成一套独立产品。

## 1. 已验证环境

- WSL2：Ubuntu 22.04.5 LTS，kernel `6.18.33.1-microsoft-standard-WSL2`
- 编译器：GCC/G++ 11.4.0
- 生成工具：CMake 4.3.3、Ninja 1.10.1
- vcpkg triplet：`x64-linux`
- 实际安装版本：oneMKL 2025.2.0、SuiteSparse 7.14.0 / CHOLMOD 5.3.5、TBB 2023.1.0、Eigen 5.0.1，以及当前 vcpkg baseline 的 METIS

脚本明确限制在 `x86_64`，因为当前生产直接求解器依赖 oneMKL。全新 Ubuntu 若缺少基础编译工具，先执行一次：

```bash
sudo apt-get update
sudo apt-get install -y build-essential cmake ninja-build git curl tar
```

当前已验证的 WSL 实例已经具备这些工具，不需要再次安装。`zip`/`unzip` 若缺失，脚本会用 `apt-get download` 将两个 `.deb` 解包到用户缓存，不执行 `sudo`，也不修改系统包状态。

## 2. 一键构建

在仓库根目录运行：

```bash
bash build.sh --headless-only
```

`headless` 是 Linux 脚本的默认值，因此 `bash build.sh` 等价。该命令依次完成：

1. 按 `--vcpkg-root`、`VCPKG_ROOT`、`PATH` 的优先级查找 Linux vcpkg；均不存在时自动 bootstrap 到 `${XDG_CACHE_HOME:-$HOME/.cache}/cpu-ipc/vcpkg`。
2. 为 `x64-linux` 安装 Eigen、TBB、METIS、oneMKL 和 CHOLMOD 的 core/Partition 依赖。
3. 从 vcpkg 已提取的同版本 SuiteSparse 源码单独构建 GPL supernodal CHOLMOD，并把 BLAS/LAPACK 直接绑定到 oneMKL static LP64/TBB。
4. 将优化 CHOLMOD 安装到 `<repo>/build-wsl/cholmod-mkl-install`。
5. 以 `CIPC_ENABLE_PARDISO=ON`、`CIPC_REQUIRE_OPTIMIZED_CHOLMOD=ON`、`CIPC_BUILD_VIEWER=OFF` 配置主工程，生成 `<repo>/build-wsl/cpu-ipc/cipc_headless`。

PARDISO 可用时仍是 Release 默认后端，线程默认 16；优化 CHOLMOD 是自动 fallback，也可用 `--linear-solver cholmod` 显式选择。Eigen-CG 仍只作为显式比较项。

## 3. 路径与平台隔离

仓库不保存开发者用户名、盘符或固定 vcpkg 目录。默认路径均从脚本位置、`HOME`/`XDG_CACHE_HOME` 和用户参数计算：

| 内容 | Windows | WSL / Ubuntu |
|---|---|---|
| 根入口 | `build.ps1` | `build.sh` |
| 自动 vcpkg | `<repo>/build/_deps/vcpkg` | `${XDG_CACHE_HOME:-$HOME/.cache}/cpu-ipc/vcpkg` |
| 优化 CHOLMOD | `<repo>/build/cholmod-mkl-install` | `<repo>/build-wsl/cholmod-mkl-install` |
| 主工程 | `<repo>/build` | `<repo>/build-wsl/cpu-ipc` |
| headless 产品 | `<repo>/build/Release/cipc_headless.exe` | `<repo>/build-wsl/cpu-ipc/cipc_headless` |

Windows 与 Linux build tree 不能互换；两边的 CMake cache、库格式和 ABI 都不同。WSL 可以直接从 `/mnt/c` 完成构建，本次验证就是这样执行的；但大量小文件的编译和 vcpkg 安装通常在 Linux 文件系统（例如 `~/src/CPU-IPC`）更快。

已有 Linux vcpkg 时可以覆盖默认值：

```bash
bash build.sh --vcpkg-root /path/to/vcpkg --headless-only
```

也可使用 `export VCPKG_ROOT=/path/to/vcpkg`。不要把 Windows vcpkg checkout 传给 WSL。

## 4. CHOLMOD 的 Linux 连接方式

vcpkg 的 CHOLMOD core/Partition 安装只用于提供同版本源码和 SuiteSparse 子库；一键脚本重新配置 CHOLMOD，并明确设置：

- `CHOLMOD_GPL=ON`、`CHOLMOD_SUPERNODAL=ON`、`CHOLMOD_PARTITION=ON`
- CUDA、CHOLMOD OpenMP、MatrixOps、Modify 关闭
- oneMKL static、LP64 接口、TBB threading layer
- 共享 `libcholmod.so`，内部包含其 MKL 实现

主程序同时直接链接 oneMKL PARDISO。若 `libcholmod.so` 导出其内嵌静态 MKL 符号，GNU linker 会看到两份实现并产生重复符号警告；脚本对该共享库使用 `-Wl,--exclude-libs,ALL` 隐藏内部静态归档符号。重新链接后警告消失，CHOLMOD 和主程序各自仍正常使用同一 oneMKL provider。

启用 supernodal 会把该 CHOLMOD 二进制纳入 GPL-2.0-or-later 分发边界，Windows 与 Linux 一致。

## 5. 可选参数

```text
--vcpkg-root PATH   使用已有 Linux vcpkg checkout 或 vcpkg 可执行文件
--build-dir PATH    覆盖默认 <repo>/build-wsl
--config TYPE       Release、RelWithDebInfo 或 Debug
--viewer            同时构建 GLUT viewer
--headless-only     只构建 benchmark/headless 产品（默认）
```

`--viewer` 会额外安装 `freeglut:x64-linux`，并要求系统具有 OpenGL/X11 开发库；vcpkg 在缺失时会列出对应的 Ubuntu `apt` 包。当前 WSL 已实际完成 viewer 配置、编译、链接，并在 WSLg 下成功启动 3 秒。普通无桌面的 Ubuntu 仍需可用的 X server/`DISPLAY` 才能显示窗口；这不影响 headless benchmark。

## 6. 实机验证结果

`bash build.sh --headless-only` 已从 vcpkg 依赖安装、优化 CHOLMOD 配置/编译/安装一直执行到 CPU-IPC GCC Release 链接成功；随后再次增量运行同一命令也成功。CMake cache 确认 `CIPC_ENABLE_PARDISO=ON`、`CIPC_REQUIRE_OPTIMIZED_CHOLMOD=ON`，编译定义包含 `CIPC_HAS_PARDISO`。

运行命令：

```bash
./build-wsl/cpu-ipc/cipc_headless \
  --scene twisting-mat --steps 1 --no-output

./build-wsl/cpu-ipc/cipc_headless \
  --scene twisting-mat --steps 1 --linear-solver cholmod --no-output
```

两条路径都完成 1 个完整时间步、2 次 Newton，`matrix_nnz=151335`；`squared_norm_sum` 均为 `560.74908434861868`，坐标和只存在约 `1e-14` 的并行浮点归约差异。`ldd` 能解析优化 `libcholmod.so.5`，没有 `not found`。这里的单次耗时只用于 smoke，不是隔离性能基准；正式 PARDISO/CHOLMOD A/B 仍应使用 `scripts/benchmark.py` 多进程重复并报告中位数。

`bash build.sh --viewer` 也已完整通过：CMake 找到 Linux OpenGL、X11 与 vcpkg FreeGLUT，生成并链接 `build-wsl/cpu-ipc/cipc`；WSLg 启动 smoke 持续运行 3 秒后由验证命令主动终止。

项目按要求没有测试源码或 CTest target；Linux 回归以 Release 构建、PARDISO/CHOLMOD headless smoke、viewer 启动 smoke、运行时依赖检查和物理/迭代指标一致性为准。
