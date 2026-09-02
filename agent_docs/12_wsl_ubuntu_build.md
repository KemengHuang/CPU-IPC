# 12 — WSL / Ubuntu 构建与配置

本页记录 CPU-IPC 的 Linux 原生构建入口及在当前 WSL 环境中的端到端验证。目标不是从 WSL 调用 Windows 的 `.exe` 或复用 Windows 库，而是用 GCC、`x64-linux` vcpkg、Linux oneMKL 和 `libcholmod.so` 生成一套独立产品。

## 1. 已验证环境

- WSL2：Ubuntu 22.04.5 LTS，kernel `6.18.33.1-microsoft-standard-WSL2`
- 编译器：GCC/G++ 11.4.0
- 生成工具：CMake 4.3.3、Ninja 1.10.1
- vcpkg triplet：`x64-linux`
- 实际安装版本：oneMKL 2025.2.0、SuiteSparse 7.14.0 / CHOLMOD 5.3.5、TBB 2023.1.0、Eigen 5.0.1，以及当前 vcpkg baseline 的 METIS

脚本明确限制在 `x86_64`，因为当前生产直接求解器依赖 oneMKL。用户通常不需要预先执行 apt 命令：`build.sh` 会检查编译器、Ninja、Git、curl、tar 和 pkg-config，缺少时自动请求一次 sudo 并安装精确缺项。CMake 要求至少 3.23；系统版本过旧或缺失时，vcpkg 会下载私有工具，脚本自动选择它。若要禁止脚本修改系统（例如 CI 或受管工作站），传 `--no-system-packages`；脚本会退出并打印等价的命令，例如：

```bash
sudo apt-get update
sudo apt-get install -y build-essential ninja-build git curl tar pkg-config
```

当前已验证的 WSL 实例已经具备主体工具，不需要再次安装。`zip`/`unzip` 若缺失，脚本优先用 `apt-get download` 将两个 `.deb` 解包到用户缓存，不执行 `sudo`，也不修改系统包状态。默认 viewer 构建还会检测 FreeGLUT 需要的 OpenGL/X11 开发包并按需安装；只有 `--headless-only` 跳过这项检查。

## 2. 一键构建

在仓库根目录运行：

```bash
./build.sh
```

viewer 与 headless 是 Linux 脚本的默认产品。该命令依次完成：

1. 自动补齐必要的 Ubuntu 基础工具；不要求用户理解数值库的系统包名。
2. 优先使用 `--vcpkg-root` / `VCPKG_ROOT` 显式选择；否则自动 bootstrap `scripts/vcpkg-revision.txt` 固定的 revision 到 `${XDG_CACHE_HOME:-$HOME/.cache}/cpu-ipc/vcpkg`，不受 PATH 上其他 vcpkg 影响。
3. 在一次 vcpkg transaction 中为 `x64-linux` 安装 Eigen、`tbb[core]`、METIS、oneMKL 和默认 viewer 所需的 FreeGLUT；`--headless-only` 时不安装 FreeGLUT。`tbb[core]` 明确排除项目未使用的 hwloc/TBBBind。
4. 按 `suitesparse-version.txt` 下载固定 SuiteSparse 7.14.0 并验证 SHA-512，只构建 Config、AMD、CAMD、CCOLAMD、COLAMD 与 GPL supernodal CHOLMOD，把 BLAS/LAPACK 绑定到 oneMKL static LP64/TBB，并安装到 `<repo>/build-wsl/cholmod-mkl-install`。不安装 OpenBLAS。
5. 以 `CIPC_ENABLE_PARDISO=ON`、`CIPC_REQUIRE_OPTIMIZED_CHOLMOD=ON` 配置主工程；默认同时生成 `<repo>/build-wsl/cpu-ipc/cipc` viewer 与 `cipc_headless`，服务器可传 `--headless-only`。

PARDISO 可用时仍是 Release 默认后端，线程默认 16；优化 CHOLMOD 是自动 fallback，也可用 `--linear-solver cholmod` 显式选择。Eigen-CG 仍只作为显式比较项。

## 3. 路径与平台隔离

仓库不保存开发者用户名、盘符或固定 vcpkg 目录。默认路径均从脚本位置、`HOME`/`XDG_CACHE_HOME` 和用户参数计算：

| 内容 | Windows | WSL / Ubuntu |
|---|---|---|
| 用户入口 | `.\build.cmd` | `./build.sh` |
| 底层入口 | `build.ps1` | `build.sh` |
| 自动 vcpkg | `%LOCALAPPDATA%/CPU-IPC/vcpkg-<revision>` | `${XDG_CACHE_HOME:-$HOME/.cache}/cpu-ipc/vcpkg` |
| 优化 CHOLMOD | `<repo>/build/cholmod-mkl-install` | `<repo>/build-wsl/cholmod-mkl-install` |
| 主工程 | `<repo>/build/cpu-ipc` | `<repo>/build-wsl/cpu-ipc` |
| headless 产品 | `<repo>/build/cpu-ipc/Release/cipc_headless.exe` | `<repo>/build-wsl/cpu-ipc/cipc_headless` |
| viewer 产品 | `<repo>/build/cpu-ipc/Release/cipc.exe` | `<repo>/build-wsl/cpu-ipc/cipc` |

Windows 与 Linux build tree 不能互换；两边的 CMake cache、库格式和 ABI 都不同。WSL 可以直接从 `/mnt/c` 完成构建，本次验证就是这样执行的；但大量小文件的编译和 vcpkg 安装通常在 Linux 文件系统（例如 `~/src/CPU-IPC`）更快。

已有 Linux vcpkg 时可以覆盖默认值：

```bash
./build.sh --vcpkg-root /path/to/vcpkg
```

也可使用 `export VCPKG_ROOT=/path/to/vcpkg`。不要把 Windows vcpkg checkout 传给 WSL。

## 4. CHOLMOD 的 Linux 连接方式

一键脚本不再安装 vcpkg 的 CHOLMOD 包，因为它会为源码供应间接拉入本项目不用的 OpenBLAS。项目以固定 URL、版本和 SHA-512 直接取得 SuiteSparse 源码，通过顶层 CMake 只启用 `cholmod`（自动带出五个必要组件），并明确设置：

- `CHOLMOD_GPL=ON`、`CHOLMOD_SUPERNODAL=ON`、`CHOLMOD_PARTITION=ON`
- CUDA、CHOLMOD OpenMP、MatrixOps、Modify 关闭
- oneMKL static、LP64 接口、TBB threading layer
- 六个同 prefix 的共享库；`libcholmod.so` 内部包含其 MKL 实现

主程序同时直接链接 oneMKL PARDISO。若 SuiteSparse 共享库导出内嵌静态 MKL 符号，GNU linker 会看到两份实现并产生重复符号警告；Linux bundle 使用 `-Wl,--exclude-libs,ALL` 隐藏内部静态归档符号。重新链接后警告消失，CHOLMOD 和主程序各自仍正常使用同一 oneMKL provider。Windows 会自动把 `suitesparseconfig/amd/camd/ccolamd/colamd/cholmod` 六个 DLL 与正确 TBB runtime 复制到每个产品目录。

bundle 安装目录包含内容签名；SuiteSparse 版本、vcpkg/关键 ports、provider/shim 或构建 recipe 均未变化且文件齐全时，第 3 阶段直接复用，不再重复配置和安装。本机二次 `--dependencies-only` 验证中，该阶段由约十余秒降为即时命中。

启用 supernodal 会把该 CHOLMOD 二进制纳入 GPL-2.0-or-later 分发边界，Windows 与 Linux 一致。

## 5. 可选参数

```text
--vcpkg-root PATH   使用已有 Linux vcpkg checkout 或 vcpkg 可执行文件
--build-dir PATH    覆盖默认 <repo>/build-wsl
--project-build-dir PATH
                    仅覆盖 CPU-IPC build tree，复用默认依赖 bundle
--config TYPE       Release、RelWithDebInfo 或 Debug
--viewer            构建 GLUT viewer（默认）
--headless-only     只构建 benchmark/headless 产品
--nonquadratic-bending
                    构建完整二面角 hinge bending 分支
--dependencies-only 只安装/构建全部库，不编译 CPU-IPC
--no-system-packages
                    不调用 sudo/apt，只报告缺失的 Ubuntu 包
```

默认构建会安装 `freeglut:x64-linux`，并要求系统具有 OpenGL/X11 开发库；脚本在缺失时会列出或自动安装对应的 Ubuntu `apt` 包。`--viewer` 是保留的显式同义开关，`--headless-only` 才跳过这些图形依赖。当前 WSL 已实际完成 viewer 配置、编译、链接，并在 WSLg 下成功启动 3 秒。普通无桌面的 Ubuntu 仍需可用的 X server/`DISPLAY` 才能显示窗口；这不影响 headless benchmark。

## 6. 实机验证结果

零参数 `./build.sh` 已从系统前置检查、固定 vcpkg 依赖安装、校验 SuiteSparse 下载、优化 bundle 编译/安装一直执行到 CPU-IPC GCC Release viewer/headless 链接成功；随后再次增量运行同一命令也成功。`./build.sh --dependencies-only --no-system-packages` 也验证通过，只准备依赖后正常退出；第二次运行命中 bundle 内容签名并即时复用。CMake cache 确认 `CIPC_ENABLE_PARDISO=ON`、`CIPC_REQUIRE_OPTIMIZED_CHOLMOD=ON`，编译定义包含 `CIPC_HAS_PARDISO`。

运行命令：

```bash
./build-wsl/cpu-ipc/cipc_headless \
  --scene twisting-mat --steps 1 --no-output

./build-wsl/cpu-ipc/cipc_headless \
  --scene twisting-mat --steps 1 --linear-solver cholmod --no-output

./build-wsl/cpu-ipc/cipc_headless \
  --scene twisting-mat-soft --steps 5 --no-output
```

前两条 hard 路径都完成 1 个完整时间步、2 次 Newton，`matrix_nnz=151335`；`squared_norm_sum` 均为 `560.74908434861868`，坐标和只存在约 `1e-14` 的并行浮点归约差异。soft 范例的 PARDISO/CHOLMOD 五步均为末步 1 Newton、`matrix_nnz=158907`、`squared_norm_sum=560.77059938845662`。`ldd` 能解析优化 `libcholmod.so.5`，没有 `not found`。这里的单次耗时只用于 smoke，不是隔离性能基准；正式 PARDISO/CHOLMOD A/B 仍应使用 `scripts/benchmark.py` 多进程重复并报告中位数。

零参数 `./build.sh` 的 viewer 路径已完整通过：CMake 找到 Linux OpenGL、X11 与 vcpkg FreeGLUT，生成并链接 `build-wsl/cpu-ipc/cipc`；默认 `twisting-mat-soft` viewer 已在 WSLg 启动并持续运行 3 秒后由验证命令主动终止。显式 `--viewer` 得到相同配置。

Gmsh 2.2 tet loader 已改用与换行类型无关的字段解析。cloth-bunny 在 Windows/WSL 都得到 `7356 tets / 2893 vertices / 4074 surface faces`；non-quadratic 模式下 PARDISO 与 CHOLMOD 均完成单步，和 Windows 同为 8 Newton、12 次能量回退、`nnz=161979`，不再出现假性非正定。

项目按要求没有测试源码或 CTest target；Linux 回归以 Release 构建、PARDISO/CHOLMOD headless smoke、viewer 启动 smoke、运行时依赖检查和物理/迭代指标一致性为准。

## 7. Windows 同等的一键体验

Windows 用户不再需要设置执行策略或输入 PowerShell 长命令，直接运行：

```bat
.\build.cmd
```

该包装会调用 `build.ps1`，检查 64-bit MSVC 的 Desktop development with C++ workload；缺失时打印可直接复制的 winget 安装命令。CMake 会依次从 PATH、Visual Studio 自带工具和 vcpkg 下载缓存解析。自动 vcpkg 固定 revision 并放在 `%LOCALAPPDATA%` 短路径，避免深层仓库触发 Windows path-length 问题；只有用户显式传参/环境变量才覆盖。随后 Eigen、`tbb[core]`、METIS、oneMKL 与可选 FreeGLUT 在一次 vcpkg transaction 中安装，固定 SuiteSparse bundle 和主工程继续自动完成。`.\build.cmd -DependenciesOnly` 只准备库，`.\build.cmd -HeadlessOnly` 跳过 FreeGLUT/viewer；两条分支均已从空 checkout 验证。
