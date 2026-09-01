# 04 — FEM 材料模型与摩擦

涉及文件：`CPU IPC/Elasticity.cpp`（逐单元内核）、`RotationAwareSVD.cpp/.h` + `ImplicitQR3x3SVD.h`、`FrictionKinematics.h`、`Friction.cpp/.h`，以及 `IPCSolver.cpp` 中的总装配端。材料常数由 `Simulator.cpp` 的 `updateMaterial` 统一换算；仅保存 π 和注释的 `fem_parameters.h` 已删除。

```
μ_Lamé = E/(2(1+ν))      λ_Lamé = Eν/((1+ν)(1−2ν))
lengthRate (μ) = 4·μ_Lamé/3        volumeRate (λ) = λ_Lamé + 5·μ_Lamé/6   (Smith et al. 2018 换算, α ≡ 1 + 3μ/4λ)
stretchStiffness = clothE/(1−ν²)
shearStiffness = shearE/(2(1+ν)·t)
plateRigidity D = bendE·t³/(12·(1−ν²))
```

## 1. Stable Neo-Hookean（四面体，活跃模型）

- **能量**（`getObjEnergy_StableNHK2_3D`，`Elasticity.cpp`，TBB parallel_reduce；`F = Ds·Dm⁻¹`，`computeRotationAwareSVD(F)`，`I2 = Σσᵢ²`，`I3 = Πσᵢ = J`）：
  ```
  ψ = ½μ(I2 − 3) + ½λ(I3 − 1 − 3μ/4λ)² − ½μ·log(I2 + 1)
  ```
  逐 tet × `volum[i]` 求和。`getObjRestEnergy_StableNHK2_3D`（`:820`）= 常数静止能量 `(½λ(3μ/4λ)² − ½μ·log4)·vol`（`restSNKE` 记账用）。
- **PK1**（`computePEPF_StableNHK3D_2_double`，`:1185-1211`，求解器实际调用）：
  ```
  P = μ(1 − 1/(I2+1))·F + λ(J − 1 − 3μ/4λ)·cof(F)
  ```
  `cof(F)` 逐元素手写在 `pI3pF`（`:1197-1207`）。
- **Hessian 投影**（`project_StabbleNHK_2_H_3D`，`:1376-1515`）：活跃的是 `else` 分支（`:1472-1513`）——闭式 9×9 Hessian：
  ```
  H = μ(1 − 1/(I2+1))·I₉ + λ(J−1−3μ/4λ)·HJ + μ/(2(I2+1)²)·g1·g1ᵀ + λ·gJ·gJᵀ
  ```
  `g1 = vec(2F)`，`gJ = vec(cof F)`（列叉积），`HJ` 由三列的反对称矩阵 `f0hat/f1hat/f2hat` 拼成；最后 **`IglUtils::makePD<double,9>` 负特征值钳零**（`:1512`）。`if(false)` 的解析特征系分支被禁用。
- 单元到全局：梯度 `f = volum·PFPXᵀ·vec(P)`，Hessian `volum·dt²·PFPXᵀ·Hq·PFPX`（装配端 `IPCSolver.cpp`）。

**PFPX**（∂vec(F)/∂x）：tet 为固定尺寸 9×12，布料为 6×9。两者只依赖静止 `Dm⁻¹`，现于 `initMesh3D` 预计算进 `tetrahedraPFPX/trianglePFPX`，Newton 内核以 const reference 使用；PK1 向量通过 Eigen Map 读取列主序数据，不再构造动态 `MatrixXd vec()`。

## 2. Baraff-Witkin 布料（三角形，拉伸/剪切）

材质坐标架由 `calculateDms2D_double`（`:615-645`）定义：把静止三角形用 Rodrigues 公式旋转到法线对齐 (0,1,0)，取 (x,z) 为 UV；各向异性方向固定 `a=(1,0), b=(0,1)`。不变量（3×2 F）：`I5u=‖F·a‖², I5v=‖F·b‖², I6=aᵀFᵀFb`。

- **能量**（`getObjEnergy_baraffwitkin_3D`，`:776-818`；`ℓu=‖Fa‖, ℓv=‖Fb‖`）：
  ```
  stretch = (ℓu−1)² + [ℓu>1]·strainRate·(ℓu−1)³ + (ℓv−1)² + [ℓv>1]·strainRate·(ℓv−1)³
  shear   = I6²
  E = (stretchS·stretch + shearS·shear)·areas[i]     // areas 含厚度 = 体积
  ```
  三次 `strainRate` 项**只在拉伸时**激活（tension-only 应变硬化，`:796-801`）。
- **PK1**（`computePEPF_baraffwitkin_double`，`:1214-1246`）：`P = stretchS·(2ucoeff·Faaᵀ + 2vcoeff·Fbbᵀ) + shearS·2(I6−a·b)(Fabᵀ+Fbaᵀ)`；`ucoeff = 1 − 1/√I5u`（拉伸时 `+1.5·strainRate·(√I5u + 1/√I5u − 2)`）。
- **Hessian**（`project_baraffwitkint_H_3D`，`:1517-1580`）：6×6，**构造即 PSD**（无钳制）——拉伸项对角块 + rank-1 `fu·fuᵀ`；剪切项用解析特征系（置换矩阵 H、`λ0 = ½(I2+√(I2²+12I6²))`、`T = ½(I+sign(I6)H)`）：`H_shear = 2[|I6|(T − TqTqᵀ/‖Tq‖²) + λ0·q0q0ᵀ]`。

## 3. 弯曲

**活跃路径（`USE_QUADRATIC_BENDING`）**——flat-shell 二次弯曲，常量 4×4 Hessian 预计算于 `prepareQuadraticBending`（`Simulator.cpp` 内部场景辅助函数）：

```
K = (c03+c04, c01+c02, −c01−c03, −c02−c04)ᵀ   // cot 权重, cot(e,f)=e·f/|e×f|
Q = 3/(a0+a1) · K·Kᵀ                          // a0/a1 = 两三角静止面积
```

按 `QuadBendingInfo{verts,Q,hessianBase=Q⊗I₃}` 存于 `mesh.quadBendingInfo`。能量 `½·D·xᵀhessianBase·x`；Newton 中局部梯度用一次固定 12×12 matvec，Hessian 直接缩放 `hessianBase`。Q 天生 PSD，已删除每轮通用 eigendecomposition。

**可选路径（关闭 `USE_QUADRATIC_BENDING`，`HingeBending.cpp`）**——完整离散 hinge：每条 interior edge 预计算

```
l0 = ‖x1−x0‖
h0 = 2A0/l0,  h1 = 2A1/l0
w_e = l0/(h0+h1) = l0²/[2(A0+A1)]
θ0 = edgeTheta(x_rest)
E_e = ½·D·w_e·(θ−θ0)²
```

角差用 `remainder(θ−θ0,2π)` 取最近分支。梯度为 `D·w_e·Δθ∇θ`，精确 Hessian 为 `D·w_e[∇θ∇θᵀ+Δθ∇²θ]`，进入 Newton 前再做 PSD 投影。能量、梯度、Hessian 共用 `HingeBendingInfo{vertices,restAngle,geometricWeight}`，不再逐轮重算静止长度/角度。实现时曾以有限差分核对一阶/二阶导数及均匀缩放不变性；测试源码现已按项目要求移除。

## 4. 摩擦（IPC 滞后摩擦，`Friction.cpp/.h` + `FrictionKinematics.h`）

- `SFCLAMPING_ORDER = 1`（`FrictionKinematics.h`）→ **C1 光滑钳制**：
  `‖u‖ < εv` 时 `f0 = x²(−√x²/3 + ε)/ε² + ε/3`，`f1/‖u‖ = (−√x² + 2ε)/ε²`，`f2 = 2(ε−√x²)/ε²`（`x = ‖u‖²`）；C0/C2 变体也在文件内（`:260-305`），分发器 `f0_SF/f1_SF_div_relDXNorm/f2_SF`（`:308-339`）。
- **过渡速度**：`εv = sqrt(Fhat·dt²)`（`IPCSolver.cpp` 传 `eps2 = mesh.Fhat * IPC_dt²`）。
- **每步冻结量**（`Friction::initialize`）：
  - λ：`λ = −κ·b'(d)·2√d`（PP/PE × 重数；EE 带符号翻转）→ `Self_lambda_lastH / Environment_lambda_lastH`；
  - 几何：最近点坐标 β/γ/η → `MMDistCoord`；3×2 切基（PT/EE 用边方向+双叉积，PE 用边方向+`v12×(v0−v1)`，PP 选较优的 `UnitX×v01`/`UnitY×v01`）→ `MMTanBasis`；
  - 活动集快照 `*_activeSet_lastH`（摩擦用**时间步开始**的活动集）。
- **梯度**（`Friction::addGradient`）：相对滑动 `relDX = Tᵀ·Δx_closest`（`computeRelDX_*` + `liftRelDXTanToMesh_*`，权重见下表）；动摩擦方向归一，静摩擦乘 `f1/‖u‖`；地面摩擦直接把 `x − V_prev` 投影到平面。
- **Hessian**（`Friction::addHessian`）：冻结 λ/最近点/切基后的 2D 特征系已解析化，不再形成“PSD矩阵减 rank-one”后调用通用 eig-project：
  - 滑动：平行滑动特征值 0，垂直特征值 `μλ/‖u‖`；只装配垂直方向 lifted rank-one。
  - C1 静摩擦：平行特征值 `μλ·f2`，垂直特征值 `μλ·f1`；零滑移用切空间 X/Y 两个 rank-one。
  - 地面 C0 clamp 同样以解析地面切基构造；所有项都是非负特征值乘 rank-one PSD 和，无 `makePD`。
  - 解析式在实现阶段曾以静/动态/零滑移样本核对旧闭式与 PSD；当前不保留测试 target。

各接触类型的 TTT 度量矩阵权重（`computeTTT_*`，`FrictionKinematics.h`）：

| 类型 | 每节点权重（⊗ 切基ᵀ） |
|---|---|
| PT（12×12） | `[1, −1+β1+β2, −β1, −β2]` |
| EE（12×12） | `[1−γ1, γ1, γ2−1, −γ2]` |
| PE（9×9） | `[1, η−1, −η]` |
| PP（6×6） | `[1, −1]` |

## 5. 旋转保持 3×3 SVD

- FEM 统一调用 `computeRotationAwareSVD(Matrix3d)`；公开结果明确命名为 `leftSingularVectors/singularValueMatrix/rightSingularVectors`。
- 底层 `ImplicitQR3x3SVD.h` 是自包含的 Gast-Fu-Jiang-Teran implicit-shifted symmetric QR 算法；原 `SVD/Tools.h` 的少量数学模板已经内联进去，不再单独建目录/辅助文件。
- σ 按绝对值降序、仅末项可负、U/V 都是 proper rotation。实现阶段曾核对零/退化/翻转/随机矩阵；当前不保留测试 target。
- 三个无人调用的 Eigen JacobiSVD 包装已删除。

## 6. 初始化与静止状态（`initMesh3D`，`Elasticity.cpp`）

- 逐 tet：`Dm`（`:575-603`）、体积（`calculateVolum`，`:540-566`，`|det|/6`）、集中质量 `V·ρ/4` 到 4 节点、`Dm.inverse()` 存入。
- 逐 tet 同时预计算固定尺寸 `TetPFPX`；逐布料三角形计算 `calculateDms2D_double`、**`areas = 面积 × clothThicness`**（体积化！）、质量和 `TrianglePFPX`。相关 vector 会按单元数 reserve。
- `type`/`scale` 参数**被忽略**。
- 其余网格伴随数据在别处建：`surface/surfVerts/surfEdges`、`quadBendingInfo` 或 `hingeBendingInfo`。固定/驱动状态只由 `boundaryTypes` 表示。

## 7. 死代码清单（本模块）

- `computePEPF_StableNHK3D_double`（v1，`:1159`）与 `project_StabbleNHK_H_3D`（解析特征系，`:1284`，内有遗留 `printf` `:1325`）——v1 模型整体未被调用。
- 各向异性肌肉模型四件套（`computePEPF_Aniostropic*`、`project_ANIOSI5*`，`:1248-1694`）——无调用方。
- `Iso*` 等距嵌入实验函数群（`:162-425, 947-1084`）、`__Inverse/__Inverse2x2`（`:427-537`）、`calculateDms3D2_double` 等。
