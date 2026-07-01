# 数值方法与正确性关键代码

> 只读总结。以下不变量必须在整理、重构、优化过程中保持，否则数值结果会改变。

## 1. 时间积分

- **方法**：隐式欧拉（Implicit Euler）。
- **预测位置**：
  ```cpp
  xTilta[i] = V_prev[i] + v[i] * dt + gravity * dt * dt;
  ```
- **总能量目标**：
  ```
  E(x) = dt² [ E_elastic(x) + E_contact(x) + E_friction(x) ]
         + Σ ½ m_i || x_i − xTilta_i ||²
         + drag term ½ drag_coeff m_i || x_i − V_prev_i ||²
  ```

## 2. 非线性牛顿求解（`solve_subIP`，`IPC_FUNC.cpp`）

```cpp
for (k = 0; k < iterCap; ++k) {
    computeGradientAndHessian(mesh, gradient, BH, gd);
    if (distToOpt_PN < sqrt(1e-4 * bboxDiagSize2 * IPC_dt * IPC_dt))
        if (k > 0) break;
    calculateMovingDirection(mesh, BH, gradient, moveDir); // CHOLMOD
    Environment_largestFeasibleStepSize(...);
    Self_largestFeasibleStepSize(...);
    Self_largestFeasibleStepSize_CCD(...);
    lineSearch(mesh, sh, gd, moveDir, gradient, alpha, 0, 0, Kappa);
    postLineSearch(mesh, gd, alpha, Kappa);
}
```

- **收敛阈值**：`sqrt(1e-4 * bboxDiagSize2 * IPC_dt²)`。
- `calculateMovingDirection` 当前硬编码调用 `cholmod_solver`。
- `moveDir` 的符号与 `stepForward` 必须保持一致。

## 3. 梯度/海森组装（`computeGradientAndHessian`）

组装顺序：

1. 惯性 + drag
2. 四面体弹性（Stable Neo-Hookean）
3. 壳弹性（Baraff-Witkin）
4. 弯曲（`USE_QUADRATIC_BENDING` 使用预计算二次弯曲矩阵）
5. 自接触/地面接触障碍海森（`compute_H_dpt`、`compute_H_dee`）
6. 摩擦（`USE_FRICTION` 下 `compute_fiction_gradient` / `compute_fiction_hessian`）

### 关键存储：`BHessian`

```cpp
class BHessian {
    vector<int>       D1Index;
    vector<Vector3i>  D3Index;
    vector<Vector4i>  D4Index;
    vector<Vector2i>  D2Index;
    vector<Matrix<double,12,12>> H12x12;
    vector<Matrix<double, 3, 3>> H3x3;
    vector<Matrix<double, 6, 6>> H6x6;
    vector<Matrix<double, 9, 9>> H9x9;
    vector<Triplet<double>> toTriplets(const vector<int>& Btype);
};
```

- 所有海森块在插入前调用 `IglUtils::makePD` 强制 PSD。
- `D1/D2/D3/D4Index` 与对应尺寸矩阵的顺序必须严格同步。

## 4. IPC 障碍公式（`IPCdistanceFuncs.cpp`）

### 障碍函数

```cpp
b(d)  = -(d - dHat)² log(d / dHat)
g_b(d) = (d - dHat) log(d / dHat) * -2 - (d - dHat)² / d
H_b(d) = -2 log(d / dHat) - 4(d - dHat)/d + (d - dHat)² / d²
```

### 距离原语

- `d_PP`, `d_PE`, `d_PT`, `d_EE`：平方距离。
- `g_*`, `H_*`：符号生成的导数，**不要手改**。
- 自接触海森：
  ```cpp
  IPHessian = (coef * H_b) * g * g.transpose() + (coef * g_b) * H;
  makePD(IPHessian);
  ```
- 边-边接触带 mollifier `e(x, eps_x)`：
  ```cpp
  PEEHessian = kappa_gradb_gradeT + kappa_gradb_gradeT.transpose()
             + (coef * b) * e_H
             + (coef * e * H_b) * grad_d * grad_d.transpose()
             + (coef * e * g_b) * H_d;
  ```

## 5. 连续碰撞检测（`ACCD.cpp`）

- `point_triangle_ccd` / `edge_edge_ccd`。
- `slackness = 0.8` → `CCDDistRatio = 0.2`。
- **搜索方向传入 CCD 时取反**：
  ```cpp
  point_triangle_ccd(..., -searchDir[...], ..., CCDDistRatio, 0);
  ```
- 线搜索会重新计算碰撞集并保证 `!isIntersected`。

## 6. 线搜索（`lineSearch`）

```cpp
bool lineSearch(mesh3D&, SpatialHash&, const Ground&,
    const vector<Vector3d>& searchDir, const vector<Vector3d>& gradient,
    double& stepSize, double armijoParam, double lowerBound, double Kappa)
```

- **重要**：`armijoParam` 被无条件覆盖为 `0`，因此线搜索退化为**简单能量下降测试**：
  ```cpp
  while (testingE > lastEnergyVal + stepSize * c1m) ...
  ```
- 先检查 `isIntersected(...)` 并折半直到无相交。
- 再检查能量下降，直到 `stepSize > 1e-3 * LFStepSize`。
- 若步长缩小，再次检查相交。

> 注意：当前行为看起来像 bug，但修改它会改变求解轨迹。重构时应默认保留。

## 7. 接触刚度 κ 自适应

- `initKappa`：平衡弹性梯度与接触梯度。
  ```cpp
  double minKappa = -gsum / gsnorm;
  ```
- `postLineSearch`：当约束距离小于 `mesh.dTol` 时，`Kappa` 翻倍。
- `upperBoundKappa`：用 `meanMass`、`Hhat`、`bboxDiagSize2` 限制 κ 上限。

## 8. 摩擦模型（`FrictionUtils.hpp`）

- 默认使用 **C1 光滑静态摩擦**：
  ```cpp
  #define SFCLAMPING_ORDER 1
  ```
- 切向基、最近点、相对切向位移、Hessian 块 lifting。
- 摩擦距离/面积：`Fhat = 1e-6 * bboxDiagSize2`。
- 摩擦求解使用 `fricDHat = Fhat * IPC_dt * IPC_dt`。

## 9. 弹性与 SVD

### Stable Neo-Hookean

```cpp
getObjEnergy_StableNHK2_3D:
    QRSVD(F) -> sigma, V
    S = V * sigma * V.transpose()
    I2 = (S*S).trace()
    I3 = S.determinant()
    energy = 0.5*lengthRate*(I2-3)
           + 0.5*volumeRate*(I3 - 1 - 3*lengthRate/(4*volumeRate))²
           - 0.5*lengthRate*log(I2+1)
```

### SVD 约定

`QRSVD` 内部对 `U`/`V` 做符号翻转使其为旋转矩阵，奇异值矩阵可能含一个负值。  
**不要改变该约定**，弹性通过 `S = V * Sigma * V.transpose()` 重建。

## 10. 参数与单位

运行时参数从 `Assets/scene/parameterSetting.txt` 加载（`LoadSettings`）。关键值：

| 参数 | 典型值 |
|------|--------|
| 体积密度 | 1000 |
| 体积杨氏模量 | 1e4 |
| 泊松比 | 0.49 |
| 摩擦系数 | 0.4 |
| 布料厚度 | 0.001 |
| 布料杨氏模量 | 5e4 |
| 弯曲杨氏模量 | 1e7 |
| `ipc_time_step` | 0.01 |
| `Newton_solver_threshold` | 1e-2 |
| `IPC_ralative_dHat` | 1e-3 |

`Hhat` 缩放：
```cpp
Hhat = IPC_ralative_dHat;   // e.g. 1e-3
Hhat *= Hhat;               // 1e-6
Hhat *= bboxDiagSize2;      // 绝对平方距离
```

## 11. 编译宏配置

当前激活：

- `USE_SNK`
- `USE_FRICTION`
- `USE_QUADRATIC_BENDING`

**未激活**：`NEWB`（因此使用原始对数障碍 `compute_g_dpt` / `compute_H_dpt`，而非 `compute_g_dpt_new` / `compute_H_dpt_new`）。

## 12. 输出与回归验证

- 每帧表面 OBJ：`saveSurface/surf_XXXXX.obj`
- 截图：`saveScreen/step_XXXXX.bmp`（当前关闭）
- 计时日志：`timeCost.txt`、`tempData/timeCost.txt`
- 每 10 帧检查点：`tempData/vertex.txt`、`tempData/vertexXtile.txt`

### 建议的回归测试

1. 对比相同场景下 `saveSurface/` 的顶点位置（紧容差）。
2. 对比 `timeCost.txt` 中的迭代数、碰撞数、总能量。
3. 每步后检查无负体积（`calculateVolum`）。
4. 线搜索后检查无自相交（`isIntersected`）。
5. 对比 κ 演进与 `Hhat`/`Fhat` 值。

## 13. 重构风险清单

| 风险点 | 影响 |
|--------|------|
| 修改 `armijoParam = 0` 行为 | 改变线搜索轨迹 |
| 改变 `BHessian` index 与矩阵的同步关系 | 海森组装错误 |
| 手改符号生成的 `g_*` / `H_*` | 导数错误 |
| 切换 `NEWB` / 未定义宏 | 使用不同障碍公式 |
| 改变 SVD 符号约定 | 弹性能量/梯度/海森错误 |
| 切换 CG/PCG 替代 CHOLMOD | 数值结果改变 |
| 全局日志/计时变量 `step_index` 等 | 重启/日志状态需保留 |
