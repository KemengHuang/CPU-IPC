# 05 — 碰撞检测：CPU LBVH、SpatialHash、Ground、ACCD

涉及文件：`CPU IPC/LBVH.cpp`、`CollisionBroadPhase.cpp`、`ContactMechanics.cpp`、`AdditiveCCD.cpp`、`FeasibleStep.cpp` 和可选 `TetInversionGuard.cpp`。默认宽阶段为 `BroadPhaseBackend::LinearBVH`，SpatialHash 保留为 `--broad-phase spatial-hash` 对照后端。

## 1. CPU Linear BVH（默认）

实现参考同仓库 `GPU_IPC/GPU_IPC/mlbvh.cu/.cuh`，保持其核心算法与数据布局语义：

1. face 与 edge 分别建立一棵树；叶 AABB 由三角形/边端点组成，Full CCD 叶盒同时覆盖 `x` 与 `x−α·searchDir`。
2. 合并全部叶盒得 scene AABB；叶盒中心归一化到 scene 后生成 30-bit Morton code。
3. key = `(morton << 32) | primitiveId`，即使中心/Morton 相同也保证唯一。
4. 排序后用 Karras `determineRange + findSplit` 建立 `2N−1` 节点的 LBVH：internal `[0,N−2]`，leaf `[N−1,2N−2]`。
5. CPU 端用显式 postorder 栈自底向上合并 internal AABB；查询使用固定栈遍历。
6. PT：点/扫掠点盒查询 face tree；EE：边/扫掠边盒查询 edge tree。结果按 primitive id 排序，再使用与 SpatialHash 共用的精确静态/扫掠 AABB 二次判定、共享顶点、全驱动和顺序过滤后才进入 partial/full ACCD。

与 GPU 版的有意差别：CPU 的 AABB overlap 包含刚好接触的边界以保持保守性；Full CCD 仍按本 CPU 求解器 thickness=0 的既有语义使用零 padding，而 GPU 查询额外传 `sqrt(dHat)` 作为更大的候选 superset。

历史验证：曾以 brute force 和不同 PT/EE sweep α 逐 pair 比较两后端，并完成 20 步轨迹审计；相关结果保留在 `09_optimization_report.md`，测试源码现已移除。

## 2. SpatialHash（回归后端，`CollisionBroadPhase.cpp`）

均匀体素网格 + `std::unordered_map<int, std::vector<int>> voxel`（只占用的格子存在）。线性 id：`ix + iy·voxelCount[0] + iz·voxelCount0x1`（x 最快）。

**格宽**：所有调用点都传 `mesh.averageEdgeLenth`（= 平均表面边长 / 3，`Simulator.cpp:317-322`）。与 dHat 无直接关系——`sqrt(Hhat)` 只作**查询半径**。

**索引编码**（`:98-99`）：
- `id < surfEdgeStartInd (= surfVerts.size())` → 表面顶点 svI
- `[surfEdgeStartInd, surfTriStartInd)` → 表面边 `id − surfEdgeStartInd`
- `≥ surfTriStartInd` → 表面三角形 `id − surfTriStartInd`

**两种构建**：

1. **静态 `build(mesh, voxelSize)`**（`:32-176`）：边界 = `surfVerts` 的精确 min/max（**无 margin**）；顶点进单格；边/三角形插入其端点体素 index 的整数 min/max 盒覆盖的所有格子（并行算、串行合并）。`pointAndEdgeOccupancy` **不填**。
2. **扫掠 `build(mesh, searchDir, curMaxStepSize, voxelSize, use_V_prev)`**（`:183-395`）：为 CCD 服务。终点 `SVt = V − curMaxStepSize·searchDir`；网格边界覆盖起点+终点整段轨迹；每个基元的插入盒是其线性扫掠的保守包围盒；**填 `pointAndEdgeOccupancy`**（仅顶点+边，不含三角形）——occupancy 类查询只在此构建后可用（`queryPointForTriangles(svI)` 内有 assert）。

退化保护：voxelCount 溢出时塌缩成单格（`:88-96`）。

**查询函数一览**（全部裁剪到网格范围）：

| 函数 | 位置 | 说明 |
|---|---|---|
| `queryPointForTriangles(pos, radius, out)` | 当前实现 | 点±半径盒 → 三角形（排序 vector 去重；LBVH/Hash 共用接口） |
| `queryPointForTriangles(svI, out)` | `:427` | occupancy 版（需扫掠构建），供全 CCD |
| `queryTriangleForPoints / ...ForEdges` | `:441/:475` | 三角形 AABB±半径 → 顶点/边 |
| `queryEdgeForEdgesWithBBoxCheck(mesh, v0, v1, r, out, eIq)` | `:539` | 活动集用：`seJ > eIq` 排序去重 + 精确静态 AABB 重叠测试 |
| `queryEdgeForEdgesWithBBoxCheck(mesh, searchDir, step, seI, out)` | `:507` | occupancy 版 + 扫掠 AABB 测试（每条边 4 角点） |
| `queryEdgeForPE` | `:587` | 边 AABB，**无半径**，同时回顶点和边 |
| `queryPointForPrimitives` | `:627` | 线段 pos→pos+dir 覆盖盒 → 三类基元 |

**当前优化**：表面局部索引、静态/扫掠 voxel range、primitive occupancy 和 point/edge occupancy 均复用成员缓冲；voxel map 不再每次销毁，而只清空上一轮 active buckets，并在累计空 key 过多时重建。主 PT/EE 查询使用线程局部排序 vector，替代逐 query `unordered_set`。静态/扫掠 build 仍有重复逻辑，`spanSize` 钳制仍未启用。

## 3. Ground

平面 `n·x = D`，默认 `n=(0,1,0), D=-1`。`calculateGapFromObj` 返回**平方** gap `(n·x[vId]−D)²`。`calculateActivateSet`：gap² < Hhat 的表面顶点（全局 id）进 `Environment_ActiveSet`。地面摩擦不存切基，解析投影 `VProj = VDiff − (VDiff·n)n`。

## 4. Additive CCD（`AdditiveCCD.cpp`，Li et al. 2021）

`point_triangle_ccd`（`:315-359`）、`edge_edge_ccd`（`:259-313`）。**无多项式求根**：每次迭代按可证线性下界步进。

算法（以 EE 为例）：
1. 减去平均位移 `mov = Σdᵢ/4`（平移不变）。
2. `max_disp_mag`：EE = `√max(|dea0|²,|dea1|²) + √max(|deb0|²,|deb1|²)`；PT = `|dp| + max|dtᵢ|`。为 0 → 返回 1.0。
3. 初始 `dFunc = dist² − thickness²`（EE 退化/已在厚度内时回退到 4 对端点-端点距离最小值）；目标 `gap = eta·(dist − thickness)`。
4. 循环：`toc_lb = (1−eta)(dist−thickness)/max_disp_mag`，全点前进 `toc_lb·d`，重算距离；`toc>0` 后当剩余 reduced distance < `gap` 即停；`toc>1` → 返回 1.0。
5. **返回 toc ∈ [0,1]** = 撞击时间占整段位移的比例（即步长缩放）；1.0 = 步内无碰撞。现有 10000 次保险上限，并对非有限/非正推进量返回当前保守 toc，防止退化输入死循环。

调用点参数（`ContactMechanics.cpp`）：`eta = CCDDistRatio = 1 − slackness = 0.2`，`thickness = 0`，位移传 `-searchDir`。

`point_triangle_ccd_broadphase / edge_edge_ccd_broadphase`（`:161-201`）：AABB 粗筛，**定义了但从未被调用**；注意其 `d*` 参数按“终点绝对位置”解释（与窄阶段的“位移”语义不同），若启用需传 `x_end = x + dx`。

## 5. 可行步长函数

- `Self_largestFeasibleStepSize`（部分 CCD）：遍历后端无关的 `Self_CCD_ActiveSet`，解码 primitive pair，ACCD 取 min；已移除无效 `SpatialHash/candidates` 参数。
- `Self_largestFeasibleStepSize_CCD`（全 CCD）：查询扫掠候选后，再用与测试共用语义的起终点 AABB overlap 过滤，跳过全驱动/共享/重复 pair，ACCD 取 min。PT 先在入口区间求得安全 α；EE 随后使用该 PT 更新后的 α 做统一扫掠 AABB 过滤。ACCD 仍接收完整 `-searchDir`，所以返回值保持原始全局 `[0,1]` 参数，无需再次乘 PT α。SpatialHash/LBVH 都经过相同 EE 精确过滤，保留后端等价性。
- `limitStepByGround`（`FeasibleStep.cpp`）：地面射线-平面解析上限；`normalMotion<0` 时取 `−distance/normalMotion·slackness`。
- `limitStepToPreventTetInversion`（`TetInversionGuard.cpp`）：保留的可选 tet 防翻路径。利用 determinant 的多线性构造三次体积多项式，求体积降到指定比例的最小正根；默认主流程不调用。
- 原 `AdditiveCCD.cpp` 内两套未调用且参数语义含糊的 broad-phase 备份已删除；统一候选入口只有 `CollisionBroadPhase/LBVH`。

## 6. 穿透检测（线搜索防护）

- 地面：`isIntersected` 检查 gap ≤ 0。
- 自碰撞：`checkEdgeTriIntersectionIfAny` 用当前 edge broad phase 查询 triangle AABB，再以 `IglUtils::segTriIntersect` 精确测试；发现相交则步长减半。线搜索与动画边界回退均有 64 次上限和 underflow 失败诊断。
