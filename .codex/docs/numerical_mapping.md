# MAS1998 数值复现、代码映射与验证状态

本文说明 Spee et al. (1998), *A numerical study for global atmospheric
transport-chemistry problems* 的数值模型在独立 MOOSE app 中如何实现，以及当前证据能够支持
什么复现结论。

## 1. 先看结论

复现代码位于：

```text
/home/fangxiaozhong/git_repo/mas1998_benchmark
```

论文和化学附录分别为：

```text
/mnt/d/ZoteroData/storage/D7JMCSVH/1-s2.0-S0378475498001554-main.pdf
/mnt/d/BaiduSyncdisk/Zotero/ZoteroData/storage/DQYC9G5U/05028D.pdf
```

当前状态必须分成三层理解：

| 层次 | 状态 | 能够得出的结论 |
| --- | --- | --- |
| 模型与算法实现 | 已完成 | 约化网格、平流、垂直扩散、化学、两类分裂、质量修正、六组 Table 3 配置和独立特征线参考生成器都有可执行实现 |
| 组件与短时全链路验证 | 已完成 | 守恒、CFL、BDF2、Gauss-Seidel、坐标/轨迹均有单元测试；Type I、Type II 和特征线参考均有 600 s smoke test |
| 14 天论文结果的定量复现 | 待长算例 | 生产/参考输入、误差范数、收敛门槛和批量运行器已就绪；六个 14 天结果及参考收敛数据尚未实际生成 |

因此，准确表述是：**已经完成论文模型、数值方法和 Section 4.5 参考解生成方法的可执行重实现，
但在 14 天任务及收敛检查结束前，仍不能声称已经定量复现论文 Figs. 4-8。** 新参考输出是
`INDEPENDENT_CHARACTERISTIC_REFERENCE`，不是尚未找回的作者 CWI 数组。

MAS1998 专用逻辑全部由独立 app 拥有。MOOSE
`modules/atmospheric_chemistry` 只保留通用 KPP 能力：加载 `Fun_SPLIT`、计算 production/loss、
设置固定物种。模块中没有 MAS1998 网格、常数、物种表、时间步或诊断逻辑。

## 2. 方程、未知量与代码数据布局

论文式 (2.1)-(2.5) 在当前假设下写为

```text
dc/dt = F0(c) + F1(c) + F2(c)

F0(c) = R(c) = P(c) - L(c)c
F1(c) = d_z [rho K d_z(c/rho)]
F2(c) = -1/(a cos(phi)) [d_lambda(uc) + d_phi(vc cos(phi))]
```

其中：

- `c` 是物种数浓度，单位 `molecules cm^-3`；
- `rho` 是空气数密度，`c/rho` 是 mixing ratio；
- `a = 6378 km`；
- 只有水平风平流和垂直湍扩散，没有垂直平流、水平扩散或地形；
- 17 个变量物种都参与水平平流和垂直扩散。

`MAS1998BenchmarkSolver` 将状态按下面的顺序展平：

```text
(horizontal_cell, vertical_layer, KPP_variable_species)
```

索引由 `MAS1998BenchmarkSolver::stateIndex()` 计算。水平单元数取决于约化网格，垂直固定为
15 层，变量物种固定为 17 个。完整运行不使用普通 libMesh 网格表达论文的约化纬圈拓扑；
MOOSE 的一维占位 `Mesh` 只用于建立 app 和 `Executioner`，实际状态及时间推进由数组求解器持有。

主要调用链是：

```text
MAS1998Executioner::execute()
  -> MAS1998BenchmarkSolver::solve()
     -> advanceTypeI() 或 advanceTypeII()
        -> advanceAdvection()
        -> advanceColumns()
           -> MAS1998ChemistryIntegrator::integrate()
              -> MAS1998VerticalOperator::solveImplicit()
              -> KPPGeneratedMechanism::computeProductionLoss()

MAS1998ReferenceExecutioner::execute()
  -> MAS1998CharacteristicReferenceSolver::solve()
     -> solidBodyTrajectory() 反向求初始脚点、正向给出每一步位置
     -> MAS1998ChemistryIntegrator::integrate()
        -> MAS1998VerticalOperator::solveImplicit()
        -> KPPGeneratedMechanism::computeProductionLoss()
```

## 3. 算子分裂：式 (3.2) 和式 (3.3)

为避免字母 `A` 同时表示“平流”和论文中的扩散矩阵，本文使用 `Adv`、`Dif`、`Chem`。

### 3.1 Type I：所有过程分开

论文式 (3.2a-e) 是对称 Strang 序列：

```text
Adv(dt_split/2)
Dif(dt_split/2)
Chem(dt_split)
Dif(dt_split/2)
Adv(dt_split/2)
```

对应 `MAS1998BenchmarkSolver::advanceTypeI()`。每个水平位置提取一根
`15 layers x 17 species` 的柱：

- `Dif` 调用 chemistry integrator 时设置 `chemistry=false, diffusion=true`；
- `Chem` 设置 `chemistry=true, diffusion=false`；
- 化学不耦合垂直层，氮修正在每个高度 box 内分别进行。

### 3.2 Type II：化学与垂直扩散耦合

论文式 (3.3a-c) 为：

```text
Adv(dt_split/2)
(Chem + Dif)(dt_split)
Adv(dt_split/2)
```

对应 `MAS1998BenchmarkSolver::advanceTypeII()`。中间子步设置
`chemistry=true, diffusion=true`，一个物种的 15 个高度值在同一个隐式三对角系统内求解；
氮修正也对层厚加权后的整根柱统一进行。

### 3.3 外层与内层步数

构造 solver 时强制检查

```text
dt_split = 2 * dt_adv
```

所以每个平流半步恰好包含一个二阶段 RK 步。六个 Type I 生产输入中
`dt_dif = dt_split/4`，故每个扩散半步包含两个隐式步：首步 implicit Euler，第二步 BDF2。
化学或耦合子步使用生产输入指定的固定 `300 s` 内步长。

`MAS1998ChemistryIntegrator::integrate()` 每进入一个新的 split 子问题都会重置 BDF 历史，
因此每个子问题的第一步都是 implicit Euler；这对应论文在每个 split interval 起点的处理。

## 4. 水平约化网格

实现类为 `MAS1998ReducedGrid`。每个纬圈保存纬度边界、中心、经度单元数、经度间距和全局
数组偏移。经度覆盖 `[-pi, pi)`，`normalizeLongitude()` 保证周期接缝闭合。

### 4.1 Table 2 的规模

下表的纬圈序列只列南极到赤道一侧，随后关于赤道镜像。例如 `12 x 64` 表示连续 12 个
纬圈各有 64 个经度单元。

| 名义网格 | 半球纬圈经度数 | 3D cells | 17 物种 unknowns |
| --- | --- | ---: | ---: |
| `64x32x15` | `4, 8, 16, 32, 12 x 64` | 24,840 | 422,280 |
| `128x64x15` | `4, 8, 16, 32, 8 x 64, 20 x 128` | 93,960 | 1,597,320 |
| `256x128x15` | `4, 8, 16, 32, 64, 17 x 128, 42 x 256` | 391,560 | 6,656,520 |

`MAS1998ReducedGrid.reproducesTable2` 对论文 Table 2 的 cell 和 unknown 总数作精确断言。

证据边界：`64x32` 的纬圈变化可从论文 Fig. 1 直接读出；论文只给出另外两组的 reduction
次数和总单元数，没有逐纬圈表。`128x64` 与 `256x128` 序列是满足对称、单调、二次幂变化
并精确匹配 Table 2 总数的重建，不能称为原始 CWI 程序逐项恢复。

### 4.2 不同分辨率纬圈的连接

在相邻纬圈经度数不同的界面，`MAS1998ReducedGrid::interfaceSegments()` 采用较细一侧的面段数：

```text
n_segments = max(n_lon_south, n_lon_north)
```

每个面段的经度中点分别映射到南、北真实单元。粗单元因而可与多个细面段连接，等价于论文
Fig. 1 所述的 piecewise-constant virtual concentration。一个真实公共面段只生成一次，后续
通量以相反符号加入两侧单元。

`MAS1998ReducedGrid.closesSeamAndMapsReductionInterfaces` 检查周期接缝、南北对称、面段范围
及粗细连接。

## 5. 水平平流：式 (3.5)-(3.13)

实现类为 `MAS1998TransportOperator`。对每个高度和每个物种独立推进一个水平数组。

### 5.1 守恒面通量

代码使用与论文式 (3.5) 一致的球面中点测度：

```text
cell_measure = a^2 cos(phi_center) d_lambda d_phi
longitude_face_length = a d_phi
latitude_segment_length = a cos(phi_face) d_lambda_segment
```

东西面通量为 `u * c_face * longitude_face_length`，南北面段通量为
`v * c_face * latitude_segment_length`。每个通量只计算一次，然后从上游侧减去、向下游侧
加入，所以按 `cell_measure` 加权后内部通量严格成对抵消。

纬度循环只组装内部纬圈界面。极点处没有额外水平边界面，和球面极点面长为零一致。经度方向
用模运算闭合 `-pi/pi` 接缝。

测试 `MAS1998TransportOperator.conservesIntegratedMass` 验证

```text
sum_h cell_measure[h] * RHS[h] = 0
```

达到浮点误差量级。

### 5.2 三阶迎风重构和 limiter

`limitedFaceValue(far_upwind, upwind, downwind)` 实现论文式 (3.7)-(3.9)：

```text
theta = (upwind - far_upwind) / (downwind - upwind)
psi(theta) = max(0, min(1, theta, 1/3 + theta/6))
c_face = upwind + psi(theta) * (downwind - upwind)
```

当分母接近零时直接返回 `upwind`，避免舍入误差放大。在约化界面跨纬圈寻找
`far_upwind` 时，`cellAtLongitude()` 给出 piecewise-constant 虚拟值；最靠近极点且不存在
更外纬圈时取 `far_upwind = upwind`，即边界侧零斜率重构。

测试 `limiterUsesPublishedFormula` 直接检查上述代数式。

### 5.3 二阶段显式梯形 RK

`MAS1998TransportOperator::advance()` 对应式 (3.10)-(3.11)：

```text
w       = c_n + dt * f(c_n)
c_(n+1) = c_n + dt/2 * [f(c_n) + f(w)]
```

风场本身与时间无关，但 limiter 使半离散算子依赖当前浓度，因此第二阶段仍需重新计算通量。
更新后只允许舍入误差量级的负数并截为零；更大的负值会报错，而不是静默修改结果。

`preservesConstantsAndSeamPulse` 检查常数场、跨周期缝 pulse 的质量守恒和非负性。

### 5.4 CFL

`maximumCFL()` 对应式 (3.12)-(3.13)，并针对约化界面做了保守扩展：先累计一个真实单元上
所有细面段的绝对南北通量率，再与东西通量率相加。solver 要求

```text
max_cell(dt_adv * (longitude_rate + latitude_rate)) <= 2/3
```

`paperStepsMeetCFLBound` 分别检查 `64x32: 2400 s`、`128x64: 1200 s` 和
`256x128: 600 s`。

## 6. 垂直扩散与式 (2.5) 零通量边界

实现类为 `MAS1998VerticalOperator`。这一节给出完整离散，避免把“反射 ghost center”和
“ghost concentration”混为一谈。

### 6.1 15 层网格

论文式 (3.15) 后给出的中心为：

```text
0.3, 1.0, 2.2, 4.3, 6.5, 8.4, 10.0, 11.3,
13.0, 15.2, 17.6, 19.8, 22.5, 27.6, 34.7 km
```

底边界 `z=0`，顶边界 `z_H=38.2 km`。为计算首末控制体宽度，代码只构造几何 ghost center：

```text
z_0     = -z_1                         = -0.3 km
z_(N+1) = 2*z_H - z_N                  = 41.7 km
Delta_z_k = (z_(k+1) - z_(k-1)) / 2
```

于是首个 ghost 与首中心的中点正好是 `0 km`，末中心与末 ghost 的中点正好是 `38.2 km`。
这只确定控制体宽度；代码没有创建或求解 ghost concentration。

### 6.2 内部面离散

定义 mixing ratio `q_k = c_k/rho_k`，内部面的扩散通量为

```text
J_(k+1/2) = (rho*K)_(k+1/2)
            * (q_(k+1) - q_k) / (z_(k+1) - z_k)
```

其中 `(rho*K)_(k+1/2)` 在两个中心的几何中点重新计算，而不是在压力层边界取值。这正是论文
式 (3.16) 为非均匀网格指定的做法。控制体内的扩散倾向是

```text
(A_z c)_k = [J_(k+1/2) - J_(k-1/2)] / Delta_z_k
```

`MAS1998VerticalOperator` 将其直接组装为 `_lower`、`_diagonal`、`_upper` 三条对角线。

### 6.3 式 (2.5) 在代码中如何成立

论文式 (2.5) 是

```text
rho*K*d_z(c/rho) = 0,  z=0 and z=z_H
```

也就是边界面通量

```text
J_(1/2) = 0
J_(N+1/2) = 0
```

代码用“边界面不进入矩阵”的方式实现：

1. `_lower`、`_diagonal`、`_upper` 初始全部为零；
2. 组装循环只遍历 `N-1` 个真实内部界面；
3. 因此第一层没有下边界面系数，最后一层没有上边界面系数；
4. `apply()` 在 `k=0` 不读取 `k-1`，在 `k=N-1` 不读取 `k+1`；
5. `solveImplicit()` 明确令首行 `lower=0`、末行 `upper=0`。

这与论文在式 (3.16) 后所说的“第一层令 `(rho K)^-=0`，最后一层令
`(rho K)^+=0`”代数完全相同。反射 ghost **位置**不表示令
`c_ghost=c_inside`；真正实施 Neumann/no-flux 条件的是外部面通量项为零。

边界两层实际使用：

```text
(A_z c)_1 =  J_(3/2) / Delta_z_1
(A_z c)_N = -J_(N-1/2) / Delta_z_N
```

没有 `J_(1/2)` 或 `J_(N+1/2)` 项。

### 6.4 这项边界条件如何验证

存在两个互补的不变量：

1. 若 `c_k = q*rho_k` 且 `q` 为常数，则所有内部面的 `q_(k+1)-q_k=0`，边界通量也为零，
   所以每层扩散倾向都应为零；
2. 对任意浓度，按层厚积分后内部面通量望远镜相消，边界通量为零，因此
   `sum_k Delta_z_k*(A_z c)_k=0`。

`MAS1998VerticalOperator.preservesConstantMixingRatioAndColumnMass` 同时检查这两点。
`implicitSolveIsNonnegative` 进一步检查无 loss 的隐式扩散保持非负，并保持层厚加权柱总量。

## 7. 垂直扩散系数与空气状态

`MAS1998::verticalDiffusivity()` 实现式 (4.2)，`z` 用 km、返回值用 `m^2 s^-1`：

```text
K = 30                              z <= 15
K = 0.2                       15 < z <= 17.5
K = 0.2 + 0.32*(z - 17.5)   17.5 < z <= 20
K = 10^(0.05*z - 1.0)              z > 20
```

垂直算子内部将中心距离和层厚转换为 m。`airNumberDensity()`、
`atmosphericTemperature()` 和 `atmosphericPressureMbar()` 使用 1976 US Standard Atmosphere
分层温度递减率近似，在 11、20、32 km 处分段。地面数密度固定为 `2.55e19 molecules cm^-3`。

这里复现的是论文所述标准大气关系的解析重建，不是已恢复的原始 CWI 大气查表数据，因此不应
期待旧程序逐位一致。

## 8. 化学机制与环境量

### 8.1 45 反应、17 个变量物种

机制源位于：

```text
test/tests/mas1998/chemistry/mas1998_methane.spc
test/tests/mas1998/chemistry/mas1998_methane.eqn
test/tests/mas1998/chemistry/mas1998_methane.def
```

`mas1998_methane.eqn` 保存恢复的最终 45 反应 benchmark 机制，包括 Arrhenius、三体反应和
光解速率。KPP 生成的变量顺序为：

```text
O1D CH4 HNO2 H2O2 N2O5 HNO3 HO2NO2 CH3OOH HCHO
CH3O2 NO3 O3P NO OH NO2 O3 HO2
```

solver 构造时要求共享库报告的名称与该顺序完全相等，否则立即报错，防止数组物种错位。

化学附录 NM-R9505 详细描述的是 46 反应、19 变量物种的前身，可用于核对速率公式和模型来源；
它不是本文最终全局 benchmark 的逐字机制定义。当前 45/17 KPP 文件仍是实际运行的权威输入。

### 8.2 Production/loss 接口

论文式 (2.2) 要求

```text
R_j(c) = P_j(c) - L_j(c)c_j
```

KPP 共享库导出的 `Fun_SPLIT` 直接返回 `P_j` 和标量 loss coefficient `L_j`。
MOOSE 通用类 `KPPGeneratedMechanism::computeProductionLoss()` 负责加载和调用该符号；app 的
`MAS1998BenchmarkSolver::evaluateRates()` 只负责提供 MAS1998 当地环境。

固定物种随高度按空气密度比例设置：

```text
M   = rho(z)
H2O = H2O_ground * rho(z)/rho(0)
CO  = CO_ground  * rho(z)/rho(0)
O2  = O2_ground  * rho(z)/rho(0)
```

此外每次速率计算都提供当地温度、mbar 压力、纬度、经度、高度和当地太阳时。必须先区分论文
绘图经度 `lambda'` 与地理经度 `lambda`：

```text
lambda = lambda' + 180 deg
local_solar_time = (GMT_seconds + lambda/360 deg*86400) mod 86400
```

`MAS1998::localSolarTime()` 统一实现这一转换，生产 solver 和特征线参考都调用它。若直接把
`lambda'` 当作地理经度，会把全球光化学昼夜相位错开 12 小时；单元测试以
`lambda'=-180 deg` 对应 Greenwich 午夜、`lambda'=0` 对应本地正午明确锁定该语义。

恢复机制中 `MAS1998_DELTA=0` 表示春分；太阳天顶角为夜间时，光解速率通过极大的 secant
衰减到接近零。

## 9. Implicit Euler、BDF2 与 Gauss-Seidel

实现类为 `MAS1998ChemistryIntegrator`，对应式 (3.20)-(3.24)。

### 9.1 变步长 BDF2 系数

令当前步长 `h_n`、前一步长 `h_(n-1)`，并令

```text
q = h_n / h_(n-1)
alpha = (1 + q)/(1 + 2q)
tau = alpha*h_n
C = [(1+q)^2*c_n - q^2*c_(n-1)]/(1+2q)
```

则一个物种的隐式关系为

```text
(I - tau*A_z + tau*diag(L_j)) c_j = C_j + tau*P_j
```

子问题第一步取 `alpha=1, C=c_n`，即 implicit Euler。后续使用上述 BDF2；若最后剩余时间
小于指定内步长，仍用不等步长系数完成末步。测试
`appliesEulerAndUnequalStepBDF2Coefficients` 对 `h={1,1,0.5}` 直接计算期望值。

### 9.2 三种退化形式

| 调用模式 | 方程与求解 |
| --- | --- |
| Type I 纯扩散 | `P=L=0`，每物种解 `(I-tau*A_z)c=C` |
| Type I 纯化学 | `A_z=0`，逐层逐物种做标量除法 `(C+tau*P)/(1+tau*L)` |
| Type II 化学+扩散 | 每物种一次解完整 15 层三对角系统 |

三对角系统由 `MAS1998VerticalOperator::solveImplicit()` 使用 Thomas 算法求解。矩阵
`I-tau*A_z+tau*diag(L)` 的边界行继承第 6 节的零外部面系数。

### 9.3 固定两次有序 sweep

每个 BDF 步的初猜对应式 (3.24)：

```text
first step: max(0, c_n)
BDF2 step : max(0, c_n + q*(c_n-c_(n-1)))
```

随后严格按 KPP 物种顺序更新。更新第 `j` 个物种时，速率计算能看到本 sweep 已更新的
`0..j-1` 和尚未更新的 `j..16`，因此是 Gauss-Seidel 而非 Jacobi。生产输入固定两次 sweep。

`updatesSpeciesInOrderForExactlyTwoSweeps` 用有顺序依赖的人工反应检查结果和速率计算次数；
`performsTwoSweepsAndRestoresNitrogen` 检查两次 sweep 后的氮修正。

## 10. 氮守恒修正与 NO 排放

17 个变量物种对应的氮原子数是：

```text
0 0 1 0 2 1 1 0 0 0 1 0 1 0 1 0 0
```

每个 sweep 后，所有含氮物种乘同一个比例

```text
correction = target_nitrogen / current_nitrogen
```

其中 N2O5 按两个氮原子计数。

- Type I chemistry-only：每个高度层单独计算比例，对应论文“each grid cell”；
- Type II chemistry-diffusion：先按 `Delta_z_k` 对整根柱积分，再用一个比例修正全柱。

论文第 4.1 节的最低层 NO 体积源为

```text
1e4 molecules/(cm^3 s)
```

实现将它加入最低层 NO 的 `production`，也把 `tau*emission` 加入 BDF 氮目标。后一项不可省略，
否则随后的氮修正会把外部排放误当作化学迭代误差消除。纯扩散子步不加入排放。

## 11. 初值、风场与论文 benchmark 参数

### 11.1 Table 1 与垂直初值

`MAS1998::groundConcentration()` 保存 Table 1 的地面背景值。除 HNO3 和 NO 外，所有物种水平
均匀。HNO3/NO 在以 `(lambda',phi')=(0,0)` 为中心的球面圆柱内分别取 `4e9` 和 `1e9`，
外部分别取 `2.55e9` 和 `1e2 molecules cm^-3`。

所有物种在垂直方向都保持初始 mixing ratio 不变：

```text
c(z) = c_ground * rho(z)/rho(0)
```

`standardAtmosphereAndInitialMixingRatio` 逐层检查 O3 的 `c/rho` 恒定，并检查 HNO3 圆柱内外值。

圆柱半径当前为显式输入 `plume_radius_degrees=20`。该值来自论文图示和恢复实现上下文；原始
CWI benchmark 数据缺失，仍存在来源不确定性，因此没有硬编码成不可修改的隐藏常数。

### 11.2 14 天固体旋转

`MAS1998::solidBodyWind()` 实现式 (4.1)：倾角 `beta=45 deg`，周期恰好 14 天，风场无散度、
与时间和高度无关。代码输入坐标为论文绘图坐标

```text
lambda' = (lambda-pi)*180/pi
phi' = phi*180/pi
```

并在计算风场前恢复 `lambda=lambda'+180 deg`。积分从 GMT 午夜 `t=0` 开始，在
`t=1,209,600 s` 结束。

## 12. 六个 Table 3 生产配置

公共配置在 `test/tests/mas1998/production/common.i`，六个输入只覆盖网格、分裂和步长：

| 网格 | `dt_adv` | `dt_split` | Type I `dt_dif` | 当前 chemistry/coupled step |
| --- | ---: | ---: | ---: | ---: |
| `64x32x15` | 2400 s | 4800 s | 1200 s | 300 s |
| `128x64x15` | 1200 s | 2400 s | 600 s | 300 s |
| `256x128x15` | 600 s | 1200 s | 300 s | 300 s |

论文 Table 3 的 `dt_che`/`dt_cvd` 是局部误差控制后的平均值：Type I 为
`313/312/300 s`，Type II 为 `301/322/311 s`；论文还说明最小步长是 `300 s`，且第 5.2 节
认为全程固定 5 min 也是合理选择。

论文没有公布足以重建控制器的全部容差和聚类细节。当前生产输入因此统一采用确定性的固定
`300 s`，不伪造缺失控制器。这是有论文依据的替代，但不是旧 CWI 自适应步序列的逐步复现。

完整 14 天时，每个 split interval 的内步数为：

| 网格 | Type I chemistry 步数 | Type II coupled 步数 | 每个 Type I 扩散半步 |
| --- | ---: | ---: | ---: |
| `64x32` | 16 | 16 | 2 |
| `128x64` | 8 | 8 | 2 |
| `256x128` | 4 | 4 | 2 |

## 13. 输出、论文曲线与参考解

每次输出产生两个 CSV：

1. 主 CSV：最低垂直中心 `z=0.3 km` 上、目标线 `lambda'=phi'/2` 的 17 个物种和
   `NOx=NO+NO2`；
2. `.summary.csv`：每个物种的全局空间加权积分和所有状态的全局最小值。

约化网格中心通常不精确落在 `lambda'=phi'/2` 上。`writeDiagnostics()` 选取包含目标经度的
真实单元，并同时写出该单元实际 `lambda_prime`；它不做经度插值。与论文曲线比较时应使用
CSV 中的实际坐标，而不要假设每一点严格位于解析对角线上。

summary 的积分权重为

```text
cell_measure * Delta_z
```

它适合检查相对守恒和不同时间的变化。代码没有把 `molecules cm^-3` 乘以 `1e6` 转成
`molecules m^-3`，所以 `_total` 是一致的加权诊断量，不应直接标作物理“分子总数”。

`scripts/plot_diagonal.py` 读取每个文件的最后时刻和最低层，绘制 O3、NOx、HNO3、HO2NO2，
对应论文 Figs. 4-8 的诊断物种。`scripts/compare_diagonal.py` 在坐标对齐通过后计算误差。

### 13.1 作者 CWI reference 是什么

论文 Section 4.5 的 reference 不是解析化学解，而是作者另外生成的高精度数值解。其关键构造
是水平固体旋转沿解析特征线精确处理，轨迹上只积分化学与论文同一套 15 层垂直半离散系统，
因此不含生产网格的水平空间离散误差。

历史页面 `http://www.cwi.nl/ftp/edwins/Ref_Sol_Benchmark_Global.html` 在 Internet Archive CDX
中有 `19970410111732`、`19980125163638`、`20000902233637`、`20010831050617` 四个时间戳，
但本次归档 replay 未能取回原始数值 payload。故当前不存在可验证字节来源和 SHA-256 的
`CWI_AUTHOR_REFERENCE`；不能把任何新计算重命名为原作者数据。

### 13.2 当前特征线参考如何生成

`MAS1998CharacteristicReferenceSolver` 对每个生产网格分别生成一个坐标完全对齐的参考。对每个
纬圈按以下步骤执行：

1. 使用与 `writeDiagnostics()` 相同的规则，取目标 `lambda'=phi'/2` 所落入真实单元的中心；
2. 从最终采样点用解析固体旋转 `-T` 反向追踪到 `t=start_time` 的特征线脚点；
3. 在脚点按 Table 1、HNO3/NO 圆柱和 `c(z)/rho(z)=constant` 初始化 `15 x 17` 状态；
4. 在每个 BDF 速率评估时，以解析旋转得到该绝对时刻的经纬度，并据此计算当地太阳时；
5. 一次性积分同一 `MAS1998ChemistryIntegrator`，设置 `chemistry=true,diffusion=true`，保留同一
   KPP、固定物种、NO 最低层排放、垂直矩阵、式 (2.5) 边界和氮修正；
6. 在 `t=end_time` 输出最低层，并写入证据标签 `INDEPENDENT_CHARACTERISTIC_REFERENCE`。

14 天恰好是一周固体旋转，故最终点与初始脚点重合；这不意味着中间化学路径可以省略，因为
太阳时、纬度和光解率沿轨迹变化。600 s smoke 故意使用非闭合时间，以实际覆盖反向脚点逻辑。

三个候选输入是：

```text
reference/characteristic_64x32.i
reference/characteristic_128x64.i
reference/characteristic_256x128.i
```

它们都用 `75 s` 内步和 8 次 Gauss-Seidel sweep，但在收敛门槛通过前只能称“候选参考”。

### 13.3 参考自身的收敛门槛

在 `64x32` 的同一组采样坐标上生成 `(150 s,4 sweeps)`、`(75 s,8 sweeps)` 和
`(37.5 s,12 sweeps)`。以最细者为 oracle，对 O3、NOx、HNO3、HO2NO2 每一个物种要求：

```text
relative_L2(75/8 vs 37.5/12)   <= 1.0e-3
relative_Linf(75/8 vs 37.5/12) <= 2.0e-3
```

而且 `75/8` 的这两个范数都必须严格小于 `150/4`。`compare_diagonal.py` 的
`--require-reference-convergence` 会自动执行该判据。若失败，不得放宽阈值；应把 `37.5/12`
提升为候选，再新增 `18.75/16` 作为 oracle 复查。

### 13.4 六个生产解如何定量比较

比较器只读取最终时刻、最低层，按 `row` 排序，并要求下列信息逐项一致：

```text
n_longitude, n_latitude, reduced_grid, final time,
row, cell, z_km, lambda_prime, phi_prime
```

经纬度容差是 `1e-12 deg`，不允许隐藏插值。对每个物种和每个 Type I/II case 输出：

```text
relative_L1   = sum(abs(case-ref)) / sum(abs(ref))
relative_L2   = sqrt(sum((case-ref)^2) / sum(ref^2))
relative_Linf = max(abs(case-ref)) / max(abs(ref))
absolute_Linf = max(abs(case-ref))
```

NaN、Inf、零参考分母、行数/网格/坐标/终止时间不一致都会使任务失败，而不是跳过坏点。
最终应得到 `errors_64x32.csv`、`errors_128x64.csv`、`errors_256x128.csv`。

### 13.5 旧的规则网格细化配方

`test/tests/mas1998/reference/type_ii_refined.i` 使用：

```text
regular 256x128x15 grid
dt_adv = 50 s
dt_split = 100 s
chemistry/coupled step = 20 s
6 Gauss-Seidel sweeps
```

这是一个较细的 Eulerian PDE 对照。它仍有水平网格、limiter 和分裂误差，既不能冒充原始
CWI reference，也不应在特征线参考收敛后作为首选 reference。

### 13.6 长算例运行与耗时

`scripts/run_quantitative_reproduction.sh` 默认用受控并发运行六个生产输入和三个对齐参考；加
`--with-convergence` 后再运行两个 `64x32` 收敛输入并强制执行第 13.3 节判据。每个任务各有
`.log` 和 `.time.txt`，终止后必须存在 `t=1209600`，才会生成误差 CSV 和叠加图。

2026-08-17 在本机单进程测得：`64x32` Type I 的一个 4800 s split 用 28.35 s，Type II 用
29.90 s，RSS 约 150 MB。据此线性外推每个生产 case：粗网格 2.0-2.1 h，中网格 7.5-7.9 h，
细网格 31-33 h；六个串行约 84 h。这些是外推而非完整计时，并发时会受 CPU 竞争影响。

后台完整命令为：

```bash
cd /home/fangxiaozhong/git_repo/mas1998_benchmark
mkdir -p runs
nohup scripts/run_quantitative_reproduction.sh \
  --jobs 4 --with-convergence \
  --output-dir runs/mas1998_quantitative_2026-08-17 \
  > runs/mas1998_quantitative_2026-08-17.launch.log 2>&1 &
```

自定义 executioner 没有中途 checkpoint/restart；中断的单个任务必须从 `t=0` 重算。

## 14. 验证证据能说明什么

最近一次完整开发验证记录如下：

| 验证层 | 结果 | 覆盖范围 |
| --- | ---: | --- |
| Atmospheric Chemistry 定向 KPP tests | 3/3 通过 | KPP 共享库通用能力 |
| app unit tests | 17/17 通过 | 网格、物理、坐标/轨迹、平流、垂直和化学测试，另含 2 个模板 sample test |
| app non-heavy tests | 3/3 通过 | 静态物理、既有 FV/基础 app 路径 |
| app heavy tests | 8/8 通过 | KPP 生成、0D box、Type I/II/参考 600 s 全链路及 CSV 检查 |
| Python 定量比较 tests | 4/4 通过 | 已知范数、元数据/坐标/非有限值拒绝和收敛判据 |
| 九个主输入与两个收敛输入 | 11/11 Syntax OK | 参数、对象构造和输入继承关系 |
| 批量运行器 dry-run | 默认 9、收敛模式 11 条命令 | 输入选择、并发参数及无输出目录副作用 |
| 14 天生产/参考运行 | 未执行 | 计算成本高，不在 CI 中；需按第 13.6 节后台执行 |

关键测试和它们真正证明的内容：

| 测试 | 证明 | 不证明 |
| --- | --- | --- |
| `reproducesTable2` | 三组总 cell/unknown 数与 Table 2 一致 | 细网格逐纬圈表来自原程序 |
| `conservesIntegratedMass` | 水平内部面通量离散守恒 | 14 天平流误差与论文一致 |
| `paperStepsMeetCFLBound` | 三个论文平流步满足实现的约化网格 CFL 上界 | 所有浓度场都绝不触发非线性 limiter 误差 |
| `preservesConstantMixingRatioAndColumnMass` | 垂直离散及式 (2.5) 边界满足两个核心不变量 | 垂直剖面与 CWI 逐位一致 |
| `appliesEulerAndUnequalStepBDF2Coefficients` | 式 (3.21) 系数和末短步正确 | 缺失的自适应控制器已复现 |
| `updatesSpeciesInOrderForExactlyTwoSweeps` | species 顺序和两次 Gauss-Seidel sweep 正确 | 两次 sweep 已收敛到完全隐式解 |
| Type I/II smoke + summary | 600 s 全链路可运行、输出有限且非负 | 14 天结果或论文 Figs. 4-8 已复现 |
| characteristic smoke + checker | 精确轨迹、化学-扩散柱和证据 CSV 在短时全链路可运行 | `75/8` 已达到 14 天收敛门槛 |
| synthetic comparison tests | 四种误差公式和拒绝条件实现正确 | 真实 14 天误差大小满足论文结果 |

## 15. 仍需保留的复现边界

1. 原始 CWI `Ref_Sol_Benchmark_Global.text` 不可用，当前独立生成结果不能验证为作者 gold 数组。
2. `128x64`、`256x128` 约化纬圈表是匹配 Table 2 的重建。
3. 生产化学步固定为 `300 s`，没有重建未公开的局部误差控制器和 C90 cluster 步序列。
4. 1976 标准大气使用解析分层重建，不保证与原 Fortran 查表逐位一致。
5. 20 度初始圆柱半径仍带有原始 benchmark 数据缺失造成的来源不确定性。
6. NM-R9505 是 46 反应/19 物种前身；运行机制是恢复的最终 45 反应/17 变量 KPP 版本。
7. 当前 solver 限制单 MPI rank，没有复现论文的 Cray C90 向量化、autotasking 和性能数据。
8. 尚未实际运行六个生产 case、三个参考候选和两个参考收敛 case；因此误差表仍待生成。

“定量复现完成”的最低闭环是：参考收敛门槛通过，九个主任务都存在合法最终行，三份误差表和
叠加图成功生成，并如实保留所有范数。作者 payload 仍缺失时，报告中必须继续使用
`INDEPENDENT_CHARACTERISTIC_REFERENCE`，而不是 `CWI_AUTHOR_REFERENCE`。

## 16. 代码索引

| 主题 | app 文件或 MOOSE 文件 |
| --- | --- |
| 顶层求解和分裂 | `src/utils/MAS1998BenchmarkSolver.C` |
| MOOSE 生产/参考输入入口 | `src/executioners/MAS1998Executioner.C`, `src/executioners/MAS1998ReferenceExecutioner.C` |
| 约化网格 | `src/utils/MAS1998ReducedGrid.C` |
| 水平平流、limiter、RK、CFL | `src/utils/MAS1998TransportOperator.C` |
| 15 层扩散、式 (2.5)、Thomas solve | `src/utils/MAS1998VerticalOperator.C` |
| Euler/BDF2、Gauss-Seidel、氮修正 | `src/utils/MAS1998ChemistryIntegrator.C` |
| Table 1、风场、K(z)、大气和初值 | `src/utils/MAS1998BenchmarkUtils.C` |
| 解析轨迹和坐标正确的当地太阳时 | `src/utils/MAS1998BenchmarkUtils.C` |
| Section 4.5 特征线参考生成 | `src/utils/MAS1998CharacteristicReferenceSolver.C` |
| 45 反应机制 | `test/tests/mas1998/chemistry/mas1998_methane.eqn` |
| 六个生产输入 | `test/tests/mas1998/production/` |
| 特征线候选、收敛及旧细化配方 | `test/tests/mas1998/reference/` |
| 定量误差、运行器与绘图 | `scripts/compare_diagonal.py`, `scripts/run_quantitative_reproduction.sh`, `scripts/plot_diagonal.py` |
| app 单元测试 | `unit/src/MAS1998*Test.C` |
| MOOSE 通用 KPP P/L API | `modules/atmospheric_chemistry/{include,src}/utils/KPPGeneratedMechanism.*` |

上表除最后一行外，路径均相对于
`/home/fangxiaozhong/git_repo/mas1998_benchmark`；最后一行相对于
`/home/fangxiaozhong/git_repo/moose`。
