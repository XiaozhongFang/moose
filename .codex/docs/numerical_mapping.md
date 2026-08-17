# MAS1998 数值实现与论文映射

## 1. 总方程与分裂

论文以物种数浓度 `c` 求解

```text
dc/dt = F0(c) + F1(c) + F2(c)
F0(c) = R(c)
F1(c) = d_z[rho K d_z(c/rho)]
F2(c) = -1/(a cos(phi)) [d_lambda(u c) + d_phi(v c cos(phi))]
```

实现中的对应关系如下。

| 论文 | 实现 | 验证 |
| --- | --- | --- |
| 式 (2.5) 零垂直通量边界 | `MAS1998VerticalOperator` 的首末外部面系数为零 | `MAS1998VerticalOperator.preservesConstantMixingRatioAndColumnMass` |
| 式 (3.2a-e) Type I | `MAS1998BenchmarkSolver::advanceTypeI()`：`A(dt/2) D(dt/2) C(dt) D(dt/2) A(dt/2)` | `mas1998_complete_type_i_smoke` 与 summary 检查 |
| 式 (3.3a-c) Type II | `MAS1998BenchmarkSolver::advanceTypeII()`：`A(dt/2) [C+D](dt) A(dt/2)` | `mas1998_complete_type_ii_smoke` 与 summary 检查 |

`dt_split = 2*dt_adv` 由 solver 构造时强制检查。风场与高度无关，两个平流半步均只需一个
`dt_adv` RK 步。

## 2. 水平约化网格和平流

| 论文 | 数值含义 | 实现与测试 |
| --- | --- | --- |
| 式 (3.5)-(3.6) | 球面单元守恒散度与面通量 | `MAS1998TransportOperator::computeRHS()`；每面等量反号累计；`conservesIntegratedMass` |
| 式 (3.7)-(3.9) | 三阶迎风重构与 limiter | `limitedFaceValue()`；`limiterUsesPublishedFormula` |
| 式 (3.10)-(3.11) | 两阶段显式梯形 RK | `MAS1998TransportOperator::advance()`；接缝 pulse 的质量/非负测试 |
| 式 (3.12)-(3.13) | `max(nu_lambda+nu_phi) <= 2/3` | `maximumCFL()`；真实东西面和南北约化面段上界；`paperStepsMeetCFLBound` |
| Fig. 1 的约化接口 | piecewise-constant virtual concentrations | `MAS1998ReducedGrid::InterfaceSegment`；接口映射和周期缝测试 |

经向面长度为 `a*dphi`，纬向面段长度为 `a*cos(phi_face)*d_lambda_segment`，单元测度为
`a^2*cos(phi_center)*d_lambda*dphi`。公共面只计算一次，所以不依赖两个相邻纬圈拥有相同经度
单元数。经度索引周期化闭合 `-pi/pi` 接缝。

Table 2 的离散规模映射为：

| 名义网格 | 极点到赤道的每纬圈经度数 | 3D cells | 17 物种 unknowns |
| --- | --- | ---: | ---: |
| `64x32x15` | `4,8,16,32,12x64` | 24,840 | 422,280 |
| `128x64x15` | `4,8,16,32,8x64,20x128` | 93,960 | 1,597,320 |
| `256x128x15` | `4,8,16,32,64,17x128,42x256` | 391,560 | 6,656,520 |

各列表再关于赤道镜像。`MAS1998ReducedGrid.reproducesTable2` 对总数作精确断言。后两组纬圈
表属于满足 Table 2 的重建，证据限制见 `architecture.md`。

## 3. 垂直扩散

式 (3.15) 的 15 个中心为：

```text
0.3, 1.0, 2.2, 4.3, 6.5, 8.4, 10.0, 11.3,
13.0, 15.2, 17.6, 19.8, 22.5, 27.6, 34.7 km
```

底和顶边界为 `0` 与 `38.2 km`。`MAS1998VerticalOperator` 用关于边界反射的 ghost center
得到控制层厚；式 (3.16) 中的 `rho*K` 在相邻中心中点求值。组装矩阵 `A` 满足式
(3.17)-(3.19)，并通过

```text
(I - tau*A + tau*diag(loss)) c = rhs
```

的 Thomas 求解用于化学-扩散和纯扩散子问题。单元测试覆盖恒定 mixing ratio 零倾向、层厚
加权守恒、隐式解非负与三对角求解守恒。

## 4. 化学、Euler/BDF2 与 Gauss-Seidel

KPP `Fun_SPLIT` 给出

```text
R_j(c) = P_j(c) - L_j(c) c_j
```

并由通用模块 `KPPGeneratedMechanism::computeProductionLoss()` 暴露。物种顺序严格使用生成机制：

```text
O1D CH4 HNO2 H2O2 N2O5 HNO3 HO2NO2 CH3OOH HCHO
CH3O2 NO3 O3P NO OH NO2 O3 HO2
```

| 论文 | 实现 | 验证 |
| --- | --- | --- |
| 式 (3.20) | `MAS1998ChemistryIntegrator::takeStep()` 的隐式 species solve | KPP 完整 smoke |
| 式 (3.21) | 子区间首步 `alpha=1,C=c_n`；其后变步长 BDF2 | `appliesEulerAndUnequalStepBDF2Coefficients` |
| 式 (3.22)-(3.23) | `tau=alpha*h`，逐物种标量 loss，逐柱三对角解 | chemistry 与 vertical unit tests |
| 式 (3.24) | 非负外推初猜，KPP 顺序更新，固定两次 sweep | `updatesSpeciesInOrderForExactlyTwoSweeps` 与 evaluation count |
| 第 3.6 节 | 每次 sweep 后修正氮守恒 | `performsTwoSweepsAndRestoresNitrogen` |

Type I 的 chemistry-only 解令 `A=0`，每个高度层独立修正氮；Type II 在整根柱上隐式耦合
扩散，并以层厚加权氮总量统一修正。最低层 NO 的 `1e4 molecules/(cm^3 s)` 排放既加入
production，也加入 BDF 氮守恒目标，避免修正步骤错误消除外部源。

论文使用未完全公开参数的局部误差控制器，且以 300 s 为最小步。生产输入采用恒定
300 s；积分器仍实现不等末步所需的式 (3.21)，但不声称还原缺失的控制器。

## 5. Benchmark 数据：式 (4.1)-(4.2) 与 Table 1

| 来源 | 实现 |
| --- | --- |
| 式 (4.1) 14 天、倾角 45 度固体旋转 | `MAS1998::solidBodyWind()` |
| 式 (4.2) 分段 `K(z)` | `MAS1998::verticalDiffusivity()` |
| 1976 标准大气近似 | `atmosphericTemperature()`、`atmosphericPressureMbar()`、`airNumberDensity()` |
| Table 1 地面浓度 | `groundConcentration()` 与 `MAS1998BenchmarkUtils.publishedMetadata` |
| HNO3/NO 柱状初值 | `initialConcentration()` 与 `MAS1998SpeciesIC` |
| 垂直恒 mixing ratio | `initialConcentration()` 乘 `rho(z)/rho(0)`；unit test 逐层检查 |

17 个变量物种的氮原子数为：

```text
0 0 1 0 2 1 1 0 0 0 1 0 1 0 1 0 0
```

固定物种 `M,H2O,CO,O2` 随空气密度比例按层设置。app 同时给 KPP 环境提供当地温度、mbar
压力、纬度、经度和当地太阳时；`MAS1998_DELTA=0` 对应基准的春分设置。

## 6. Table 3 生产配置

六个输入位于 `test/tests/mas1998/production/`：

| 网格 | `dt_adv` | `dt_split` | Type I `dt_dif` | chemistry/coupled step | 输入 |
| --- | ---: | ---: | ---: | ---: | --- |
| `64x32x15` | 2400 s | 4800 s | 1200 s | 300 s | `type_{i,ii}_64x32.i` |
| `128x64x15` | 1200 s | 2400 s | 600 s | 300 s | `type_{i,ii}_128x64.i` |
| `256x128x15` | 600 s | 1200 s | 300 s | 300 s | `type_{i,ii}_256x128.i` |

全部从 GMT 午夜运行恰好 14 天并使用两次 Gauss-Seidel sweep。Table 3 中约 300 s 的数值是
自适应步长平均值而非三个算例都严格等于 300 s；此处的恒定步选择是论文第 5.2 节明确认为
可行的确定性替代，并在 README 中公开。

## 7. 输出与参考解状态

`MAS1998BenchmarkSolver` 输出最低层 `lambda'=phi'/2` 剖面，包含 17 个物种与
`NOx=NO+NO2`；summary 输出各物种全局体积加权总量和全局最小值。`plot_diagonal.py` 可叠加
Type I、Type II 与 independent refined 的 O3、NOx、HNO3、HO2NO2，分别对应论文 Fig. 4-8
的诊断量。

`reference/type_ii_refined.i` 使用规则 `256x128x15`、`dt_adv=50 s`、`dt_split=100 s`、
`chemistry_step=20 s` 和六次 sweep。它只定义独立现代参考配方；原始
`Ref_Sol_Benchmark_Global.text` 缺失，因此没有论文 gold 数组可作数值误差断言。
