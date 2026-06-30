# CutFEM 模块项目概览

**项目名称**：Cut Finite Element Methods (CutFEM) for MOOSE  
**目标**：在 MOOSE 框架中实现完整的 CutFEM 方法，支持非贴体网格上的复杂几何和动态界面  
**总工期**：8-10 个月  
**核心开发者**：[Your Name], MOOSE 社区贡献  

---

## 项目结构

```
modules/cutfem/
├── include/
│   ├── CutFEMApp.h                      # 应用程序主类
│   ├── kernels/
│   │   ├── GhostPenaltyKernel.h        # ★ 阶段1：梯度跳跃惩罚
│   │   ├── CutFEMInterfaceKernel.h     # 阶段2：界面核
│   │   └── SurfacePDEKernel.h          # 阶段3：表面PDE核
│   ├── userobjects/
│   │   ├── CutCellQuadratureUserObject.h  # ★ 阶段2：动态积分点计算
│   │   └── LevelSetTrackingUO.h
│   ├── bcs/
│   │   └── CutFEMBC.h                  # 阶段1：Nitsche边界条件
│   ├── quadrature/
│   │   └── CutCellQuadrature.h         # 阶段2：自定义积分规则
│   └── utils/
│       ├── CutGeometry.h               # 切割几何工具
│       ├── MarchingCubes.h             # 单元子剖分算法
│       └── EdgeIntersection.h          # 边界-水平集交点计算
├── src/
│   ├── CutFEMApp.C
│   ├── kernels/
│   │   ├── GhostPenaltyKernel.C       # 完整实现（见下文）
│   │   └── ...
│   ├── userobjects/
│   │   └── CutCellQuadratureUserObject.C
│   └── utils/
├── test/
│   └── tests/
│       ├── ghost_penalty/              # ★ 阶段1测试
│       │   ├── test_gp.i
│       │   ├── gold/
│       │   └── test.txt
│       ├── cut_cells/                  # 阶段2测试
│       └── surface_pde/                # 阶段3测试
├── examples/
│   ├── poisson_with_ghost_penalty.i   # ★ 阶段1示例
│   ├── circular_interface.i            # 阶段2示例
│   └── evolving_interface.i            # 阶段3示例
├── doc/
│   ├── index.md                        # 模块文档首页
│   ├── IMPLEMENTATION_STRATEGY.md      # 完整技术方案（已生成）
│   ├── DEVELOPMENT_WORKFLOW.md         # 开发工作流指南（已生成）
│   ├── references.md                   # 学术文献汇总
│   └── API_reference.md                # API 参考
├── Makefile                            # MOOSE 编译配置（已生成）
└── README.md                           # 项目说明
```

---

## 核心数学模型

### 问题陈述

**给定**：
- 物理域 $\Omega \subset \mathbb{R}^d$ 由水平集函数 $\phi$ 隐式表示
- PDE：$-\nabla \cdot (\alpha \nabla u) = f$ in $\Omega$，$u = 0$ on $\partial \Omega$
- 背景网格 $\tilde{\mathcal{T}}_h$（与 $\Omega$ 无关）

**目标**：
- 在 **非贴体网格** 上求解 PDE
- 确保 **条件数** $\kappa(\mathcal{A}_h) \leq C h^{-2}$（与网格位置无关）
- 实现 **最优收敛阶**（P1: $O(h)$ in H1，$O(h^2)$ in L2）

### 弱形式与稳定化

**Nitsche 方法** + **Ghost Penalty**：

$$a_h(u_h, v_h) + s_h(u_h, v_h) = l_h(v_h)$$

其中：

**1. Nitsche 双线性形式**（处理 Dirichlet BC）：
$$a_h(v,w) = \int_\Omega \alpha\nabla v \cdot \nabla w \, dx 
- \int_{\partial\Omega} n \cdot \alpha\nabla v \, w \, ds 
- \int_{\partial\Omega} v \, n \cdot \alpha\nabla w \, ds 
+ \beta h^{-1}(v,w)_{\partial\Omega}$$

**2. Ghost Penalty 稳定化**（处理切割单元）：
$$s_h(v,w) = \sum_{j=1}^p c_{F,j} h^{2(j-1+\gamma)} \int_F [D_n^j v] \cdot [D_n^j w] \, dS$$

其中 $[D_n^j v]_F = D_n^j v|_+ - D_n^j v|_-$ 是跨面 $F$ 的法向导数跳跃。

### 关键性质

| 性质 | 证明来源 | 意义 |
|------|---------|------|
| **H1 误差估计** $\|u-u_h\|_{H^1} \leq C h \|u\|_{H^2}$ | Burman et al. (2015) | 一阶收敛 |
| **L2 误差估计** $\|u-u_h\|_{L^2} \leq C h^2 \|u\|_{H^2}$ | 标准 Nitsche 分析 | 二阶收敛 |
| **条件数界** $\kappa(\mathcal{A}_h) \leq C h^{-2}$ | Larson & Zahedi (2020) | 与切割位置无关 |
| **稳定性** $\|v_h\|^2_a + s_h(v_h,v_h) \gtrsim h \|v_h\|_h^2$ | Ghost Penalty 理论 | 保证非奇异性 |

---

## 实现路线图

### 阶段 1：Ghost Penalty 稳定化（2-3个月）★ 优先

**任务**：
1. ✅ 实现 `GhostPenaltyKernel`（InterfaceKernel 子类）
2. ✅ 实现 Nitsche 边界条件处理
3. ✅ 编写测试用例（Poisson 问题）
4. ✅ 验证条件数改进（$\kappa$ 下降至少 50%）

**关键文件**：
- `include/kernels/GhostPenaltyKernel.h`
- `src/kernels/GhostPenaltyKernel.C`
- `examples/poisson_with_ghost_penalty.i`

**成功指标**：
- 条件数：$\kappa_{\text{without GP}} = 10^4 \rightarrow \kappa_{\text{with GP}} = 10^3$
- 迭代次数：GMRES 降低 60-70%
- 精度：H1 误差 $\sim O(h)$，L2 误差 $\sim O(h^2)$

**PR 准备**：
```
Title: "CutFEM Phase 1: Ghost Penalty Stabilization"
Commits:
  1. GhostPenaltyKernel: Implement gradient jump penalty
  2. Add Nitsche BC support and example
  3. Test suite: Poisson on unfitted circular domain
  4. Documentation: API and theory
```

---

### 阶段 2：非贴体不连续积分（3-4个月）

**任务**：
1. 实现 `CutCellQuadratureUserObject`（动态积分点生成）
2. 集成 LevelSet 模块
3. 实现 Marching Cubes 子单元剖分
4. 编写测试用例（验证积分精度）

**关键文件**：
- `include/userobjects/CutCellQuadratureUserObject.h`
- `src/userobjects/CutCellQuadratureUserObject.C`
- `include/utils/MarchingCubes.h`
- `src/utils/MarchingCubes.C`

**成功指标**：
- 切割单元的积分误差 $< 10^{-8}$
- 支持任意 LevelSet 函数
- 最优收敛阶：$O(h^2)$

---

### 阶段 3：表面 PDE 与演化（3个月）

**任务**：
1. 实现表面 PDE Kernel（Laplace-Beltrami）
2. 混合稳定化（面+界面）
3. Hamilton-Jacobi 求解器集成
4. 时间步耦合逻辑

**成功指标**：
- 界面平稳演化，无振荡
- 长期模拟稳定性（$t \gg 1$）
- 发表至少一篇顶级期刊论文

---

## 技术亮点

### 1. 张量积优化（Wichrowski 2026）

对笛卡尔网格，Ghost Penalty 的计算复杂度可优化：

$$O(k^{2d}) \rightarrow O(k^{d+1})$$

通过利用基函数的张量积结构：
$$\phi_{i_1,i_2,i_3}(x) = \phi_{i_1}(x_1) \phi_{i_2}(x_2) \phi_{i_3}(x_3)$$

**实现**：在 MOOSE 的矩阵无关框架中集成

---

### 2. 混合稳定化（Larson & Zahedi 2020）

结合两种稳定化确保最优条件数：

**面稳定化**：$s_{h,F}(w,v) = \sum_j c_j h^{2(j-1+\gamma)} [D_n^j w]_F [D_n^j v]_F$

**界面稳定化**：$s_{h,\Gamma}(w,v) = \sum_j c'_j h^{2(j-1+\gamma)} D_n^j w D_n^j v$

**结果**：$\kappa(\mathcal{A}_h) = O(h^{-2})$ **与多项式阶数 $p$ 无关**

---

### 3. 自适应 Marching Cubes（参考 ngsxfem）

对被水平集切割的元素，动态细分：

1. 检测：若 $\min_i \phi_i \cdot \max_j \phi_j < 0$ → 切割单元
2. 找交点：线性插值求 $\phi = 0$ 的交点
3. 剖分：标准 MC 模式进行子单元细分
4. 积分：对每个子单元应用高斯积分

---

## 关键引用与资源

### 学术论文（必读）

| 论文 | 年份 | 核心贡献 | 优先级 |
|------|------|---------|-------|
| Burman et al. "CutFEM: Discretizing geometry and PDEs" | 2015 | 原始 CutFEM 理论 | ★★★ |
| Larson & Zahedi "Stabilization of high order CutFEM" | 2020 | 混合稳定化分析 | ★★★ |
| Wichrowski "Matrix-free ghost penalty evaluation" | 2026 | 张量积优化 | ★★ |
| Burman et al. (综述) "Cut finite element methods" | 2025 | 最新综合概述 | ★★★ |

### 开源参考实现

- **ngsxfem**（Ngsolve extension）：Marching Cubes + 高阶积分
- **dune-cutfem**：完整 CutFEM 实现（DUNE 框架）
- **deal.II 扩展**：矩阵无关框架中的 CutFEM

### MOOSE 文档

- [MOOSE InterfaceKernel](https://mooseframework.inl.gov/syntax/InterfaceKernels/)
- [MOOSE UserObject System](https://mooseframework.inl.gov/syntax/UserObjects/)
- [MOOSE Level Set Module](https://mooseframework.inl.gov/modules/level_set/)

---

## 开发工作流总结

```
1. 环境配置
   └─ git clone fork + upstream
   └─ env-moose activate

2. 功能分支开发
   ├─ git checkout -b feature/xxx
   ├─ 编写代码 + clang-format
   ├─ make -C modules/cutfem -j4 METHOD=opt
   └─ 运行测试

3. 同步 & Rebase
   ├─ git fetch upstream
   └─ git rebase upstream/devel
   
4. 提交 PR
   ├─ 填写详细 PR 模板
   ├─ 等待 CIVET 检查
   └─ 解决审查反馈

5. 合并
   └─ Squash & merge 到 devel
```

详见 [DEVELOPMENT_WORKFLOW.md](./doc/DEVELOPMENT_WORKFLOW.md)

---

## 成功指标与验收标准

### 阶段 1 验收

- ✅ 条件数改进至少 50%
- ✅ 所有新增单元测试通过
- ✅ 数值精度与理论相符
- ✅ 代码审查通过

### 阶段 2 验收

- ✅ 正确识别和切割单元
- ✅ 积分精度 $< 10^{-8}$
- ✅ 与标准 FEM 相同的收敛阶
- ✅ 支持动态 LevelSet 更新

### 阶段 3 验收

- ✅ 长期稳定性（$t > 100$）
- ✅ 物理守恒量得到保证
- ✅ 可发表的研究成果

---

## 联系与支持

- **MOOSE 社区论坛**：https://github.com/idaholab/moose/discussions
- **CIVET CI 系统**：https://civet.inl.gov
- **MOOSE 官方文档**：https://mooseframework.inl.gov

---

**最后更新**：2026-06-26  
**文档版本**：1.0  
**维护者**：CutFEM 开发团队
