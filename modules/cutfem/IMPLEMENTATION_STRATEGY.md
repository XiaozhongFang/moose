# CutFEM in MOOSE: 完整实现方案

## 执行摘要

本文档为在 MOOSE 框架中实现完整的 CutFEM（切割有限元方法）提供了详细的技术路线。根据数学复杂度和工程难度，采用**三层递进式开发**策略，避免从零重新写有限元引擎的陷阱，而是通过"高价值嫁接"充分利用 MOOSE 现有机制。

---

## 第一阶段：鬼惩罚项实现（Ghost Penalty）

### 时间估计：2-3 个月
### 难度等级：★★☆☆☆ 中等

### 1.1 数学基础

**鬼惩罚项（Ghost Penalty）定义：**

$$G(u,v) = \sum_{F \in \mathcal{F}_G} \gamma h_F^{2k-1} \int_F [\nabla^k u] \cdot [\nabla^k v] \, dS$$

其中：
- $[\cdot]$ 表示梯度跳跃（gradient jump）
- $F$ 是内部边界（internal faces）
- $\gamma$ 是惩罚参数（通常取 $\gamma = \mathcal{O}(1)$）
- $h_F$ 是边界特征尺度
- $k$ 是导数阶数（CutFEM 中通常 $k=1$）

### 1.2 MOOSE 切入点与 Nitsche 方法

**核心组件：`InterfaceKernel` 或 `InternalSideIntegralKernel`**

在 MOOSE 中，Ghost Penalty 本质上是一个作用在**内部边界**上的边界积分，与 DG（不连续 Galerkin）方法的惩罚项类似。完整的 Nitsche 方法弱形式为：

$$a_h(u_h, v) = \int_\Omega \alpha \nabla u_h \cdot \nabla v \, dx - \int_{\partial \Omega} n \cdot \alpha \nabla u_h \cdot v \, ds$$
$$- \int_{\partial \Omega} u_h \cdot n \cdot \alpha \nabla v \, ds + \beta h^{-1}(u_h, v)_{\partial \Omega}$$

其中惩罚参数 $\beta = \beta_0 \alpha$，$\beta_0 = 30k(k+1)$（$k$ 为多项式阶数）。

#### 1.2.1 关键 API 与数学对应

| MOOSE 组件 | 数学含义 | 公式对应 |
|-----------|---------|--------|
| `InterfaceKernel` | 内部面积分 | $\int_F [\nabla^k u] \cdot [\nabla^k v]$ |
| `_grad_u_current` | $\nabla u_+$ | 当前单元梯度 |
| `_grad_u_neighbor` | $\nabla u_-$ | 邻近单元梯度 |
| `_normal` | $\mathbf{n}_F$ | 面法向（从当前到邻近） |
| `_h` | $h_F$ | 面特征尺度 |
| `_JxW_face` | 积分权重 | 面积加权 |

#### 1.2.2 Ghost Penalty 的精确数学定义

根据 Larson & Zahedi (2020) 和 Wichrowski (2026) 的最新工作，CutFEM 的 Ghost Penalty 定义为：

**面上的梯度跳跃项：**
$$s_{h,F}(w,v) = \sum_{j=1}^{p} c_{F,j} \, h_F^{2(j-1+\gamma)} \int_F [D_n^j w] \cdot [D_n^j v] \, dS$$

**界面上的法向导数项（表面 PDE 的情况）：**
$$s_{h,\Gamma}(w,v) = \sum_{j=1}^{p} c_{\Gamma,j} \, h^{2(j-1+\gamma)} \int_{\Gamma_h} D_n^j w \cdot D_n^j v \, dS$$

其中：
- $D_n^j v = (\mathbf{n}_h \cdot \nabla)^j v$ 是 $j$ 阶法向导数
- $[D_n^j v]_F = D_n^j v|_+ - D_n^j v|_-$ 是跨面 $F$ 的梯度跳跃
- $h_F = \sqrt{|F|}$ 是面的特征尺度
- $\gamma \in [0,1]$ 是稳定化参数
- $c_{F,j}$ 和 $c_{\Gamma,j}$ 是 $\mathcal{O}(1)$ 的稳定化常数
- $p$ 是有限元多项式阶数

**条件数上界（证明见 Larson & Zahedi）：**
$$\kappa(\mathcal{A}_h) \leq C h^{-2}$$

其中条件数与网格参数 $h$ 的关系为 $O(h^{-2})$，与标准 FEM 相同。

#### 1.2.3 实现步骤

**步骤1：自定义 `InterfaceKernel`**

```cpp
// include/kernels/GhostPenaltyKernel.h
#pragma once

#include "InterfaceKernel.h"

class GhostPenaltyKernel : public InterfaceKernel
{
public:
  static InputParameters validParams();
  GhostPenaltyKernel(const InputParameters & params);

protected:
  virtual Real computeQpResidual(Moose::DG_RESIDUAL_TYPE type) override;
  virtual Real computeQpJacobian(Moose::DG_JACOBIAN_TYPE type) override;

private:
  // 惩罚参数
  const Real & _gamma;
  
  // 导数阶数（通常为1）
  const unsigned int & _k;
  
  // 梯度跳跃计算函数
  RealVectorValue computeGradientJump() const;
};
```

**步骤2：C++ 实现（基于论文算法）**

```cpp
// src/kernels/GhostPenaltyKernel.C
#include "GhostPenaltyKernel.h"

InputParameters
GhostPenaltyKernel::validParams()
{
  InputParameters params = InterfaceKernel::validParams();
  params.addParam<Real>("gamma", 1.0, "Penalty parameter γ ∈ [0,1]");
  params.addParam<unsigned int>("k", 1, "Order of derivative j (typically 1)");
  params.addParam<Real>("c_F", 1.0, "Face stabilization constant c_{F,j}");
  params.addParam<Real>("c_Gamma", 1.0, "Surface stabilization constant c_{Γ,j}");
  return params;
}

GhostPenaltyKernel::GhostPenaltyKernel(const InputParameters & params)
  : InterfaceKernel(params),
    _gamma(getParam<Real>("gamma")),
    _k(getParam<unsigned int>("k")),
    _c_F(getParam<Real>("c_F")),
    _c_Gamma(getParam<Real>("c_Gamma"))
{
  // 参数检查
  if (_gamma < 0.0 || _gamma > 1.0)
    mooseError("Penalty parameter gamma must be in [0,1]");
}

Real
GhostPenaltyKernel::computeQpResidual(Moose::DG_RESIDUAL_TYPE type)
{
  // 计算梯度跳跃：[∇u] = ∇u_+ - ∇u_-
  // 在 MOOSE 中：_grad_u_current = ∇u_+，_grad_u_neighbor = ∇u_-
  RealVectorValue grad_jump = _grad_u_current[_qp] - _grad_u_neighbor[_qp];
  
  // 计算法向导数的跳跃：[D_n^k u]_F = [∇u] · n
  // 其中 n = _normal[_qp]
  Real grad_jump_normal = grad_jump * _normal[_qp];
  
  // 计算面特征尺度：h_F = sqrt(|F|)
  // _JxW_face[_qp] 已经包含了 |F| 权重
  // 对于标准高斯积分，|F| ≈ sqrt(_JxW_face[_qp])
  Real h_F_squared = _JxW_face[_qp];
  
  // 计算稳定化系数：c_F * h_F^(2(k-1+γ))
  // 对于 k=1, γ=1：h_F^(2(0+1)) = h_F^2
  Real exponent = 2.0 * (_k - 1.0 + _gamma);
  Real penalty_coeff = _c_F * std::pow(h_F_squared, exponent / 2.0);
  
  // 对试函数计算法向导数：D_n^k(phi_j) = grad_phi_j · n
  Real test_grad_normal = _grad_phi[_j][_qp] * _normal[_qp];
  
  if (type == Moose::Element)
  {
    // 当前单元的残差贡献：
    // r^+ = penalty_coeff * [D_n^k u] * (D_n^k phi_i) 
    return penalty_coeff * grad_jump_normal * test_grad_normal;
  }
  else // type == Moose::Neighbor
  {
    // 邻近单元的残差贡献（相反符号）：
    // r^- = -penalty_coeff * [D_n^k u] * (D_n^k phi_i)
    Real neighbor_test_grad_normal = _grad_phi_neighbor[_j][_qp] * _normal[_qp];
    return -penalty_coeff * grad_jump_normal * neighbor_test_grad_normal;
  }
}

Real
GhostPenaltyKernel::computeQpJacobian(Moose::DG_JACOBIAN_TYPE type)
{
  // Jacobian 矩阵计算：∂r/∂u_h
  
  Real h_F_squared = _JxW_face[_qp];
  Real exponent = 2.0 * (_k - 1.0 + _gamma);
  Real penalty_coeff = _c_F * std::pow(h_F_squared, exponent / 2.0);
  
  Real test_grad_normal = _grad_phi[_j][_qp] * _normal[_qp];
  Real basis_grad_normal = _grad_phi[_i][_qp] * _normal[_qp];
  Real neighbor_test_grad_normal = _grad_phi_neighbor[_j][_qp] * _normal[_qp];
  Real neighbor_basis_grad_normal = _grad_phi_neighbor[_i][_qp] * _normal[_qp];
  
  switch (type)
  {
    case Moose::ElementElement:
      // ∂r^+_current / ∂u_+
      return penalty_coeff * basis_grad_normal * test_grad_normal;
      
    case Moose::ElementNeighbor:
      // ∂r^+_current / ∂u_-
      return -penalty_coeff * neighbor_basis_grad_normal * test_grad_normal;
      
    case Moose::NeighborElement:
      // ∂r^-_neighbor / ∂u_+
      return -penalty_coeff * basis_grad_normal * neighbor_test_grad_normal;
      
    case Moose::NeighborNeighbor:
      // ∂r^-_neighbor / ∂u_-
      return penalty_coeff * neighbor_basis_grad_normal * neighbor_test_grad_normal;
      
    default:
      mooseError("Unsupported Jacobian type");
  }
}
```

**稳定化参数选择指南**（基于 Wichrowski 2026）：

| 参数 | 推荐值 | 含义 |
|-----|-------|------|
| `gamma` | 1.0 | 最强稳定化（$\gamma=1$）；值越小稳定化越强 |
| `c_F` | 1.0 | 面惩罚系数，通常取 $\mathcal{O}(1)$ |
| `c_Gamma` | 1.0 | 界面惩罚系数，通常取 $\mathcal{O}(1)$ |
| `k` | 1 | 对于二阶 PDE，使用 $k=1$；对于高阶问题 $k \geq 2$ |



### 1.3 表面 PDE 的稳定化（Trace FEM）

对于嵌入在 $\mathbb{R}^d$ 中的 $n$ 维光滑流形 $\Gamma$ 上的椭圆型 PDE，Laplace-Beltrami 问题的稳定化形式为：

$$a_h(u_h, v_h) + s_h(u_h, v_h) = l_h(v_h)$$

其中：

**原始双线性形式（仅切向梯度）：**
$$a_h(w,v) = \int_{\Gamma_h} \nabla_{\Gamma_h} w \cdot \nabla_{\Gamma_h} v \, dS$$

**混合稳定化形式（Larson & Zahedi, 2020）：**
$$s_h(w,v) = s_{h,F}(w,v) + s_{h,\Gamma}(w,v)$$

- **边界面稳定化**：
$$s_{h,F}(w,v) = \sum_{j=1}^{p} c_{F,j} h^{2(j-1+\gamma)} \int_F [D_{n_F}^j w] \cdot [D_{n_F}^j v] \, dS$$

- **界面稳定化**（直接作用在 $\Gamma_h$ 上）：
$$s_{h,\Gamma}(w,v) = \sum_{j=1}^{p} c_{\Gamma,j} h^{2(j-1+\gamma)} \int_{\Gamma_h} D_{n_h}^j w \cdot D_{n_h}^j v \, dS$$

这种组合稳定化确保：
$$\|v_h\|_{L^2(\Gamma_h)}^2 + h^{2\gamma}\|\nabla v_h\|_{L^2(\Gamma_h)}^2 \lesssim h(a_h(v_h, v_h) + s_h(v_h, v_h))$$

从而得到条件数 $\kappa(\mathcal{A}_h) \leq C h^{-2}$（与多项式阶数 $p$ 无关）。

### 1.4 矩阵无关实现的张量积优化

根据 Wichrowski (2026)，对于笛卡尔网格上的 Ghost Penalty，可以利用张量积结构进行高效计算：

**基函数张量表示**（在 2D/3D 直角坐标系中）：
$$\phi_{i_1, i_2, i_3}(x_1, x_2, x_3) = \phi_{i_1}(x_1) \cdot \phi_{i_2}(x_2) \cdot \phi_{i_3}(x_3)$$

**梯度的张量分解**：
$$\partial_x \phi_{i_1, i_2, i_3} = \phi'_{i_1}(x_1) \cdot \phi_{i_2}(x_2) \cdot \phi_{i_3}(x_3)$$

因此 Ghost Penalty 的计算复杂度可以从 $O(k^{2d})$ 优化到 $O(k^{d+1})$（$k$ 为多项式阶数，$d$ 为维数）。

MOOSE 中可以通过 deal.II 的矩阵无关框架（matrix-free framework）来实现这一优化。

### 1.3 验证策略

#### 1.3.1 理论验证

**目标：验证条件数从 $\mathcal{O}(h^{-2})$ 下降到 $\mathcal{O}(h^{-2})$（稳定）**

| 指标 | 不加GP | 加GP | 预期改进 |
|-----|-------|------|--------|
| 条件数 | $10^4$ | $10^3$ | 1个数量级下降 |
| 迭代次数 | 150-200 | 50-80 | 减少 60-70% |
| 求解时间 | 2.5s | 0.8s | 加速 3 倍 |

#### 1.3.2 测试用例

**测试1：Poisson 方程 + 内部边界约束**

```
Domain: [0,1]^2
PDE: -Δu = 1, u|_∂Ω = 0
约束: u(x,y) 在直线 y=0.5 两侧值有跳跃（模拟界面）
```

**测试2：条件数监控**

使用 PETSc 的 `-ksp_monitor_true_residual` 和 `-ksp_compute_eigenvalues` 来实时监控：

```bash
mpirun -n 4 moose_app-opt -i input.i \
  -ksp_compute_eigenvalues \
  -ksp_type gmres \
  -pc_type ilu
```

### 1.4 MOOSE 输入文件示例（MOOSE input syntax）

```ini
[Mesh]
  type = GeneratedMesh
  dim = 2
  nx = 20
  ny = 20
  xmin = 0
  xmax = 1
  ymin = 0
  ymax = 1
[]

[Variables]
  [u]
    order = FIRST
    family = LAGRANGE
  []
[]

[Kernels]
  [diffusion]
    type = Diffusion
    variable = u
  []
  [source]
    type = BodyForce
    variable = u
  []
[]

[InterfaceKernels]
  [ghost_penalty]
    type = GhostPenaltyKernel
    variable = u
    neighbor_var = u
    boundary = internal_surface  # 指定内部边界
    gamma = 1.0
    k = 1
  []
[]

[BCs]
  [left]
    type = DirichletBC
    variable = u
    boundary = 'left right bottom'
    value = 0
  []
  [top]
    type = NeumannBC
    variable = u
    boundary = top
    value = 0
  []
[]

[Executioner]
  type = Steady
  solve_type = NEWTON
  petsc_options_iname = '-pc_type -pc_sub_type -ksp_type'
  petsc_options_value = 'asm ilu gmres'
  nl_rel_tol = 1e-8
[]

[Outputs]
  exodus = true
  csv = true
[]
```

---

## 第二阶段：非贴体不连续积分实现

### 时间估计：3-4 个月
### 难度等级：★★★★☆ 困难

### 2.1 核心挑战

CutFEM 的关键在于**单元被几何界面切割**，需要在运行时（runtime）识别并正确积分被切割的区域 $T \cap \Gamma_h$。

MOOSE 的默认积分使用固定的高斯积分点在整个单元上进行，无法处理这种非标准的积分域。

### 2.2 解决方案架构

#### 2.2.1 关键数据结构

**LevelSet 表示：**
$$\phi(\mathbf{x}) = \text{signed distance function}$$

- $\phi > 0$：域 $\Omega_+$
- $\phi < 0$：域 $\Omega_-$
- $\phi = 0$：界面 $\Gamma$

**切割单元检测：**
```cpp
bool isCutElement(const Elem * elem, const std::vector<Real> & phi_values)
{
  Real phi_min = *std::min_element(phi_values.begin(), phi_values.end());
  Real phi_max = *std::max_element(phi_values.begin(), phi_values.end());
  return (phi_min * phi_max) < 0; // 符号相反 = 单元被切割
}
```

#### 2.2.2 子单元切割算法

**参考实现：ngsxfem 的 Cell Cutting 算法**

对于被 LevelSet 切割的四面体/六面体，需要：

1. **找交点**（Edge-Level Set 零等位面交点）
2. **剖分子单元**（使用标准单纯剖分）
3. **计算新的积分点和权重**

**伪代码：**

```cpp
struct CutCellQuadrature
{
  std::vector<Point> quad_points;  // 积分点
  std::vector<Real> quad_weights;  // 积分权重
  
  void computeForCutCell(const Elem * elem, const LevelSetAux * phi_aux)
  {
    // 第1步：提取元素节点处的 Level Set 值
    std::vector<Real> phi_at_nodes = evaluateLevelSetAtNodes(elem, phi_aux);
    
    // 第2步：检测是否为切割单元
    if (!isCutElement(elem, phi_at_nodes)) {
      // 使用标准积分
      return;
    }
    
    // 第3步：找交点
    std::vector<Point> intersect_points = findEdgeIntersections(elem, phi_at_nodes);
    
    // 第4步：子单元剖分（例如使用 Marching Cubes）
    std::vector<std::vector<Point>> sub_elements = marchingCubes(
      elem->nodes(), intersect_points, phi_at_nodes);
    
    // 第5步：对每个子单元应用标准高斯积分
    for (const auto & sub_elem : sub_elements) {
      auto [sub_quad_pts, sub_quad_wts] = gaussQuadrature(sub_elem, order);
      quad_points.insert(quad_points.end(), sub_quad_pts.begin(), sub_quad_pts.end());
      quad_weights.insert(quad_weights.end(), sub_quad_wts.begin(), sub_quad_wts.end());
    }
  }
};
```

#### 2.2.3 MOOSE 集成策略

**方案A：自定义 UserObject**

```cpp
// include/userobjects/CutCellQuadratureUserObject.h
#pragma once

#include "ElementUserObject.h"

class CutCellQuadratureUserObject : public ElementUserObject
{
public:
  static InputParameters validParams();
  CutCellQuadratureUserObject(const InputParameters & params);

  virtual void initialize() override;
  virtual void execute() override;
  virtual void finalize() override;
  virtual void threadJoin(const UserObject & y) override;

  // 查询接口
  const std::vector<Point> & getQuadraturePoints(const Elem * elem) const;
  const std::vector<Real> & getQuadratureWeights(const Elem * elem) const;
  bool isCutElement(const Elem * elem) const;

private:
  // 存储每个单元的积分点和权重
  std::map<const Elem *, std::vector<Point>> _quad_points_map;
  std::map<const Elem *, std::vector<Real>> _quad_weights_map;
  std::set<const Elem *> _cut_elements;

  // LevelSet 变量名
  const VariableName & _level_set_var_name;

  // 算法参数
  const unsigned int & _quad_order;
  const Real & _level_set_tolerance;
};
```

**方案B：继承 Quadrature 基类（更深度集成）**

```cpp
// include/quadrature/CutCellQuadrature.h
#pragma once

#include "QBase.h"

class CutCellQuadrature : public QBase
{
public:
  CutCellQuadrature(unsigned int dim, unsigned int order);

  virtual void init_2D(const ElemType type, unsigned int p_order) override;
  virtual void init_3D(const ElemType type, unsigned int p_order) override;

private:
  // 子单元切割实现
  void cutElement(const Elem * elem, const LevelSetFunction & level_set);
};
```

### 2.3 LevelSet 与切割集成

#### 2.3.1 LevelSet Initialization

```ini
[AuxVariables]
  [phi]
    order = FIRST
    family = LAGRANGE
  []
[]

[AuxKernels]
  [compute_phi]
    type = FunctionAux
    variable = phi
    function = level_set_function
    execute_on = 'INITIAL TIMESTEP_BEGIN'
  []
[]

[Functions]
  [level_set_function]
    type = ParsedFunction
    # 示例：圆形界面 r=0.3，中心在 (0.5, 0.5)
    value = 'sqrt((x-0.5)^2 + (y-0.5)^2) - 0.3'
  []
[]
```

#### 2.3.2 CutCellQuadrature UserObject

```ini
[UserObjects]
  [cut_cell_quadrature]
    type = CutCellQuadratureUserObject
    level_set_variable = phi
    quadrature_order = 4
    level_set_tolerance = 1e-12
    execute_on = 'INITIAL TIMESTEP_BEGIN'
  []
[]
```

### 2.4 自定义 Kernel：表面 PDE 项

```cpp
// include/kernels/CutFEMKernel.h
#pragma once

#include "Kernel.h"

class CutFEMKernel : public Kernel
{
public:
  static InputParameters validParams();
  CutFEMKernel(const InputParameters & params);

protected:
  virtual Real computeQpResidual() override;
  virtual Real computeQpJacobian() override;

private:
  // 指向 CutCellQuadratureUserObject 的指针
  const CutCellQuadratureUserObject * _cut_cell_quad_uo;
  
  // 法向量计算
  RealVectorValue computeNormalVector() const;
};
```

---

## 第三阶段：演化流形耦合

### 时间估计：3 个月
### 难度等级：★★★★★ 极困难

### 3.1 动态 LevelSet 演化

#### 3.1.1 Hamilton-Jacobi 方程

$$\frac{\partial \phi}{\partial t} + V |\nabla \phi| = 0$$

其中 $V$ 是界面速度（通常从表面 PDE 计算得出）。

#### 3.1.2 MOOSE 时间步集成

```ini
[Executioner]
  type = Transient
  dt = 1e-3
  num_steps = 100
  
  # 在每个时间步执行的顺序
  # 1. 求解表面 PDE（获得界面速度）
  # 2. 更新 LevelSet
  # 3. 重新计算切割单元与积分点
  # 4. 下一个时间步
[]

[AuxKernels]
  [update_level_set]
    type = TimeDerivativeAux
    variable = phi_dot
    v = phi
    execute_on = 'TIMESTEP_END'
  []
  
  [hj_equation]
    type = CoupledForceAux
    variable = phi_dot
    v = velocity  # 从表面 PDE 计算的速度
    coefficient = 1.0
    execute_on = 'TIMESTEP_END'
  []
[]
```

### 3.2 耦合表面 PDE 求解器

```cpp
// include/kernels/SurfacePDEKernel.h
// 计算表面 PDE 的弱形式
// 涉及拉普拉斯-贝特拉米算子等
```

---

## 第四部分：完整工程实现清单

### 4.1 文件结构

```
modules/cutfem/
├── include/
│   ├── kernels/
│   │   ├── GhostPenaltyKernel.h           # 阶段1
│   │   ├── CutFEMKernel.h                 # 阶段2
│   │   ├── CutFEMInterfaceKernel.h        # 阶段2
│   │   └── SurfacePDEKernel.h             # 阶段3
│   ├── bcs/
│   │   └── CutFEMBC.h
│   ├── userobjects/
│   │   ├── CutCellQuadratureUserObject.h  # 阶段2 核心
│   │   └── LevelSetTrackingUO.h
│   ├── quadrature/
│   │   └── CutCellQuadrature.h            # 阶段2 深度集成
│   └── utils/
│       ├── CutGeometry.h                  # 切割算法
│       ├── MarchingCubes.h
│       └── EdgeIntersection.h
├── src/
│   ├── kernels/                          # 对应 include/
│   ├── bcs/
│   ├── userobjects/
│   ├── quadrature/
│   └── utils/
├── test/                                 # 单元测试
│   ├── tests/ghost_penalty/
│   ├── tests/cut_elements/
│   └── tests/surface_pde/
├── examples/                             # 示例用例
│   ├── poisson_with_gp.i
│   ├── circular_interface.i
│   └── evolving_interface.i
└── doc/
    └── index.md
```

### 4.2 关键开发清单

#### 阶段1 清单

- [ ] `GhostPenaltyKernel` 实现与测试
- [ ] 条件数验证测试
- [ ] 与 XFEM 模块兼容性检查
- [ ] 文档编写

#### 阶段2 清单

- [ ] LevelSet 数据结构设计
- [ ] `CutCellQuadratureUserObject` 实现
- [ ] 2D 切割算法实现与验证
- [ ] 3D 切割算法实现与验证
- [ ] 集成测试：单元被切割的 Poisson 问题

#### 阶段3 清单

- [ ] Hamilton-Jacobi 求解器集成
- [ ] 表面 PDE 弱形式实现
- [ ] 时间步耦合逻辑
- [ ] 长期稳定性测试

---

## 参考资源

### 学术论文
1. **Hansbo et al. (2002)**：Ghost Penalty 原始理论
2. **Burman & Zunino (2014)**：稳定化 CutFEM 标准形式
3. **Olshanskii et al. (2009)**：表面 PDE 的稳定有限元方法

### 开源库参考
- **ngsxfem**：Cut cell 算法 C++ 实现
- **dune-cutfem**：另一个 CutFEM 参考实现
- **scikit-geometry**：几何计算 Python 库

### MOOSE 文档
- [MOOSE InterfaceKernel](https://mooseframework.inl.gov/syntax/InterfaceKernels/)
- [MOOSE Quadrature](https://mooseframework.inl.gov/source/quadrature/)
- [MOOSE UserObject](https://mooseframework.inl.gov/syntax/UserObjects/)

---

## 成功指标

### 阶段1 验收标准
- ✓ 条件数下降至少 50%
- ✓ 与标准有限元方法数值精度一致（$\mathcal{O}(h^2)$）
- ✓ 所有单元测试通过

### 阶段2 验收标准
- ✓ 精确处理被切割单元的积分
- ✓ 支持任意 LevelSet 函数
- ✓ 与 MOOSE 并行框架兼容

### 阶段3 验收标准
- ✓ 能够模拟动态演化的流形
- ✓ 长期模拟稳定性（$t \gg 1$）
- ✓ 发表至少一篇顶刊论文

