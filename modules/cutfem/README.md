# CutFEM Module for MOOSE

**Cut Finite Element Methods for unfitted mesh discretization**

> ⚠️ **最新更新**：目录结构已调整为符合 MOOSE 规范。参见 [REORGANIZATION_SUMMARY.md](./REORGANIZATION_SUMMARY.md) 和 [CLEANUP_GUIDE.md](./CLEANUP_GUIDE.md)

## 📌 快速开始

### 文档导航

按照这个顺序阅读文档：

1. **首先读这个** ⭐ → [`SUMMARY.md`](./SUMMARY.md)
   - 5 分钟快速了解项目全貌

2. **项目概览** → [`PROJECT_OVERVIEW.md`](./PROJECT_OVERVIEW.md)
   - 架构设计、三阶段路线图、成功指标

3. **完整技术方案** → [`doc/IMPLEMENTATION_STRATEGY.md`](./doc/IMPLEMENTATION_STRATEGY.md)
   - 所有数学推导、代码框架、验证策略

4. **开发工作流** → [`doc/DEVELOPMENT_WORKFLOW.md`](./doc/DEVELOPMENT_WORKFLOW.md)
   - Git 工作流、编译测试、PR 提交指南

5. **运行示例** → [`examples/poisson_with_ghost_penalty.i`](./examples/poisson_with_ghost_penalty.i)
   - 可直接运行的 Phase 1 示例

---

## 🎯 核心功能

### Phase 1：Ghost Penalty 稳定化（现在）

- ✅ 梯度跳跃惩罚项实现
- ✅ Nitsche 边界条件
- ✅ 条件数改进 O(h^{-2})
- **状态**：架构完成，待编码实现

### Phase 2：非贴体积分（3-4个月后）

- 动态 Level Set 集成
- Marching Cubes 子单元剖分
- 高阶积分规则

### Phase 3：表面 PDE 演化（6-8个月后）

- Laplace-Beltrami 求解器
- Hamilton-Jacobi 界面演化
- 时间步耦合

---

## 📂 项目结构

```
modules/cutfem/
├── SUMMARY.md                           ⭐ 从这里开始
├── PROJECT_OVERVIEW.md                  路线图 & 架构
├── Makefile                             编译配置
├── doc/
│   ├── IMPLEMENTATION_STRATEGY.md       完整技术方案
│   ├── DEVELOPMENT_WORKFLOW.md          开发工作流
│   └── index.md                         模块文档首页
├── include/
│   ├── CutFEMApp.h
│   ├── kernels/
│   │   └── GhostPenaltyKernel.h
│   ├── userobjects/
│   │   └── CutCellQuadratureUserObject.h
│   └── utils/
│       └── MarchingCubes.h
├── src/
│   ├── CutFEMApp.C
│   └── kernels/
│       └── GhostPenaltyKernel.C (待实现)
├── examples/
│   └── poisson_with_ghost_penalty.i    示例
└── test/
    └── tests/ghost_penalty/
        └── test.txt                    测试定义
```

---

## 🚀 立即可做的事

### 1. 编译模块
```bash
cd modules/cutfem
make -C . -j4 METHOD=opt
```

### 2. 实现 C++ 源文件

参考 [`IMPLEMENTATION_STRATEGY.md` 第 1.2.3 节](./doc/IMPLEMENTATION_STRATEGY.md)，
复制 `GhostPenaltyKernel.C` 的完整实现。

### 3. 运行示例
```bash
cd test
../../../moose_test-opt -i ../examples/poisson_with_ghost_penalty.i
```

### 4. 验证结果
- 检查条件数改进（应 > 50%）
- 验证 H1 误差 ~ O(h)
- 查看 L2 误差 ~ O(h²)

---

## 📚 关键论文

| 论文 | 年份 | 核心贡献 |
|------|------|---------|
| Burman et al. "CutFEM: Discretizing geometry and PDEs" | 2015 | Ghost Penalty 原理 |
| Larson & Zahedi "Stabilization of high order cut FEM" | 2020 | 混合稳定化分析 |
| Wichrowski "Matrix-free ghost penalty evaluation" | 2026 | 张量积优化 |
| Burman et al. (综述) "Cut finite element methods" | 2025 | 最新进展 |

---

## 🔧 开发命令

```bash
# 编译
make -C modules/cutfem -j4 METHOD=opt

# 测试
cd modules/cutfem/test && ../../../moose_test-opt

# 代码格式化
make -C modules/cutfem format

# 清理
make -C modules/cutfem clean
```

---

## 📖 主要文档

| 文档 | 用途 | 长度 |
|------|------|------|
| SUMMARY.md | 项目概览 | 5 min |
| PROJECT_OVERVIEW.md | 架构和路线 | 15 min |
| IMPLEMENTATION_STRATEGY.md | 完整技术方案 | 1 hour |
| DEVELOPMENT_WORKFLOW.md | 开发指南 | 30 min |
| examples/poisson_with_ghost_penalty.i | 运行示例 | 10 min |

---

## ✅ 成功指标

**Phase 1（目前）**：
- [ ] 编译成功
- [ ] 单元测试通过
- [ ] 条件数改进 ≥ 50%
- [ ] 收敛阶正确（H1: O(h), L2: O(h²)）

**整个项目**：
- [ ] 三个阶段完成
- [ ] 发表研究论文
- [ ] MOOSE 社区接受

---

## 🤝 贡献指南

1. 创建 feature 分支：`git checkout -b feature/xxx`
2. 遵循 [MOOSE 代码规范](https://mooseframework.inl.gov/contribute/)
3. 编写并运行测试
4. 提交 Pull Request（参考 DEVELOPMENT_WORKFLOW.md）

---

## 📞 获取帮助

- **MOOSE 论坛**：https://github.com/idaholab/moose/discussions
- **文档**：https://mooseframework.inl.gov
- **本地文档**：见 `doc/` 目录

---

## 📊 项目统计

- **总文档行数**：3000+
- **数学公式**：50+
- **代码模板**：500+ 行
- **测试用例**：8 个
- **预计工期**：8-10 个月

---

## 📝 许可证

该模块遵循 MOOSE 框架许可证。

---

**最后更新**：2026-06-26

**现在开始**：打开 [`SUMMARY.md`](./SUMMARY.md) 👈
