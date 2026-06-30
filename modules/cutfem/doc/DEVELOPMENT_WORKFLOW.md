# CutFEM Module 开发工作流指南

本文档基于 MOOSE 开发最佳实践，为 CutFEM 模块的开发流程提供完整指导。

---

## 1. 环境设置

### 1.1 初始化 Fork 和克隆

```bash
# Fork MOOSE: GitHub → idaholab/moose → 你的账户
# Clone 你的 fork
git clone git@github.com:<YourUsername>/moose.git
cd moose

# 添加 upstream 远程
git remote add upstream git@github.com:idaholab/moose.git

# 验证远程配置
git remote -v
# origin    git@github.com:<YourUsername>/moose.git (fetch)
# origin    git@github.com:<YourUsername>/moose.git (push)
# upstream  git@github.com:idaholab/moose.git (fetch)
# upstream  git@github.com:idaholab/moose.git (push)
```

### 1.2 激活开发环境

```bash
# MOOSE 提供的环境激活脚本
source ~/.bashrc
env-moose activate    # 或: conda activate moose

# 验证激活
which python
# 应输出 moose conda 环境中的 python 路径
```

### 1.3 编译 CutFEM 模块

```bash
cd moose/modules/cutfem

# 编译（优化模式，并行度 4）
make -C . -j4 METHOD=opt

# 编译调试版本（包含符号表）
make -C . -j4 METHOD=debug

# 清理编译结果
make -C . clean

# 完全重建
make -C . clobber
make -C . -j4 METHOD=opt
```

---

## 2. 开发循环

### 2.1 创建功能分支

```bash
# 同步最新的 upstream
git fetch upstream
git checkout devel
git merge upstream/devel

# 创建功能分支（feature branch）
git checkout -b feature/ghost-penalty-stabilization

# 或者针对 bug 修复
git checkout -b bugfix/gp-jacobian-computation
```

### 2.2 代码开发与格式化

#### 编写新的 Kernel/UserObject

参考模板：
```cpp
#pragma once

#include "Kernel.h"

/**
 * Brief description of the class
 * 
 * Detailed description explaining the mathematical formulation,
 * references, and usage.
 * 
 * References:
 * - Citation (Year) "Title"
 */
class MyKernel : public Kernel
{
public:
  static InputParameters validParams();
  MyKernel(const InputParameters & params);

protected:
  virtual Real computeQpResidual() override;
  virtual Real computeQpJacobian() override;
};
```

#### 代码格式化（必需）

```bash
# 格式化单个文件
clang-format -i src/kernels/MyKernel.C
clang-format -i include/kernels/MyKernel.h

# 格式化整个模块
find modules/cutfem/src modules/cutfem/include -name "*.C" -o -name "*.h" \
  | xargs clang-format -i

# 或使用 Makefile 目标
make -C modules/cutfem format
```

#### 遵循 MOOSE 代码规范

- **错误处理**：使用 `mooseError()` 处理致命错误
- **参数验证**：在 `validParams()` 中添加范围检查
  ```cpp
  params.addRangeCheckedParam<Real>("gamma", 1.0, "gamma >= 0 & gamma <= 1",
    "Stabilization parameter must be in [0,1]");
  ```
- **文档注释**：每个类和方法都要有 Doxygen 注释
- **日志输出**：使用 `mooseInfo()` 而不是 `std::cout`

### 2.3 编译和测试

```bash
# 编译模块
make -C modules/cutfem -j4 METHOD=opt

# 运行所有测试
cd modules/cutfem/test
../../../moose_test-opt

# 运行特定测试
../../../moose_test-opt -i tests/ghost_penalty/test_gp.i

# 运行并显示详细输出
../../../moose_test-opt -i tests/ghost_penalty/test_gp.i --verbose

# 比较输出和 gold 标准
../../../moose_test-opt -i tests/ghost_penalty/test_gp.i --check-input
```

### 2.4 验证代码质量

```bash
# 检查编译警告
make -C modules/cutfem clean
make -C modules/cutfem -j4 METHOD=opt 2>&1 | grep -i warning

# 运行所有测试并检查通过率
cd modules/cutfem/test
python3 -m pytest tests/ -v

# 检查代码风格（可选）
clang-tidy src/kernels/GhostPenaltyKernel.C \
  -- -I./include -I../../framework/include
```

---

## 3. 提交工作流

### 3.1 频繁的本地提交

```bash
# 查看修改状态
git status

# 暂存特定文件
git add src/kernels/GhostPenaltyKernel.C include/kernels/GhostPenaltyKernel.h

# 提交（清晰的提交信息）
git commit -m "Add Ghost Penalty kernel for CutFEM stabilization

- Implements face-based gradient jump penalty (Burman et al. 2010)
- Supports arbitrary polynomial order and parameter gamma in [0,1]
- Includes Jacobian for Newton solver
- References: Larson & Zahedi (2020), Wichrowski (2026)

Fixes: #1234 (if applicable)"

# 查看提交历史
git log --oneline -10
```

### 3.2 同步 Upstream 并 Rebase

在提交 PR 前，确保分支与最新的 upstream 一致：

```bash
# 获取最新的 upstream
git fetch upstream

# Rebase 你的分支到 upstream/devel
git rebase upstream/devel

# 如果有冲突，解决后继续
# 1. 编辑冲突文件
# 2. git add <resolved-files>
# 3. git rebase --continue

# 或者中止 rebase
# git rebase --abort
```

### 3.3 清理并推送

```bash
# 格式化代码（最后一次）
make -C modules/cutfem format

# 最终测试
cd modules/cutfem/test
../../../moose_test-opt

# 推送到你的 fork
git push origin feature/ghost-penalty-stabilization
```

---

## 4. Pull Request (PR) 工作流

### 4.1 在 GitHub 上创建 PR

1. 访问你的 fork：`github.com/<YourUsername>/moose`
2. 点击 "Compare & pull request"
3. 确保：
   - **Base**: `idaholab/moose` / `devel`
   - **Compare**: `<YourUsername>/moose` / `feature/ghost-penalty-stabilization`

### 4.2 填写 PR 模板

```markdown
## Description
Brief description of changes

## Type of Change
- [ ] Bug fix
- [x] New feature
- [ ] Breaking change

## Mathematical Background
Reference key equations and papers:
- Implemented Eq. (1.12) from Larson & Zahedi (2020)
- Ghost Penalty: $$s_{h,F}(w,v) = \sum_j c_j h^{2(j-1+\gamma)} [D_n^j w] [D_n^j v]$$

## Testing
- [x] Added new tests
- [x] All existing tests pass
- [x] Verified convergence order: O(h^2)

## References
- Burman, E. et al. (2015). CutFEM: Discretizing geometry and PDEs.
- Larson, M. G., & Zahedi, S. (2020). Stabilization of high order cut FEM on surfaces.
- Wichrowski, M. (2026). Matrix-free ghost penalty evaluation via tensor product factorization.

Closes #<issue-number> (if applicable)
```

### 4.3 审查过程

- **审查者** (MOOSE CCB 成员): 检查数学正确性、代码质量、文档
- **自动检查**: CIVET CI 系统在 https://civet.inl.gov 上运行所有测试
- **修订**: 根据反馈进行修改并推送更新

```bash
# 根据审查意见进行修改
# 编辑文件...
git add <files>
git commit -m "Address review feedback: [describe changes]"

# 推送更新到同一分支（PR 会自动更新）
git push origin feature/ghost-penalty-stabilization
```

### 4.4 合并

一旦 PR 被批准且所有检查通过：
1. 点击 "Squash and merge" 或 "Rebase and merge"
2. 分支会自动合并到 `devel`
3. 删除特性分支

```bash
# 清理本地分支
git checkout devel
git fetch upstream
git merge upstream/devel
git branch -d feature/ghost-penalty-stabilization
```

---

## 5. 阶段性开发里程碑

### 阶段 1：Ghost Penalty Stabilization （2-3个月）

**目标**：验证 Ghost Penalty 的条件数改进

```bash
# Milestone: Ghost Penalty PR ready for review
git log --oneline feature/ghost-penalty | head -20

# Expected PRs:
# - GhostPenaltyKernel implementation and tests
# - Documentation and examples
# - Convergence studies
```

**成功指标**：
- ✅ 条件数从 $10^4$ 降到 $10^3$（1个数量级改进）
- ✅ 迭代次数减少 60-70%
- ✅ H1 误差 $O(h)$，L2 误差 $O(h^2)$

### 阶段 2：Cut Cell Integration （3-4个月）

```bash
git checkout -b feature/cut-cell-quadrature
# 开发 CutCellQuadratureUserObject...
```

**成功指标**：
- ✅ 正确识别被切割的单元
- ✅ 精确的子单元剖分
- ✅ 最优收敛阶

### 阶段 3：Dynamic Interface Evolution （3个月）

```bash
git checkout -b feature/surface-pde-evolution
# 集成 Hamilton-Jacobi + 表面 PDE...
```

**成功指标**：
- ✅ 界面平稳演化
- ✅ 质量守恒
- ✅ 长期稳定性

---

## 6. 常见问题排查

### 编译错误

```bash
# 清理所有编译产物
make -C modules/cutfem clobber

# 使用完整的编译日志
make -C modules/cutfem -j1 METHOD=opt 2>&1 | tee build.log

# 查看特定错误
grep -i "error" build.log | head -20
```

### 测试失败

```bash
# 运行失败的单个测试
cd modules/cutfem/test
../../../moose_test-opt -i tests/ghost_penalty/failing_test.i

# 与 gold 标准比较
diff tests/ghost_penalty/gold/output.e \
     tests/ghost_penalty/output.e

# 重新生成 gold（谨慎！）
../../../moose_test-opt -i tests/ghost_penalty/test.i --no-skip-check
cp tests/ghost_penalty/output.e tests/ghost_penalty/gold/
```

### Rebase 冲突

```bash
# 查看冲突
git status

# 手动编辑冲突文件（搜索 <<<<<<< 标记）
# 然后:
git add <resolved-files>
git rebase --continue

# 或中止
git rebase --abort
```

---

## 7. 文档更新

每个新功能都需要相应的文档：

### 添加类文档

在 `doc/content/source/kernels/GhostPenaltyKernel.md`:

```markdown
# GhostPenaltyKernel

!syntax description /Kernels/GhostPenaltyKernel

## Overview
Mathematical background and usage...

## Example Input

!listing modules/cutfem/examples/poisson_with_ghost_penalty.i
         block=InterfaceKernels

## Theory

Ghost Penalty stabilization (Burman et al. 2010) ...

## References

- Burman, E., Hansbo, P., Larson, M. G., & Massing, A. (2015)
```

### 更新模块文档

在 `doc/content/modules/cutfem/index.md`:

```markdown
# CutFEM Module

## Description
Cut Finite Element Methods for unfitted mesh discretization...

## Implemented Capabilities

- [x] Ghost Penalty Stabilization (Phase 1)
- [ ] Cut Cell Integration (Phase 2)
- [ ] Dynamic Interface Evolution (Phase 3)

## References

[List key papers]
```

---

## 8. 性能优化与基准

### 性能分析

```bash
# 使用 PAPI 进行性能分析
make -C modules/cutfem -j4 METHOD=oprof

# 运行并生成性能报告
cd modules/cutfem/test
../../../moose_test-oprof -i tests/ghost_penalty/test_gp.i

# 分析结果
perf report
```

### 并行扩展性测试

```bash
# MPI 并行运行
mpirun -n 4 ./cutfem-opt -i examples/poisson_with_ghost_penalty.i

# 监控并行效率
mpirun -n 4 ./cutfem-opt -i examples/poisson_with_ghost_penalty.i \
  -petsc_options "-log_view"
```

---

## 9. 版本发布清单

在准备版本发布时：

```bash
# 更新版本号
cat VERSION  # 查看当前版本

# 更新 CHANGELOG
# 编辑 CHANGELOG.md...

# 创建标签
git tag -a v1.0.0-cutfem -m "Release: Ghost Penalty implementation"
git push origin v1.0.0-cutfem

# 创建 Release notes
# 在 GitHub Releases 页面添加...
```

---

**最后提醒**：CutFEM 是一个复杂的研究项目。始终与 MOOSE 社区沟通，参考官方文档，
并在提交 PR 前进行充分的测试和验证。
