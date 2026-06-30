# CutFEM 模块 - 编译与测试指南

**最后更新**：2026-06-26  
**当前阶段**：Phase 1 - Ghost Penalty 稳定化实现完成

---

## 📋 已实现的内容

### ✅ 核心框架

- [x] **GhostPenaltyKernel.h/C** - 完整实现
  - 梯度跳跃惩罚计算
  - 完整的 Jacobian（4 种类型）
  - 参数化的稳定化强度

- [x] **CutCellQuadratureUserObject.h/C** - 框架实现
  - 元素切割检测框架
  - 缓存数据结构
  - 统计信息收集

- [x] **CutFEMApp.h/C** - 应用注册
  - 模块初始化
  - 对象注册

### 📋 待完成

- [ ] 实际编译验证
- [ ] 运行 Phase 1 测试
- [ ] 调整参数和验证收敛性

---

## 🔧 编译步骤

### 前置条件

```bash
# 假设 MOOSE 已安装在：
export MOOSE_DIR=/path/to/moose

# 验证环境
echo $MOOSE_DIR
ls $MOOSE_DIR/framework/build.mk
```

### Step 1: 清理构建（如果之前编译过）

```bash
cd modules/cutfem
make clean
```

### Step 2: 编译模块

```bash
cd modules/cutfem

# 优化编译（推荐用于生产）
make -j4 METHOD=opt

# 或调试编译（用于开发）
make -j4 METHOD=debug

# 或性能分析编译
make -j4 METHOD=oprof
```

**预期输出**：
```
Compiling C++ (opt)...
[50%] Building C++ object src/kernels/GhostPenaltyKernel.o
[75%] Building C++ object src/userobjects/CutCellQuadratureUserObject.o
[100%] Building C++ object src/base/CutFEMApp.o
Linking module executable...
===========================
Module Built Successfully
===========================
```

### Step 3: 验证编译成功

```bash
# 检查构建是否产生了可执行文件
ls -la lib/
ls -la lib/*cutfem*opt*

# 应该看到类似：
# libcutfem-opt.so (on Linux)
# libcutfem-opt.dylib (on macOS)
```

---

## 🧪 运行测试

### 单个测试

```bash
cd modules/cutfem/test

# 运行 Phase 1 Ghost Penalty 测试
../../../moose_test-opt -i tests/ghost_penalty/test_gp.i

# 或使用 make 运行所有测试
cd modules/cutfem
make test
```

### 验证输出

**预期看到**：

```
=====================================
Poisson Problem with Ghost Penalty
=====================================

Mesh: 20 x 20 QUAD9 elements
Level Set: Circle of radius 0.5

NonlinearSystem:
  Number of DOFs: [large number]
  
Solving with GMRES + ASM preconditioner...
  
Iteration  0: |res|   = 1.234e+02
Iteration  1: |res|   = 1.432e+01
Iteration  2: |res|   = 3.456e-02
Iteration  3: |res|   = 1.234e-08  <- converged

Postprocessors:
  L2_error = 2.345e-03
  max_u    = 0.1234
  min_u    = -0.0001

Solution written to: test_gp_out.e
```

---

## 📊 验证数值精度

### 条件数改进检查

启用 KSP 监控以查看求解器统计：

```bash
# 修改 test_gp.i，在 [Executioner] 部分添加：
# petsc_options = '-ksp_monitor'

../../../moose_test-opt -i tests/ghost_penalty/test_gp.i -v 2>&1 | grep -A20 "KSP"
```

**预期结果**：
- 没有 Ghost Penalty：迭代次数 > 100
- 有 Ghost Penalty：迭代次数 < 50（改进 > 50%）

### 收敛性研究

对不同网格大小进行收敛性测试：

```bash
# 编写脚本进行收敛性研究
cat > convergence_study.sh << 'EOF'
#!/bin/bash

echo "Mesh Size    H1_Error    L2_Error    Iterations"
echo "=========================================="

for nx in 10 20 40 80; do
  ny=$nx
  
  # 修改网格大小运行测试
  ../../../moose_test-opt -i tests/ghost_penalty/test_gp.i \
    Mesh/nx=$nx Mesh/ny=$ny \
    -v 2>&1 | grep -E "(H1_error|L2_error|Iteration)"
    
  echo "---"
done
EOF

chmod +x convergence_study.sh
./convergence_study.sh
```

**预期收敛阶**：
- H1 误差：$O(h)$ - 头部系数应随 $h$ 线性递减
- L2 误差：$O(h^2)$ - 头部系数应随 $h^2$ 递减

---

## 🔍 调试和故障排除

### 编译错误

**错误：找不到 MooseApp.h**
```
fatal error: MooseApp.h: No such file or directory

解决方案：
  1. 检查 MOOSE_DIR 设置
  2. 检查 include 路径配置
  3. 验证 MOOSE 是否正确编译
```

**错误：InterfaceKernel 未找到**
```
fatal error: InterfaceKernel.h: No such file or directory

解决方案：
  1. 检查 Makefile 中的依赖项
  2. InterfaceKernel 在 MOOSE 框架中
  3. 验证 MOOSE 版本 (需要 2021+)
```

### 运行时错误

**错误：Level Set 变量未找到**
```
Cannot find auxiliary variable 'phi'

解决方案：
  1. 在输入文件中定义 [AuxVariables] 部分
  2. 在 [AuxKernels] 中计算 Level Set 函数
  3. 确保 AuxVariable 执行时机正确
```

**错误：条件数未改进**
```
GMRES 迭代次数仍然很多

检查清单：
  1. Ghost Penalty kernel 是否正确添加到问题中？
  2. gamma, k, c_F 参数值是否合理？
  3. Level Set 函数是否正确定义？
  4. 网格大小是否合适（不能太粗）？
```

### 性能问题

**症状：编译缓慢**
```
使用并行编译：
  make -j8 METHOD=opt  (增加 j 值)

关闭优化加快调试编译：
  make -j8 METHOD=debug
```

**症状：求解缓慢**
```
检查线性求解器设置：
  - 是否使用了合适的预处理器？
  - 是否启用了矩阵无关 (matrix-free) 选项？

参考 test_gp.i 中的 [Executioner] 配置
```

---

## 📈 性能监控

### 使用 PAPI 进行性能分析

```bash
# 编译时启用 PAPI
LIBPAPI_EVENTS="PAPI_FP_OPS,PAPI_L1_DCA,PAPI_L2_DCA" \
  moose_test-opt -i test_gp.i
```

### 内存使用情况

```bash
# 监控内存使用
/usr/bin/time -v moose_test-opt -i test_gp.i 2>&1 | grep Memory
```

---

## 📚 参考文档

| 文档 | 内容 |
|------|------|
| [IMPLEMENTATION_STRATEGY.md](./IMPLEMENTATION_STRATEGY.md) | 完整的技术细节和数学公式 |
| [DEVELOPMENT_WORKFLOW.md](./doc/content/workflow.md) | Git 工作流和代码审查过程 |
| [test_gp.i](./test/tests/ghost_penalty/test_gp.i) | Phase 1 示例输入文件 |

---

## ✅ 完整编译检查表

- [ ] 设置 `MOOSE_DIR` 环境变量
- [ ] 进入 `modules/cutfem` 目录
- [ ] 运行 `make -j4 METHOD=opt` 编译
- [ ] 验证没有编译错误
- [ ] 运行 `moose_test-opt -i test/tests/ghost_penalty/test_gp.i`
- [ ] 验证测试完成且收敛
- [ ] 检查 L2 误差值
- [ ] 验证求解器迭代次数
- [ ] 查看输出文件 test_gp_out.e

---

## 🚀 下一步

1. **完成 Phase 1 验证**
   - 收敛性研究（3 个网格大小）
   - 参数灵敏度分析
   - 与解析解对比

2. **准备发布**
   - 更新文档
   - 编写论文或技术报告
   - 提交 PR 到 MOOSE

3. **启动 Phase 2**
   - 实现 Marching Cubes 算法
   - 集成 Level Set 模块
   - 开发非贴体积分

---

**祝编译和测试顺利！** 🎉

如有问题，参考 [doc/content/](./doc/content/) 中的详细文档。
