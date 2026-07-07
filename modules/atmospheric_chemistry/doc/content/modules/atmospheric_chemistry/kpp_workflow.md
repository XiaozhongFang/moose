# KPP Workflow for Atmospheric Chemistry Box Models

## Overview

The Kinetic Pre-Processor (KPP) generates optimized C code for chemical
ODE systems from a high-level `.kpp` mechanism description. This module uses
KPP-generated shared libraries (`.so`) as the runtime backend for box-mode
chemistry integration.

Execution flow:

```
.kpp mechanism file
        │
        ▼
  kpp/build/Makefile ──► kpp_build_<mech>/
  (KPP tool)                 ├── libkpp_<mech>.so   (shared library)
                             └── kpp_<mech>.json    (species metadata)
        │
        ▼
  KppBoxIntegrator (dlopen .so, dlsym adapter functions)
    ├── kpp_init()        — 初始化 KPP 全局状态
    ├── kpp_set_conc()    — 设置浓度数组 C[]
    ├── kpp_integrate()   — 调用 Rosenbrock/SDIRK/Runge-Kutta 积分器
    └── kpp_get_conc()    — 读取积分后的浓度
        │
        ▼
  MCMBoxModel::execute() — 每个时间步：
    1. 从 ScalarVariable 读取浓度
    2. 调用 KppBoxIntegrator::solve(t0, t1, C)
    3. 写回 ScalarVariable + 解向量
```

## Quick Start (手动测试)

### 1. 确认 KPP 可用

```bash
# 检查 KPP 工具
which kpp || echo "需要安装 KPP"
# 检查编译支持
grep KPP_ENABLED modules/atmospheric_chemistry/Makefile
```

如果 KPP 工具未安装：
```bash
git clone https://github.com/KineticPreProcessor/KPP.git ~/KPP
cd ~/KPP && make
export PATH=$PATH:~/KPP/bin
```

### 2. 编译带 KPP 支持的 MOOSE

```bash
cd modules/atmospheric_chemistry
KPP_ENABLED=1 make -j4
```

确认编译成功：
```bash
# 验证可执行文件能加载 KPP 库
./atmospheric_chemistry-opt --help 2>&1 | head -1
```

### 3. 构建 KPP 共享库

```bash
cd modules/atmospheric_chemistry
make -f kpp/build/Makefile MECH=test/tests/kpp/kpp_small_strato/small_strato.kpp
```

成功输出：
```
=== Done: test/tests/kpp/kpp_small_strato/kpp_build_small_strato/libkpp_small_strato.so ===
```

验证产物：
```bash
nm -D test/tests/kpp/kpp_small_strato/kpp_build_small_strato/libkpp_small_strato.so | grep -E 'kpp_|Jac_SP'
# 应有: kpp_get_nspec, kpp_init, kpp_integrate, kpp_set_conc, kpp_get_conc, Jac_SP
```

### 4. 手动运行 KPP box 测试

```bash
cd modules/atmospheric_chemistry
cd test/tests/kpp
../../atmospheric_chemistry-opt -i kpp_small_strato.i --no-gdb-backtrace
```

观察输出：
- `KppBoxIntegrator: loaded ... libkpp_small_strato.so` — 库加载成功
- `Solve Converged!` — 每个时间步收敛
- 最终应生成 `kpp_small_strato_out.csv`

### 5. 检查输出

```bash
cd modules/atmospheric_chemistry/test/tests/kpp
# 浓度应有变化（非初始值）
head -5 kpp_small_strato_out.csv 2>/dev/null || head -5 kpp_small_strato.csv
```

比较 gold：
```bash
diff <(head -3 kpp_small_strato_out.csv) <(head -3 gold/kpp_small_strato.csv)
```

## 常见问题

### "Loaded 0 species from KPP mechanism"

MCMBoxModel::initialize() 的 KPP 路径从 ScalarVariable 获取 _n_species。
如果 ScalarVariable 尚未初始化（execute_on 时机问题），会报 0 species。
不影响实际积分——积分器直接操作 KPP 全局数组。

### "KPP integration failed"

积分器报错。先试更宽松的容差：
```
chem_solver_rtol = 1e-4
chem_solver_atol = 1e-6
```
或换求解器类型（kpp_sdirk / kpp_runge_kutta）。

### 浓度没有变化（所有时间步值相同）

KPP 积分器可能未正确初始化。检查：
1. `.so` 是否正确构建（`nm -D` 确认符号）
2. `kpp_init` 是否被调用（日志中应有 `KppBoxIntegrator: loaded`）
3. 时间步长是否合理（dt 太大/太小）

### "requires KPP_ENABLED"

编译时未加 `-DKPP_ENABLED`：
```bash
# 检查当前二进制
strings atmospheric_chemistry-opt | grep -c kpp_rosenbrock
# 如果返回 0，需重新编译
KPP_ENABLED=1 make -j4
```

## 环境变量

| 变量 | 用途 | 示例 |
|---|---|---|
| `KPP_HOME` | KPP 安装根目录（用于 #MODEL 包含） | `/opt/KPP` |
| `KPP_LIB` | 覆盖 KPP .so 路径 | `/path/to/libkpp_mech.so` |
| `KPP_ENABLED` | 编译时启用 KPP | `1` |

## 测试文件

- `test/tests/kpp/kpp_small_strato.i` — KPP box 验证测试
- `test/tests/kpp/kpp_small_strato/` — KPP 机制文件（.kpp/.spc/.eqn）
- `test/tests/kpp/gold/kpp_small_strato.csv` — gold 参考输出

## 更新 Gold 文件

```bash
cd modules/atmospheric_chemistry
make -f kpp/build/Makefile MECH=test/tests/kpp/kpp_small_strato/small_strato.kpp
cd test/tests/kpp
../../atmospheric_chemistry-opt -i kpp_small_strato.i --no-gdb-backtrace
cp kpp_small_strato_out.csv gold/kpp_small_strato.csv
```
