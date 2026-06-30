# 目录结构检查清单

最后更新：2026-06-26

## 核心目录验证

### ✅ include/ 目录结构

```
include/
├── base/                           ✅ 新建
│   └── CutFEMApp.h                ✅ 存在
├── kernels/                        ✅ 存在
│   └── GhostPenaltyKernel.h       ✅ 存在
├── userobjects/                    ✅ 存在
│   └── CutCellQuadratureUserObject.h ✅ 存在
└── utils/                          ✅ 存在
    └── (待实现)
```

**需要验证**：
- [ ] `include/base/CutFEMApp.h` 文件存在
- [ ] `include/kernels/GhostPenaltyKernel.h` 文件存在
- [ ] `include/userobjects/CutCellQuadratureUserObject.h` 文件存在

### ✅ src/ 目录结构

```
src/
├── base/                           ✅ 新建
│   └── CutFEMApp.C                ✅ 存在
├── kernels/                        ✅ 新建
│   └── (待实现)
├── userobjects/                    ✅ 新建
│   └── (待实现)
└── utils/                          ✅ 新建
    └── (待实现)
```

**需要验证**：
- [ ] `src/base/CutFEMApp.C` 文件存在
- [ ] `src/kernels/` 目录存在
- [ ] `src/userobjects/` 目录存在
- [ ] `src/utils/` 目录存在

### ✅ test/ 目录结构

```
test/
└── tests/
    ├── ghost_penalty/
    │   ├── test_gp.i              ✅ 存在
    │   ├── tests.txt              ✅ 存在
    │   └── test_gp_gold/          ✅ 目录（空）
    ├── cut_cells/
    └── surface_pde/
```

**需要验证**：
- [ ] `test/tests/ghost_penalty/test_gp.i` 存在
- [ ] `test/tests/ghost_penalty/tests.txt` 存在

### ✅ doc/ 目录结构

```
doc/
└── content/
    ├── index.md                   ✅ 存在
    ├── overview.md                ✅ 存在
    ├── implementation.md          ✅ 存在
    └── workflow.md                ✅ 存在
```

**需要验证**：
- [ ] `doc/content/` 目录存在
- [ ] 所有 markdown 文件存在

### ✅ examples/ 目录

```
examples/
└── poisson_with_ghost_penalty.i   ✅ 存在
```

**需要验证**：
- [ ] `examples/poisson_with_ghost_penalty.i` 存在

### 📋 根级文件

```
✅ README.md
✅ STRUCTURE.md
✅ FILE_INVENTORY.md
✅ REORGANIZATION_SUMMARY.md       ⭐ 新增
✅ CLEANUP_GUIDE.md                ⭐ 新增
✅ Makefile
✅ .gitignore
✅ .clang-format
✅ IMPLEMENTATION_STRATEGY.md
✅ PROJECT_OVERVIEW.md
✅ SUMMARY.md
```

## 🗑️ 需要清理的旧文件

⚠️ **以下文件应该被删除**：

```
❌ include/CutFEMApp.h              （已移到 include/base/CutFEMApp.h）
❌ src/CutFEMApp.C                  （已移到 src/base/CutFEMApp.C）
```

参见 [CLEANUP_GUIDE.md](./CLEANUP_GUIDE.md)

## 📝 文档更新状态

| 文档 | 更新状态 | 更新内容 |
|------|---------|---------|
| STRUCTURE.md | ✅ 已更新 | 反映 base/ 目录 |
| FILE_INVENTORY.md | ✅ 已更新 | 反映新位置 |
| REORGANIZATION_SUMMARY.md | ✅ 已创建 | 说明调整原因 |
| CLEANUP_GUIDE.md | ✅ 已创建 | 删除旧文件步骤 |
| README.md | ✅ 已更新 | 指向新指南 |

## 🔍 快速验证

运行以下命令验证目录结构：

```bash
# 验证 base 目录存在且有正确的文件
ls -la include/base/
ls -la src/base/

# 验证旧文件（应该显示为空或不存在）
ls include/CutFEMApp.h
ls src/CutFEMApp.C

# 验证文档
ls doc/content/
```

## ✅ 完成检查表

### 结构验证
- [ ] `include/base/CutFEMApp.h` 存在
- [ ] `src/base/CutFEMApp.C` 存在
- [ ] 所有其他 include/ 子目录存在
- [ ] 所有其他 src/ 子目录存在
- [ ] test/tests/ 结构正确
- [ ] doc/content/ 存在所有文档

### 清理步骤
- [ ] 根据 CLEANUP_GUIDE.md 删除旧文件
- [ ] 验证编译成功：`make -j4 METHOD=opt`
- [ ] 运行测试：`make test`

### 提交步骤
- [ ] 执行 `git add -A`
- [ ] 执行提交：`git commit -m "Reorganize: Move application class to base/ directory"`
- [ ] 推送到分支

## 📚 相关文档

- [REORGANIZATION_SUMMARY.md](./REORGANIZATION_SUMMARY.md) - 为什么进行调整
- [CLEANUP_GUIDE.md](./CLEANUP_GUIDE.md) - 如何删除旧文件
- [STRUCTURE.md](./STRUCTURE.md) - 完整的目录结构说明
- [FILE_INVENTORY.md](./FILE_INVENTORY.md) - 文件清单

## ✨ 现在的目录结构

现在完全符合 MOOSE Framework 标准：

```
modules/cutfem/          ← 模块根目录
├── include/base/        ← ⭐ 应用程序类
├── include/kernels/     ← 核头文件
├── include/userobjects/ ← 用户对象头文件
├── include/utils/       ← 工具类头文件
├── src/base/            ← ⭐ 应用程序实现
├── src/kernels/         ← 核实现
├── src/userobjects/     ← 用户对象实现
├── src/utils/           ← 工具类实现
├── test/tests/          ← 测试
├── examples/            ← 示例
└── doc/content/         ← 文档
```

---

状态：✅ 结构调整完成，⏳ 待清理和编译验证
