# 目录结构调整 - 清理指南

## 背景

根据 MOOSE Framework 的标准（参考 chemistry、heat_conduction 等标准模块），应用程序类应该放在 `include/base/` 和 `src/base/` 目录下。

## 变更

### ✅ 已完成

- ✅ 创建了 `include/base/` 目录
- ✅ 创建了 `src/base/` 目录  
- ✅ 在新位置创建了正确的文件：
  - `include/base/CutFEMApp.h`
  - `src/base/CutFEMApp.C`
- ✅ 更新了所有文档（STRUCTURE.md, FILE_INVENTORY.md）

### 🗑️ 需要删除的旧文件

以下文件现在已过时，应该被删除：

```
modules/cutfem/include/CutFEMApp.h        ← 删除！已移到 include/base/CutFEMApp.h
modules/cutfem/src/CutFEMApp.C            ← 删除！已移到 src/base/CutFEMApp.C
```

## 删除步骤

### 选项 1：使用 Git（推荐）

```bash
cd modules/cutfem

# 删除旧文件
rm include/CutFEMApp.h
rm src/CutFEMApp.C

# 确认删除
git status

# 提交删除
git add -u
git commit -m "Reorganize: Move application class to base/ directory

- Move include/CutFEMApp.h → include/base/CutFEMApp.h
- Move src/CutFEMApp.C → src/base/CutFEMApp.C
- Follow MOOSE Framework standard directory layout
- Update documentation and structure files"
```

### 选项 2：手动删除

在文件浏览器中：
1. 导航到 `modules/cutfem/include/`
2. 删除 `CutFEMApp.h`
3. 导航到 `modules/cutfem/src/`
4. 删除 `CutFEMApp.C`

## 验证

删除后，目录结构应该是这样的：

```
modules/cutfem/
├── include/
│   ├── base/
│   │   └── CutFEMApp.h          ✅ 这里
│   ├── kernels/
│   ├── userobjects/
│   └── utils/
│
├── src/
│   ├── base/
│   │   └── CutFEMApp.C          ✅ 这里
│   ├── kernels/
│   ├── userobjects/
│   └── utils/
│
└── ...
```

**不应该有**：
- ❌ `include/CutFEMApp.h` (在根 include 下)
- ❌ `src/CutFEMApp.C` (在根 src 下)

## 编译验证

删除后编译以确保没有问题：

```bash
cd modules/cutfem
make clean
make -j4 METHOD=opt
```

应该编译成功，没有错误。

## 标准参考

参考 MOOSE 框架的其他标准模块：

- `modules/chemistry/include/base/ChemistryApp.h`
- `modules/heat_conduction/include/base/HeatConductionApp.h`
- `modules/level_set/include/base/LevelSetApp.h`

所有这些都遵循相同的 `include/base/` 和 `src/base/` 结构。

---

**完成此清理后，模块目录结构将完全符合 MOOSE Framework 标准。** ✅
