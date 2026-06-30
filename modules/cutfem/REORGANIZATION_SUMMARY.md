# 目录结构调整总结

**调整时间**：2026-06-26  
**原因**：根据 MOOSE Framework 标准，应用程序类应该在 `base/` 目录下  
**参考**：chemistry、heat_conduction、level_set 等标准模块

---

## 前后对比

### ❌ 旧结构（不符合规范）

```
modules/cutfem/
├── include/
│   ├── CutFEMApp.h               ← 应用类不应该在这里
│   ├── kernels/
│   ├── userobjects/
│   └── utils/
└── src/
    ├── CutFEMApp.C               ← 应用类不应该在这里
    ├── kernels/
    ├── userobjects/
    └── utils/
```

### ✅ 新结构（符合 MOOSE 规范）

```
modules/cutfem/
├── include/
│   ├── base/                     ← ⭐ 新增
│   │   └── CutFEMApp.h          ← 正确位置
│   ├── kernels/
│   ├── userobjects/
│   └── utils/
└── src/
    ├── base/                     ← ⭐ 新增
    │   └── CutFEMApp.C          ← 正确位置
    ├── kernels/
    ├── userobjects/
    └── utils/
```

---

## 已完成的操作

### 1️⃣ 创建新目录
✅ `include/base/`  
✅ `src/base/`

### 2️⃣ 复制文件到新位置
✅ `include/base/CutFEMApp.h`  
✅ `src/base/CutFEMApp.C`

### 3️⃣ 更新文档
✅ STRUCTURE.md  
✅ FILE_INVENTORY.md  
✅ CLEANUP_GUIDE.md (新增)

### 4️⃣ 待做：删除旧文件
⏳ 删除 `include/CutFEMApp.h`  
⏳ 删除 `src/CutFEMApp.C`

（见 CLEANUP_GUIDE.md）

---

## MOOSE 框架标准

标准的 MOOSE 模块结构：

```
modules/<module_name>/
├── include/
│   ├── base/             # 应用类、基础类
│   ├── kernels/
│   ├── bcs/
│   ├── ics/
│   ├── materials/
│   ├── userobjects/
│   ├── postprocessors/
│   ├── outputs/
│   ├── functions/
│   └── utils/
│
├── src/
│   ├── base/
│   ├── kernels/
│   ├── bcs/
│   ├── ics/
│   ├── materials/
│   ├── userobjects/
│   ├── postprocessors/
│   ├── outputs/
│   ├── functions/
│   └── utils/
│
├── test/
│   └── tests/
│
├── examples/
│
└── doc/
    └── content/
```

**关键点**：
- ✅ 所有应用程序类都在 `base/` 目录中
- ✅ 所有其他物理对象按类型组织
- ✅ 头文件和源文件镜像结构
- ✅ 文档在 `doc/content/`

---

## 参考模块的例子

### chemistry 模块
```
modules/chemistry/
├── include/base/ChemistryApp.h
└── src/base/ChemistryApp.C
```

### heat_conduction 模块
```
modules/heat_conduction/
├── include/base/HeatConductionApp.h
└── src/base/HeatConductionApp.C
```

### level_set 模块
```
modules/level_set/
├── include/base/LevelSetApp.h
└── src/base/LevelSetApp.C
```

---

## 下一步

### 立即完成
1. 按照 [CLEANUP_GUIDE.md](./CLEANUP_GUIDE.md) 删除旧文件

### 编译验证
```bash
cd modules/cutfem
make clean
make -j4 METHOD=opt
```

### 提交
```bash
git add -A
git commit -m "Reorganize: Move application class to base/ directory"
```

---

## 验收检查表

- ✅ 创建了 `include/base/` 目录
- ✅ 创建了 `src/base/` 目录
- ✅ CutFEMApp.h 在 `include/base/`
- ✅ CutFEMApp.C 在 `src/base/`
- ✅ 更新了所有相关文档
- ⏳ 删除旧的 `include/CutFEMApp.h`
- ⏳ 删除旧的 `src/CutFEMApp.C`
- ⏳ 验证编译成功

---

**状态**：✅ 结构调整完成，⏳ 待清理旧文件

参考 CLEANUP_GUIDE.md 了解删除步骤。
