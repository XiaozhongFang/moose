# Atmospheric Chemistry 标准 App 接入设计

## 目标

让 `atmospheric_chemistry` 与其他标准 MOOSE 模块一样，可由 Stork 生成的外部 App 通过
`ATMOSPHERIC_CHEMISTRY := yes` 完整启用。完整启用包括模块及其传递依赖的构建、MOOSE
对象与 App 类型注册，以及 KPP/SUNDIALS 构建配置的一致传播。

## 设计选择

采用 MOOSE 已有的模块同名构建片段机制。新增
`modules/atmospheric_chemistry/atmospheric_chemistry.mk`，由 `framework/app.mk` 现有的
`-include $(APPLICATION_DIR)/$(APPLICATION_NAME).mk` 自动读取。模块自身构建和外部 App
构建因此使用同一份配置，不修改 framework，也不复制 KPP/SUNDIALS 逻辑。

KPP 保持当前行为，启用 atmospheric chemistry 时默认定义 `KPP_ENABLED`。SUNDIALS
继续根据 `SUNDIALS_DIR`（默认 `CONDA_PREFIX`）中的头文件自动检测。

## 修改范围

1. `modules/atmospheric_chemistry/atmospheric_chemistry.mk`
   保存 KPP include 路径、`KPP_ENABLED`、SUNDIALS 检测、链接参数和
   `HAVE_SUNDIALS`。路径以稳定的 atmospheric chemistry 模块目录为基准，不能依赖
   `app.mk` 返回后会恢复的 `APPLICATION_DIR`。
2. `modules/atmospheric_chemistry/Makefile`
   删除已经迁入同名 `.mk` 的重复配置，仅保留模块顶层构建入口。
3. `modules/modules.mk`
   当 `ATMOSPHERIC_CHEMISTRY=yes` 时传递启用 `NAVIER_STOKES`。将 atmospheric
   chemistry 的注册段移动到 Navier-Stokes 及其依赖之后，并保留
   `DEPEND_MODULES := navier_stokes`。
4. `modules/module_loader/include/ModulesApp.h` 与
   `modules/module_loader/src/ModulesApp.C`
   在 `ATMOSPHERIC_CHEMISTRY_ENABLED` 条件下包含 `AtmosphericChemistryApp.h`，并在
   `registerAllObjects`、`registerAll` 和 `registerApps` 中注册该模块。
5. `stork/Makefile.app`
   在标准模块开关列表中加入 `ATMOSPHERIC_CHEMISTRY := no`。
6. `scripts/stork.sh`
   补充新模块接入提示，明确 Stork 模板和 module loader 也是标准模块接入清单的一部分。

## 数据与构建流程

Stork App Makefile 设置 `ATMOSPHERIC_CHEMISTRY=yes` 后包含 `modules/modules.mk`。
`modules.mk` 先启用 Navier-Stokes 及其传递依赖，再依序解析依赖模块和 atmospheric
chemistry。解析 atmospheric chemistry 的 `app.mk` 时自动加载
`atmospheric_chemistry.mk`，使模块库和最终外部 App 链接获得相同的 KPP/SUNDIALS
设置。编译宏 `ATMOSPHERIC_CHEMISTRY_ENABLED` 随后使 `ModulesApp` 注册 atmospheric
chemistry 的对象、Action、语法和 App 类型。

## 验证策略

先增加并运行失败的集成回归检查，覆盖以下不变量：

- Stork 模板公开 atmospheric chemistry 开关。
- 启用 atmospheric chemistry 后，Make 变量数据库同时启用 Navier-Stokes 及其依赖。
- 外部 App 构建变量包含正确的 atmospheric chemistry KPP include 路径和
  `KPP_ENABLED`；检测到 SUNDIALS 时还包含对应编译与链接参数。
- module loader 的模板注册、兼容注册和 App 注册路径都包含 atmospheric chemistry。

实现后重新运行回归检查，并进行 atmospheric chemistry 模块构建、combined 构建，
以及临时 Stork 外部 App 的构建和对象语法查询。现有 atmospheric chemistry 回归测试
用于确认模块行为没有退化。

## 非目标

- 不迁移 MAS1998 benchmark 源码或输入数据。
- 不改变 KPP 默认启用或 SUNDIALS 自动检测策略。
- 不重构其他模块，也不修改 `framework/app.mk`。
