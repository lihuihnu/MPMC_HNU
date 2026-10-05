# CI 编写与治理规范

本文件补充根 AGENTS.md，适用于 .github 内所有自动化。规则约束现有中央入口及专项 workflow；执行状态以实际 DAG、作业和日志为准，不得把规范本身当作完成证据。

## 自动触发与测试所有权

1. 完整功能的 Ready PR 验证只允许一个中央自动入口：`pr_incremental_ci.yml`。不得新增独立的 pull_request、push 或 workflow_run 测试入口；不使用“先创建几十个 workflow，再让测试 job skipped”冒充单入口。已有专项 workflow 迁移时保留必要的手动入口、参数和特殊事件语义，不得机械删除唯一入口。
2. 新测试接入现有 Gate，不为每个切片新建 workflow。必要下游由中央路由计算集合并集。同一 `(Gate, 平台, 编译器, 构建配置, sanitizer/验证模式)` 每轮只有一个执行 owner；不得在多个上层 job 里重复 configure/build/ctest 同一下游。相同参数的 oracle 生成也应去重，不同验证模式不得混同。
3. 影响分析必须覆盖源码、公共头、CMake、共享 fixture、参考生成脚本、数据文件、协议以及必要的传递依赖。删除与重命名同时考虑旧、新路径。未知影响明确失败或保守扩大相关验证，禁止静默不测。
4. PR 云端验收始终基于整个功能的累计差异 `merge-base(base, HEAD)..HEAD`，包括 Ready 后的集中修复；不再查询 PR 增量 checkpoint，Draft 的全 skipped 调度记录不能成为基线。原非 PR push checkpoint 语义保留。治理脚本必须覆盖 Draft 不调度、Ready/修复累计选测及原失败传播。
5. 纯治理文档不触发无关科学编译；执行命令、测试集合或公共依赖变更必须验证适用矩阵。不以“CI 改动”名义统一跳过实际受影响测试。

## 验收与安全

- 保留测试名、目标、断言、参考值、数值容差、严格告警及适用 GCC sanitizer / Clang / MSVC / macOS 组合。不通过删测试、改变 golden 或忽略失败换取速度。
- Linux 编译/数值测试默认使用 GitHub 官方 `ubuntu-24.04`；Windows 使用 `windows-2022`、macOS 使用既有官方托管 runner。私有 `mpmc_hnu` 只允许进入显式白名单：当前白名单为 `flow_discretization_petsc.yml` 的 PETSc 3.19.6/MPI 2-rank 集成/求解 Gate，以及 `cpa_performance_audit.yml` 的经审计长时 paired performance audit。前者属于已授权 PETSc/MPI Gate，后者属于非 PETSc 长时 Gate；除此之外的普通编译、单元测试、矩阵回归和 Android 工具链验证继续使用官方 runner。新增私有-runner owner 必须先审计运行时间、依赖和可移植性并更新本白名单；workflow 必须显式安装/核验所需工具，不能依赖机器偶然状态。
- 默认最小只读权限，固定 Action 提交，合理 timeout 与 superseded-run cancellation。临时迁移/预检入口完成或失败后清理，不作为长期产品代码保留。
- 必须报告 selected、未选中、实际执行、失败、取消和未验证状态；workflow success 不等于其 CTest 执行通过。最终汇总不得把失败、取消或缺失验收计为成功；不擅改分支保护。
- 开发切片在本地形成完整差异并验证，Draft 下批量提交同步；完整功能满足本地验收后才转 Ready 触发云端多平台/科学性验收。失败后收齐证据集中修复，不逐文件 push，不在无关 pending 验证期间追加清理提交。
- CI 变更必须补充单入口、直接/间接命中、无关路径、未知路径、删除/重命名、去重、checkpoint 与失败传播的回归；以实际 DAG 和 job 日志验收，不能只数 YAML 中出现一次的字符串。

## 本地模式

本地是开发期默认验证场所；小切片通过后同步 Draft PR，完整功能完成才转 Ready 进行云端验收。遵守同样的断言、容差和测试 ownership；使用 [本地指南](../docs/development.md)。现有 GitHub runner 白名单约束云端自动化，不禁止开发机运行原有测试。本地能力不足项明确列为功能级云端验收待办，不为每个切片自动启动云端补验；不删测试、不降低科学判据或绕过保护。

本地可以承担大算例实际运行；云端按资源条件保留正确性与可行性所需的代表性部分。规模不同不降低断言、容差或必要覆盖，具体规则统一见 [根 AGENTS 的资源分工](../AGENTS.md#本地与云端的资源分工)。
