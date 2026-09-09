# PT 闪蒸：顶层未确定原因摘要

## 审阅与范围

基线为 PR #13 的 `20fbb52477b0868b5a7ec54927573161b73b9551`，tree `35cff1861ea2cf62fdaa5c998ba16d5d2895fac7`。先读取根 AGENTS、PR13 差异、独立露点参考/测试及相分裂/稳定性结果路径，提交原实现者的 COMMENT review `5153138414` 后才编写本增量；不是第三方批准，不合并或修改前序分支。

PR13 说明露点固定状态的液相尝试已满足守恒/逸度，只因小相量低于原门槛而被拒绝。原顶层统一文本不能呈现这一事实，也不能区分初始/最终 TPD 和分裂预算。此增量只增强 `PtSplitResult::diagnostic`，不新增状态枚举、结果字段或机器协议。

## 摘要规则

- 初始 TPD 未确定：保留阶段前缀，并汇总各 trial 的终止状态及次数。参考物性失败单独表述为 `reference property failure; trial searches not started`，不把驱动填充的失败占位 trial 说成已运行的搜索。
- 无可接受两相候选：固定顺序列出全部 attempt 终止类型及次数；`phase_disappearance` 另说明 `balance and fugacity tolerances met; phase fraction at or below minimum_phase_fraction`。混合失败同时列出，不让第一/最后一条尝试覆盖其他原因。
- 既有 `attempt_limit_reached` 标志下，分别检查“尝试次数配额阻止继续”和“分裂物性调用配额阻止继续”。零次尝试不能自动归类为初始化失败。
- 最终 TPD 未确定：保留 `split equations converged; final phase-set stability is indeterminate` 前缀，再列出最终 trial 终止计数。已选候选继续保留；摘要不把它升级为已接受平衡相集合。

摘要使用已有终止状态，而非历史 `property_issue` 或被拒步的历史字符串。历史拒绝、具体根错误和全部物性信息仍在原低层结构中；不会为了拼接摘要复制回调可提供的任意长文本。只扫描有限种类的状态，输出固定标签和计数，不拼接组成数组或每步轨迹。

例如，当一个尝试被小相门槛拒绝、另一个回溯失败时，摘要包含：

```text
no acceptable two-phase candidate; attempts=2;
phase_disappearance=1 (balance and fugacity tolerances met;
phase fraction at or below minimum_phase_fraction);
line_search_failed=1; final phase-set stability not evaluated
```

上例仅为可读性换行，实际字符串为单行。字段名指向调用时记录的 options，文本不硬编码默认门槛。`phase_disappearance` 不等于物理相严格不存在，不证明单相稳定。`no_interior_rr_root` 是固定 K 求解结果，不自动解释成热力学单相。

## 行为不变量与兼容性

生产变化仅为 `pt_split.hpp` 的两个内部只读格式化函数及三个诊断赋值。没有改变输入校验、数值表达式、收敛/失稳/相分率门槛、初值、回溯、预算、物性调用、attempt 顺序、候选选择或任何状态转换。非表示尺度与较高 Gibbs 等原本已有具体原因的路径保持原文；成功路径和额外相失稳路径也不改。

`candidate()` 与 `equations_converged()` 保持原语义：没有被选择的可接受候选时，即使底层相消失尝试已满足方程，顶层仍无候选/返回 false。数值可通过 `attempt.point` 查看，不能以新诊断替代已有验收流程。

诊断文本供人阅读，不是版本稳定的解析协议。机器分流继续读取状态、attempt/trial 终止类型、配额字段和可选候选；不建议匹配整条英文句子。字符串采用既有标准库，不新增依赖；分配异常仍按原方式传播，不吞掉程序错误。

## 验证范围（实际通过情况以具体提交 Actions 为准）

新增独立 `diagnostic_summary_test.cpp` 的五项 CTest：所有尝试原因与混合/无工作计数、TPD 终止/历史拒绝/参考失败/未开始搜索、PR13 露点及门槛上方/外侧实际结果、分阶段预算与回调计数、参考失败与顺序复用。制造终止记录只测试文本，不冒充流体解。真实露点和预算调用复用原模型参数与只读参考，没有新参考数值或拟合。

本轮 `pt_split.hpp` 为下游协调层变化：运行原 52 项相分裂/边界/停滞/露点集成及新增五项，保留原全部 Decimal 参考步骤。使用既有官方 GCC Debug+ASan/UBSan、Clang Release、MSVC Release，检查公共头模板、字符串和告警的可移植性。共享 `pt_stability.hpp`、PR76、RR、旧测试与独立参考均不变，不重跑上游独立22项稳定性或无关AD/热力学旧套件。

本次不是新的数值修复，也没有三相、含水、全局稳定认证、实验验证、性能或并发结论。露点科学依据与门槛含义沿用 [PR13 审计](dew_limit_audit.md)，完整汽液接口见 [PT 相分裂契约](pt_split.md)。
