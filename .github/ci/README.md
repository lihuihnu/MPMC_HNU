# 单入口测试调度与旧入口映射

`workflow_map.json` 是逐 workflow 的审计清单：记录迁移前事件和路径/分支条件、完整手动输入及默认值、runner/平台矩阵、权限、environment、timeout、原 job 条件、中央 job 映射及语义指纹。不得为了匹配指纹而盲目更新清单；须先审计行为变化。

## 入口语义

- 常规 PR 只有 `pr_incremental_ci.yml` 自动运行。旧独立 PR job 进入这个 run，不再生成各自的空 workflow run。原先已由中央路由调用的 reusable workflow 不改执行体。
- 旧 `workflow_dispatch` 输入（包含 required/default/type/options）、手动签名/发布限制、平台矩阵、artifact、安装包验证、Clapeyron/ThermoPack 固定来源均保留。中央自动 job 不注入手动参数默认值，避免把手动签名/发布操作变成 PR 自动操作。
- 原来只有 PR 入口的门禁显式补参数为空的手动入口，不删除其唯一可执行路径。中央手动入口维持既有核心套件集合；特殊手动验证仍使用对应命名 workflow。
- 不把各旧 workflow 全部变成新增 reusable 调用，避免超过 GitHub 对单调用树唯一 reusable workflow 数量的限制；保留 legacy job 的 needs/output、defaults/env 和每 job 权限。

## 增量与失败

`plan.py` 仅采用同 PR、同 base、祖先可达且中央新版本完整成功的 checkpoint。失败、取消或排队的提交不前移基线；无证据、首次运行、Ready-for-review 回退累计 PR 差异。删除/重命名按旧新路径共同选测。未知可执行路径明确失败，不能静默漏测。首次迁移没有新版本 checkpoint，会验证累计受影响范围，不能为节省本轮时间伪造成功基线。

Linux 科学测试默认使用 GitHub 官方 `ubuntu-24.04`，Windows/macOS 保留官方平台。私有 `mpmc_hnu` 仅用于明确白名单中的长时 Gate；当前仅 CPA performance audit 使用私有 runner。新增私有 runner 消费必须同时更新治理白名单并给出运行时间或稳定硬件需求证据。`Required CI result` 汇总失败/取消，不将 skipped 当作已跑过测试。分支保护不在此次修改范围内。

## 后续编写

先读 `.github/AGENTS.md` 与 `tests/AGENTS.md`；新增测试接现有 owner Gate，同步依赖规则、入口映射与反例。运行 `python3 .github/ci/verify_workflows.py`（PyYAML 6.0.2）。手动参数或矩阵变更必须更新映射并验证特殊路径。此次入口迁移不等同于全部测试计算已经去重；SW92 跨 workflow 的重复计算仍须单独审计，不能以单入口验收替代。
