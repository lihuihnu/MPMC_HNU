# 单入口测试调度与旧入口映射

`workflow_map.json` 是逐 workflow 的审计清单：记录迁移前事件和路径/分支条件、完整手动输入及默认值、runner/平台矩阵、权限、environment、timeout、原 job 条件、中央 job 映射及语义指纹。不得为了匹配指纹而盲目更新清单；须先审计行为变化。

## 入口语义

- 常规 PR 只有 `pr_incremental_ci.yml` 自动运行。旧独立 PR job 进入这个 run，不再生成各自的空 workflow run。原先已由中央路由调用的 reusable workflow 不改执行体。
- 旧 `workflow_dispatch` 输入（包含 required/default/type/options）、手动签名/发布限制、平台矩阵、artifact、安装包验证、Clapeyron/ThermoPack 固定来源均保留。中央自动 job 不注入手动参数默认值，避免把手动签名/发布操作变成 PR 自动操作。
- 原来只有 PR 入口的门禁显式补参数为空的手动入口，不删除其唯一可执行路径。中央手动入口维持既有核心套件集合；特殊手动验证仍使用对应命名 workflow。
- 不把各旧 workflow 全部变成新增 reusable 调用，避免超过 GitHub 对单调用树唯一 reusable workflow 数量的限制；保留 legacy job 的 needs/output、defaults/env 和每 job 权限。

## 增量与失败

`plan.py` 仅采用同 PR、同 base、祖先可达且中央新版本完整成功的 checkpoint。失败、取消或排队的提交不前移基线；无证据、首次运行、Ready-for-review 回退累计 PR 差异。删除/重命名按旧新路径共同选测。未知可执行路径明确失败，不能静默漏测。首次迁移没有新版本 checkpoint，会验证累计受影响范围，不能为节省本轮时间伪造成功基线。

根 `README.md` 是项目概览、能力矩阵与导航，不被 SW92 科学构建或测试消费。仅修改该文件不选择科学 Gate；中央 impact 中的 workflow/注册测试 inventory/选测检查和最终 Required CI result 仍执行。此前它直接命中 `sw92_profile_c_phase_set`，继而触发共享 topology 依赖闭包；现已移除这一条项目概览绑定。科学专题文档仍保持原 owner，例如 `modules/flash/sw92_profile_c_phase_set.md`；不全局排除 Markdown，也不因一次提交包含 README 而跳过同行源码、共享 fixture、CMake 或 workflow 的必要验证。若未来在根 README 引入可执行科学示例，须重新审计其执行 owner。

`verify_workflows.py` 的 `assert_root_readme_selection` 覆盖四种 PR action、纯概览/治理文档、源码/共享夹具/科学专题/CMake/workflow 混合提交，以及删除/重命名两侧路径的选测并集；混合提交的选择应与其实际科学输入单独变更一致，AD suite 也不能因 README 同行而丢失。该回归由现有 impact job 执行，不新增 workflow 或科学测试入口。

Linux 科学测试默认使用 GitHub 官方 `ubuntu-24.04`，Windows/macOS 保留官方平台。私有 `mpmc_hnu` 仅用于明确白名单中的长时 Gate，以及项目负责人已显式授权的 PETSc/MPI Linux 集成/求解 Gate；当前白名单包含 CPA performance audit 与 `flow_discretization_petsc`。新增私有 runner 消费必须同时更新治理白名单，并给出长时/稳定硬件需求或明确的项目负责人授权。`Required CI result` 汇总失败/取消，不将 skipped 当作已跑过测试。分支保护不在此次修改范围内。

## 后续编写

根可选库入口的消费方契约由 `ad.yml` 中 arithmetic matrix 唯一执行，覆盖默认 AD、thermodynamics/flash/mesh/discretization/flow_discretization/well_discretization 全部 64 种选项组合、库专用 preset 和重新关闭选项。`tests/build/root_libraries/`、相关 CMake 与六个消费程序的公共头闭包进入 arithmetic；不重新运行模块专项 CTest。空间、flow_thermodynamics、flow_discretization、well、well_discretization 消费程序分别链接自己的生产 target，不借用其他消费者的 include 路径。`verify_workflows.py` 从消费源码递归发现项目头，明确将 flow_discretization、well_discretization 命名空间映射到各自的 flow/discretization、well/discretization 模块目录，检查直接/传递选测、无关路径、删除/重命名并集和唯一 owner。现有 checkpoint、未知可执行路径拒绝和 Required result 失败传播机制保持不变。`workflow_map.json` 中 `ad.yml` 指纹对应已审计的 arithmetic 条件消费步骤；扩展消费契约无需再改 workflow，原平台、配置、sanitizer、手动输入、测试命令和结果汇总均保持。

先读 `.github/AGENTS.md` 与 `tests/AGENTS.md`；新增测试接现有 owner Gate，同步依赖规则、入口映射与反例。运行 `python3 .github/ci/verify_workflows.py`（PyYAML 6.0.2）。手动参数或矩阵变更必须更新映射并验证特殊路径。此次入口迁移不等同于全部测试计算已经去重；SW92 跨 workflow 的重复计算仍须单独审计，不能以单入口验收替代。
