# 免登录 PR76 自定义模型计算

## 当前开发重点

软件不需要账号、注册、登录、OAuth/OIDC 或 IdP。当前路径直接使用桌面本地计算，
不再以 Web identity 或部署平台作为热力学工作的前置条件。

- C++ 负责 PR76 物性/逸度、初始相稳定性搜索、两相/三相分裂与最终相集合复核。
- 前端只编辑模型定义、求解设置和 P/T/z，不包含 EOS、TPD、RR/RR3 或闪蒸公式。
- Electron 主进程拥有本地模型 session、原生 model reference 和 transport；renderer
  只获得版本化 `apply / solve / release / cancel` workbench 能力、权威模型快照及结果。
- renderer 不获得 model handle/token、session ID、connect/reconnect、原生端口、
  子进程凭据或 transport 对象。
- 本地子进程仍使用私有临时 token 保护同机 IPC；这是进程传输保护，不是用户登录。

普通 Web/Android 继续保留已有 PT Flash 入口；没有新增登录页、公开原生服务端口或
身份依赖。

## 实际使用

1. 打开桌面软件，进入 **PR76 Expert**；直接进入本地工作台，无需账号或建立用户会话。
2. 在左侧填写模型标识、参数来源与组分。支持新增、删除、重排组分；填写摩尔质量
   [kg/mol]、Tc [K]、Pc [Pa]、偏心因子及每个无序组分对的显式 kij。缺失 kij
   不会被当成零。选择已提供的 preset，或在完整后端设置快照上修改 CUSTOM 设置。
3. 点击 **Create immutable model** 应用模型。主进程创建新的不可变后端模型，并在
   成功接管后释放旧模型；renderer 不接触任何后端引用。
4. 在右侧输入 P [Pa]、T [K]、按当前模型组分顺序排列的进料摩尔分数 z，点击
   **Compute PR76 flash**。初始值均为空，不填充猜测数据，不静默归一化进料。
5. 阅读已接受相数、各相摩尔相分率 beta、各相组分摩尔分数 x，以及可用的 Z、
   原生诊断与相数转移证据。应用的模型参数和来源仍可展开查看。

修改参数、kij、组分顺序或设置后，旧结果立即隐藏；先应用新模型，再计算。
应用新模型会重新建立空的温压/进料输入，避免新增/删除/重排组分后沿用错位的 z。
修改 P/T/z 也会清除旧结果。每个工作台最多一个在途模型调用，不排队、不自动重试。
取消、窗口关闭、renderer 崩溃或所有权不明确时，主进程按保守策略回收整个模型 epoch，
而不是猜测某个后端引用仍然安全。

## 稳定性与结果语义

`ExpertOwnedModel` 在 renderer 中只是“当前主进程模型”的本地代际包装，不再持有
backend token。`snapshot` 来自成功 `apply` 的权威 `ModelSnapshot`；`solve()` 通过
workbench 让主进程调用当前不可变模型，并保留完整 `FullPtResult`。

`candidatePhaseSet` 存在或含三相，不代表该相集合已经被接受。只有原生 outcome
为 ACCEPTED 才显示“Accepted phase count”；INDETERMINATE/PHASE_SET_UNSTABLE 的
相数据放在明确标识的诊断候选区。前端不根据 Z、根序号或数组顺序推测油/气/水相。

PR76 capability 声明执行初始稳定性搜索和最终复核，但 `globalStabilityProven=false`
时不能称为数学全局稳定性证明。界面展示后端 diagnostic 与 transitionReport；
workbench 边界不会把候选收敛、根数或小残差重新解释成相集合接受。

本地输入检查仅处理十进制/科学记数法、有限正温压、组分身份/顺序和单个摩尔分数
范围。归一化容差、适用范围、active support、求解门槛由 C++ 模型判断。字段错误
使用已有 typed code/field；任意原生异常文本、session/token/handle 不进入 renderer。

## 所有权与回归

成功应用模型会提升 renderer 本地 generation。旧 React 对象随后执行 `release()`
只做本地失效，不会误释放已经由主进程接管的新模型。当前对象释放才会调用 workbench
`release`。失败 apply 会使旧 renderer ownership 一并失效，因为主进程已经保守清理
该模型 epoch；不会恢复可能已经陈旧的对象。

生产 `preload.cjs` 只包含通用 PT bridge 和无 handle 的 model workbench，源码中不存在
低级 create/describe/model-token/connect/reconnect bridge。旧的低级 model desktop bridge
被隔离到专用 `modelDebugPreload.cjs`，只由 native IPC 回归 harness 临时加载；产品打包
只复制生产 preload，且产品主进程不注册该低级 IPC。

聚焦前端/桌面回归覆盖：动态组分对齐、空/非法输入、原样传递进料、单请求限制、
迟到结果/错误、模型替换、旧对象释放、窗口清理、完整结果克隆、字段级错误映射和
preload API 表面。UI/序列化测试不是物性验证。

科学侧另有一条配置链路回归：复用已冻结的 Li–Firoozabadi 六组分 sour-gas PR76
文献基准，通过运行时 `ThermodynamicModelDefinition -> Pr76ExecutableModel` 重建模型，
再执行稳定性/最多三相求解，并与直接 PR76 backend 的完整结果逐字段比较。它验证
“用户可编辑参数/组分 -> 动态 C++ 模型 -> 三相结果”没有绕开科学内核；独立三相数值
oracle 仍由原 #94 Decimal/物理混合物回归负责。

## 数值边界

本 workbench 增量不改 PR76 EOS、根选择、TPD、两相/三相方程、默认容差、相分率门槛
或 `indeterminate` 语义。科学算法的后续改进必须由独立物理/数值证据驱动，不能为了
前端成功率删除小相、放宽门槛或把有限稳定性搜索描述为全局证明。
