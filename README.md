# MPMC_HNU

面向多相、多组分计算的模块化高性能计算平台，采用可移植的 C++20 计算后端与独立 React/TypeScript 前端。计算核心可以作为库无界面运行。

当前已有独立 AD、PR76/SW92/CPA PT 闪蒸后端、局部灵敏度与部分 physics closure，以及已合并 PR #118 建立的可变相数 fully implicit 非等温多组分流动链：自然变量 1/2/3 相、组分/能量守恒、TPFA Darcy/重力、相渗/毛细压契约、PETSc SNES/GMRES/ASM、相变重启、自适应时间步和 Peaceman 单井控制；同时已有 PT 服务、Web/Electron 前端和 Android 工程应用。安装包候选和工程预览仍不等于正式发行。具体能力与证据范围见下表，测试结果以对应提交的 [GitHub Actions](https://github.com/lihuihnu/MPMC_HNU/actions) 为准。

## 阅读入口与文档分工

本页维护项目概览、当前能力与证据矩阵、依赖边界和里程碑；接口、公式、参数、构建命令及验证详情在所属模块维护，历史审计保留在专题文档中。能力变化时同步更新所属模块与本页对应行，避免在多个章节重复粘贴整段状态。

| 需要了解的内容 | 入口 |
| --- | --- |
| 开发、审计、测试与提交规则 | [AGENTS.md](AGENTS.md) |
| AD 类型、初等函数、固定及运行期 Jacobian | [AD 模块](modules/ad/README.md) |
| 组分/参数契约、PR76、SW92、CPA 物性 | [热力学模块](modules/thermodynamics/README.md) |
| 相稳定性、两相/三相、统一后端、验证与历史审计 | [闪蒸模块](modules/flash/README.md) |
| 局部物性快照、导数与每总流体体积组分库存 | [Physics closure](modules/physics/README.md)、[component inventory](modules/physics/component_inventory.md) |
| 网格 topology/geometry/field/DoF、Gmsh/VTU/GRDECL、PETSc 与 TPFA ownership 边界 | [Mesh 模块](modules/mesh/README.md) |
| 自然变量多相流、守恒装配、PETSc/SNES、相变与时间推进 | [Flow 模块](modules/flow/README.md) |
| Peaceman 井模型、fixed-BHP / total-rate 控制与控制状态机 | [Well 模块](modules/well/README.md) |
| 服务、通信、进程与部署 | [Runtime](modules/runtime/README.md)、[API](api/README.md)、[gRPC adapter](modules/runtime_grpc/README.md)、[process host](modules/pt_process/README.md)、[Envoy edge](deploy/pt-grpc-web/README.md) |
| Web、桌面和 Android 使用/构建边界 | [前端](frontend/README.md)、[桌面/原生产品](products/pt/README.md)、[Android](products/pt_android/README.md) |
| 项目级工程参考与原始文献 | [参考资料](docs/references.md)；具体公式、参数来源仍以模型专题为准 |

## 当前计算能力

以下“当前能力—证据—限制—里程碑”矩阵以 **2026-09-30 的 [main@10ea62aa](https://github.com/lihuihnu/MPMC_HNU/commit/10ea62aabe658eedd542bdbf363adfd735f70663)** 为本轮核对基线。未合并 PR 不计入主线能力；M0–M5 是下文定义的后续验收方向，不是完成标记。根 README 维护项目级矩阵，模块文档维护具体契约与来源。

证据按用途区分：**契约/数值**验证接口、不变量、解析导数或独立数值参考；**物理**对照明确体系与工况的实验或文献数据；**工程**验证平台、通信、装配、求解和产品集成；**性能**需要同环境、同工作量的重复测量。它们不能相互替代。制造解与 synthetic fixture 必须显式标注；从待测实现生成的期望值不能称为独立参考。

| 层次 | 当前能力 | 证据与复现入口 | 限制／未验证范围 | 对应里程碑 |
| --- | --- | --- | --- | --- |
| AD | 独立 `Dual<T,N>`、初等函数、固定与运行期分块 Jacobian；嵌套 forward AD 与固定维 `value_gradient_hessian` | 契约/数值：[模块说明](modules/ad/README.md)、[Jacobian/Hessian 注册入口](tests/ad/jacobian/CMakeLists.txt)、[解析 Hessian 回归](tests/ad/jacobian/nested_hessian_test.cpp)；算术、math、runtime 各有独立入口 | 仅依赖标准库；未提供反向模式、动态 Hessian API 或稀疏导数容器；AD 类型支持嵌套不表示每个 EOS 都支持嵌套 | M1：按实际内核核对导数证据；M4：有测量后优化 |
| PR76 PT | 有序参数、纯/混合/PT 物性、有限 TPD、汽液 PT、最多三相编排、continuation 和统一 backend | 契约/独立数值/物理：[三相说明](modules/flash/pr76_three_phase.md)、[测试入口](tests/flash/pr76_three_phase/CMakeLists.txt)；Li–Firoozabadi 六组分独立数值三相参考、Heringer 3→2 边界及二元/三元参考 | classical vdW one-fluid mixing、显式常数对称 `kij`；有限搜索不是全局稳定性证明；[盲测 #102](https://github.com/lihuihnu/MPMC_HNU/pull/102) 的模型偏差保留，不证明宽工况实验精度或通用相形态分类 | M1：审计模型适用性；M5：有来源时扩模型 |
| SW92 PT | corrected-original；fixed-family stability/VLE、Whitson dual-model observables、Xu-style max2、Profile-C authoritative 1/2/3-phase PT | 数值/物理：[分 profile 入口](modules/flash/README.md)、[Sample-6 审计](modules/flash/sw92_authoritative_three_phase_audit.md)、[实验三相线](modules/flash/sw92_three_phase_experimental_validation.md)及其[测试入口](tests/flash/sw92_three_phase_experimental_line/CMakeLists.txt) | 各 profile 的假设与结论独立；固定 NaCl molality 不等于盐库存守恒；H 相仍为 `nonaqueous_unclassified`，family/root 不自动给出跨单元物理相身份 | M1：保持 profile/参考一致；M2：闭合流动 |
| CPA PT | explicit-site-pair 参数、缔合、PT 物性、minimum-Gibbs stability、VLE、max3 和统一 backend | 物理：[甲醇(2B)/水(4C) 两相验证](modules/flash/cpa_physical_validation.md)及[测试入口](tests/flash/cpa_physical_validation/CMakeLists.txt)；结构：[max3 测试](tests/flash/cpa_max3/CMakeLists.txt) | max3 仍为 synthetic structural validation；缺兼容的物理 VLLE 参数与独立三相 oracle，不能据两相结果宣布物理三相已验证 | M1：标清证据缺口；M5：资料齐备后扩展 |
| 灵敏度与 physics | PR76 内部两相、SW92 Profile-C 固定相集合的局部隐式导数及 closure；model-neutral component inventory/local Jacobian | 契约/数值：[Physics](modules/physics/README.md)、[库存契约](modules/physics/component_inventory.md)、[库存测试入口](tests/physics/component_inventory/CMakeLists.txt) | 导数绑定接受的光滑相/根/family 分支；PR76 standalone 单相 closure、CPA flash sensitivity/closure 尚未实现；`mol/m³ fluid` 库存不等于孔隙体积累积项 | M1：单位、分支与导数核对；M5：按模型需求补齐 |
| Mesh / discretization | 声明的 1D/2D/3D topology/geometry、fields、DoF/partition；Gmsh 4.1 ASCII、VTU ASCII、GRDECL 子集；DMPlex/PetscSection/PetscSF；TPFA admissibility/transmissibility、owned connection schedule 和 MPIAIJ 符号预分配 | 契约/工程：[Mesh 范围](modules/mesh/README.md)、[core](tests/mesh/core/CMakeLists.txt)、[外部格式兼容](tests/mesh/external_compatibility/CMakeLists.txt)、[PETSc/MPI](tests/mesh/petsc/CMakeLists.txt)、[离散测试](tests/discretization/core/CMakeLists.txt) | Mesh core 不依赖 PETSc/MPI；TPFA 只处理 admissibility 允许的面；GRDECL fault/pinch/NNC、MPFA/非正交处理及任意复杂网格不属于已验证能力 | M1：离散与网格证据；M4：规模；M5：按算例扩展 |
| Flow / PETSc 通用链 | 1/2/3 相自然变量；Backward Euler 组分/能量守恒；TPFA Darcy/重力、上风组分/焓通量与导热；ragged residual/Jacobian；Newton/BT → GMRES/restricted ASM；相变外层 rebuild/rebind、自适应时间步与 accepted history/time commit | 契约/数值/工程：[计算链](modules/flow/README.md)、[core](tests/flow/core/CMakeLists.txt)、[离散](tests/flow_discretization/core/CMakeLists.txt)、[PETSc owning CTest](tests/flow_discretization/petsc/CMakeLists.txt)；PR76 有真实物性非零通量 transient；[sour-gas 回归](tests/flow_discretization/petsc/pr76_li_firoozabadi_sour_gas_short_step.md)是单元、无内部面、零井率的定常短步 | 小规模集成不等于场尺度验证；不同 EOS 的 property/transition 覆盖不对等；未提供扩散/弥散、反应、地质力学、裂缝或 AMR 的完整生产链 | M1：收敛与独立判据；M3：正式算例；M4：规模效率 |
| SW92 flow property / restart | frozen selected-phase 1P/2P/3P value/AD bridge；零盐 CO₂/H₂O 有来源的输运/热量 provider；authoritative sidecar 驱动同 dt 重建与重求解 | 数值/工程：[契约 §49–52](modules/flow/contracts.md)、[物性回归](tests/flow/core/sw92_co2_water_properties_test.cpp)、[CO₂/H₂O restart](tests/flow_discretization/petsc/sw92_transactional_phase_transition_restart_test.cpp)、[Sample-6 restart](tests/flow_discretization/petsc/sw92_transactional_phase_transition_sample6_test.cpp)，由现有 Flow core/PETSc 入口拥有 | `pc=none`；CO₂/H₂O provider 为零盐、300–1200 K，来源完整不等于全域实验精度；restart 拒绝 authoritative faces 与 ghost overlap。Sample-6 使用真实 SW92 平衡、synthetic molar masses 与 manufactured transport/caloric properties，且显式启用 storage anchor；不是物理输运验证，也不自动获得井耦合或 accepted physical-time commit | **M2：family/root-aware 跨相数输运、ghost 与相变闭环** |
| Well | Cartesian/diagonal-K Peaceman WI；fixed-BHP 单/多 completion owner-only source 与可变相数 rebind；单井 total-molar-rate 全局 BHP unknown/四块 Jacobian；minimum-BHP、跨步 accepted control、滞回及相变事务重启 | 契约/数值/工程：[Well](modules/well/README.md)、[core](tests/well/core/CMakeLists.txt)、[离散](tests/well/discretization/CMakeLists.txt)、共享 [PETSc owning CTest](tests/flow_discretization/petsc/CMakeLists.txt) | 单井控制基线；无 multi-well network/priority、井筒压降/热损失、surface/phase-rate、maximum-BHP、completion schedule 或设施模型；不能自动套用到尚未闭合的 SW92 restart 路径 | M1：守恒与控制判据；M3：配置运行；M5：多井等扩展 |

最近一轮相关执行证据为 [CI #257](https://github.com/lihuihnu/MPMC_HNU/actions/runs/36689406584)，对应 #127 的提交 `915ae0ad`：Flow core 与 Flow discretization 的 GCC/Clang/MSVC 作业、PETSc 两进程作业及 Required CI result 成功；已核对 GCC 日志的 52/52、68/68 和 PETSc 1/1。PETSc 的一个 CTest 内聚合多个回归场景，不能把 CTest 数当作科学覆盖数量。该 run 只证明选中的受影响范围；合并提交 `10ea62aa` 没有另跑整仓 CI，矩阵列出的其余测试入口也不表示本轮已重新执行。

三个预置 PT 后端均通过 [统一能力与结果契约](modules/flash/pt_flash_backend.md) 暴露 1/2/3 相能力、版本与有序组分身份。预置 backend inventory 本身仍是冻结快照；组分增减、替换、重排不能原地修改该快照。Electron 的免登录 **PR76 Expert** 工作台另行支持用完整的新参数快照创建不可变运行时 PR76 模型，可新增、删除、重排组分并显式给出 `Tc`、`Pc`、偏心因子、摩尔质量、完整 `kij` 与 solver settings，再由同一 C++ 相稳定/最多三相内核求解。这个能力不是通用物性数据库；内建 PR76 甲烷/乙烷/丙烷、SW92 CO₂/淡水、CPA 甲醇/水快照仍只覆盖声明的窄文献体系。

### PR76 验证边界与历史证据

严格 PR76 当前生产路径继续使用 classical vdW one-fluid mixing 与显式常数对称 `kij`。其 pure/mixing/root/`ln(phi)`/TPD/two-phase/max-three-phase/boundary/continuation/publication 链路已经在已合并的 [PR #112](https://github.com/lihuihnu/MPMC_HNU/pull/112) 做系统审计；当时在所覆盖的数学与数值契约内未保留未修复的 production implementation defect，该历史结论不排除后续发现缺陷。完整公式、接口和当前测试入口由 [热力学模块](modules/thermodynamics/README.md) 与 [闪蒸模块](modules/flash/README.md) 维护，根 README 不再复制逐 PR 的历史过程。

实现正确性不等于模型具有宽体系、宽温压范围的实验预测精度。独立盲测 [PR #102](https://github.com/lihuihnu/MPMC_HNU/pull/102) 仍保留为严格 PR76 的模型能力受阻证据：其 M-40 三元相变压力相对实验偏差约 `13.97%`，不能通过放宽生产容差、拟合三元目标或静默切换模型形式来“修复”。高级混合规则探索仍限定在测试/审计侧；任何生产级 WS/NRTL 等扩展都必须重新冻结公式、导数、参数来源、有效范围和独立验证。

## 应用与产品

| 路径 | 当前状态与入口 |
| --- | --- |
| Hosted Web | React 经版本化 Protobuf/gRPC-Web、Envoy 和 C++ adapter 调用 `PtService`；部署使用显式 HTTPS origin、mTLS 与外部身份配置。详见 [部署说明](deploy/pt-grpc-web/README.md)。该可选部署路径的边缘安全配置不是本地桌面工作台的登录前置条件。 |
| Windows/macOS/Linux native staging | 锁定 Conan binary graph，以 dependency seed + restore-only staging 生成可搬移 host。vcpkg 路径保留用于本地 source build。详见 [产品构建](products/pt/README.md)。 |
| Electron desktop | 复用 React 与 native staging，经 sandbox preload/IPC、ephemeral loopback transport token 和原生 gRPC 调用后端；PR76 Expert 直接本地使用，无账号、注册或登录。renderer 只得到 `apply/solve/release/cancel` 工作台能力与权威快照/结果，不得到原生 model handle、session ID、connect/reconnect 或 transport。跨平台产物仍为未签名工程预览。详见 [桌面前端](frontend/README.md)。 |
| 原生 host 安装包 | 有 Linux DEB、Windows MSI、macOS PKG 候选打包 gate；其载荷为 native host。详见 [host installer](products/pt/installer/README.md)。 |
| Windows desktop 安装器 | 有完整 Electron MSI 安装/启动/卸载 gate，以及固定发布身份、手动受保护的签名 RC workflow；真实签名材料缺失时停止，公开发行仍未完成。详见 [desktop installer](products/pt/desktop_installer/README.md)。 |
| Android | 复用 React，经 Capacitor/JNI 直接调用 C++ `PtService`；universal debug APK 包含 `arm64-v8a` 与 `x86_64`，已有 emulator gate。详见 [Android Product Shell](products/pt_android/README.md)。 |

各适配层传输和展示后端结果，保留 `accepted`、`phase_set_unstable`、`indeterminate` 与服务/通信错误的区别，不重复求解 EOS、守恒、逸度或相数。真实签名、许可、发行、更新与部署条件按产品文档独立验收，不能用工程 smoke 代替科学验证。

## 架构与技术边界

- **计算核心：** C++20、目标级 CMake 与 CTest；`ad` 只依赖标准库。热力学不得反向依赖 flash，公共领域接口不得泄漏 Protobuf、HTTP、UI 或 PETSc/MPI 类型。
- **依赖方向：** `runtime -> flash -> thermodynamics`；`physics` 消费 flash/sensitivity。`runtime_grpc`、`pt_process`、产品壳和前端位于外层，配置后端对象图或适配传输。参数和模型选择放在应用边界，不进入逐标量内循环。
- **前端与通信：** React + TypeScript + Vite；Protobuf + gRPC，浏览器经 gRPC-Web。大体量网格/场结果拟采用独立 HTTP 二进制分块，vtk.js 仅为待评估可视化方案。任务/进度、重试幂等、背压和新增双向通信须在实际实现前单独审计；当前 PT RPC 不代表完整任务系统。
- **网格与求解链：** `core` 最小基础层与 `numerics` 领域无关算法层是架构规划，目前尚无这两个独立模块目录；`mesh` 负责 topology/geometry/field/DoF/I/O/partition，`discretization` 负责 TPFA admissibility/transmissibility 与 PETSc symbolic structure；`flow`/`flow_discretization` 在其上拥有组分/能量守恒、Darcy/热通量、variable-cardinality 自然变量和时间推进，PETSc bridge 负责 MPIAIJ 数值装配与 SNES/KSP/PC 求解。基础模块公共领域接口仍不得泄漏 PETSc/MPI 类型。
- **可移植性：** 目标覆盖 Linux、Windows、macOS，Android 有独立产品路径；具体架构和编译器以实际 CI 证据为准。OpenMP、MPI、PETSc、GPU 和 `SolverBackend` 扩展须独立审计，不成为基础数值测试的前置条件。不默认启用 `fast-math`、`-march=native`，不承诺跨编译器 C++ ABI 兼容或未经测量的性能收益。

## 科学与结果解释

输入采用明确 SI 单位（压力 Pa、温度 K）、摩尔组成与摩尔相分率。模型版本、参数来源、单位、适用域和水处理方式必须可追溯；缺失参数不自动设零，不静默切换 EOS。

| 水处理模式 | 约定 |
| --- | --- |
| `immiscible_water` | 水独立处理，禁止其与其他组分跨该相交换；分别保持物性和质量守恒。 |
| `mutual_solubility` | 允许迁移的水和其他组分在同一守恒与平衡问题中联合求解，不能在无水闪蒸后拼接固定水相。 |

这两种模式是不同物理假设，不代表每个当前后端都支持任意模式。相稳定性、组分守恒、组成归一化和适用的逸度/化学势平衡必须联合验收；不能用小残差、根数、相分率裁剪或两个独立两相结果替代相集合判断。水是组分，水富相是相态，平衡分配不等于有限速率传质。

所有当前有限搜索结果保持 `global_stability_proven=false`。`candidate` 或方程收敛不等于结果已接受；`indeterminate` 必须保留。局部 Jacobian 不跨相出现/消失、根/family/topology 切换或临界退化外推。详细停止阈值、状态枚举和失败处理仅在对应模块契约维护。

## 构建、验证与待完成方向

根 [CMakeLists.txt](CMakeLists.txt) 和 [CMakePresets.json](CMakePresets.json) 默认仍只配置 AD，可按需启用 thermodynamics、flash、mesh、discretization，以及 flow_discretization / well_discretization 桥接库。flash 自动带入 thermodynamics 和 AD，discretization 自动带入 mesh；flow_discretization 带入 flow、thermodynamics、discretization 与 mesh，并提供 flow_thermodynamics 消费目标；well_discretization 带入 well 和完整 flow 桥接依赖链。各选项可同时启用。选项、preset 与下游消费示例见 [根入口用法](tests/build/root_libraries/README.md)。根 CTest 仍只注册 AD 算术，启用库不会自动执行模块专项测试；其他模块、PETSc 桥接、服务与产品保留各自独立入口。配置日志明确报告库集合和测试范围，不能把根目录构建误认为整个项目的验证。

开发流程以 [AGENTS.md](AGENTS.md) 为唯一规则入口：先审计，做可回退增量，再按受影响依赖选择测试。正式自动化中，普通 Linux 编译、单元测试和常规矩阵默认使用 GitHub 官方 `ubuntu-24.04`，Windows/macOS 使用对应官方托管 runner；当前私有 `mpmc_hnu` 白名单仅覆盖 `flow_discretization_petsc.yml` 的 PETSc/MPI 2-rank 集成/求解 Gate，以及 `cpa_performance_audit.yml` 的经审计长时 paired performance audit。workflow 在私有 runner 上仍必须显式安装/核验固定依赖，不依赖持久机器的偶然环境。纯文档核对内容、相对引用和差异范围，代码变更运行必要增量及受影响下游。解析/数值回归、synthetic fixture、实验验证与性能基准分别报告，不把已有工作流当作已经通过的证据。

| 工程层次 | 当前能力 | 证据与复现入口 | 限制／缺口 | 对应里程碑 |
| --- | --- | --- | --- | --- |
| 构建与 CI | 模块独立构建；一个常规自动入口按影响范围选测，保留手动/特殊 Gate | 工程：[根构建](CMakeLists.txt)、[CI 入口映射与规则](.github/ci/README.md)、[workflow verifier](.github/ci/verify_workflows.py) | 根构建只覆盖 AD；尚无统一的可选模块构建入口；当前根 README 路由仍选择 `sw92_profile_c_phase_set`，纯文档选测需另行审计；skipped/未执行不是通过，测试量不是覆盖率 | M0：明确入口、审计选测并另行补齐构建编排 |
| 可维护性与效率 | Flow 已分层并共享校验/坐标算法，混合相数装配已有职责重构 | 工程/局部性能：[已合并 #130](https://github.com/lihuihnu/MPMC_HNU/pull/130)、[已合并 #127](https://github.com/lihuihnu/MPMC_HNU/pull/127) | #130 的局部库存组装与分配测量不是端到端求解加速；核对基线时 #126/#128/#129 尚未合并，不计入能力；缺统一完整算例的时间、内存和并行效率基线 | M0：审查待合并工作；M4：按实际瓶颈优化 |
| PT 应用与运行时模型 | `PtService`、Web/Electron/Android 工程链；Electron PR76 Expert 用新参数快照重建不可变模型，支持组分增删/重排 | 工程：[前端与测试命令](frontend/README.md)、[运行时模型](modules/model_configuration/README.md)、[native 产品](products/pt/README.md)、[Android](products/pt_android/README.md)；产品路径见上表 | 预置 backend 仍是冻结快照；专家输入不是通用物性数据库；当前 PT RPC/UI 不代表完整流动任务、场结果分析或正式发行 | M3：计算任务基础；M5：流动产品与发行 |
| 研究算例与恢复 | 已有库级时间推进、相变重建及测试夹具；网格格式可读写 | 工程：[Flow 调用链](modules/flow/README.md)、[PETSc 回归入口](tests/flow_discretization/petsc/CMakeLists.txt)、[Mesh I/O](modules/mesh/README.md) | 尚未形成配置驱动的完整流动算例产品与经过验收的持久化 checkpoint/restart；相变 same-dt restart 不等于进程中断恢复，网格 I/O 不等于场结果产品 | M3：可复现输入、运行、输出与恢复 |

### 里程碑与验收

开发目标是形成**可独立配置运行、物理依据可追溯、支持恢复计算并能测量并行性能的研究级多相多组分模拟器**。下面列出依赖和完成判据，不承诺未经评估的日期或完成百分比。当前仅推进 M0 的文档矩阵切片；建立矩阵不等于完成整个 M0。后续编码、合并、算例运行或发布仍按 [AGENTS.md](AGENTS.md) 在具体任务授权下推进。

| 里程碑 | 前置条件 | 工作与交付 | 完成判据 | 当前状态 |
| --- | --- | --- | --- | --- |
| **M0 清晰基线** | 核对主分支、适用规则和未合并工作 | 维护本矩阵；审查处理已有独立 PR；整理构建/验证入口，后续单独实现可选模块构建编排 | 每项能力能追溯契约、证据类型、验证 owner 与限制；保留 AD 等基础模块独立构建，项目入口能明确选择和报告验证范围 | **进行中：本次仅文档矩阵**；构建编排和待合并 PR 另行验收 |
| **M1 可信度基线** | M0 的范围和证据入口明确 | 选择代表算例，分别核对公式/单位、导数、离散/求解与模型适用性；在运行前确定误差判据 | 对适用算例分别报告守恒误差、Jacobian 独立核对、网格/时间步收敛和实验或独立数值偏差；保留失败与适用域，#102 不靠放宽容差变绿 | 待完成项目级基线；已有模块回归不能代替全部判据 |
| **M2 SW92 实际流动闭环** | M1 中所选零盐 CO₂/H₂O 物性与算例证据足够；明确相身份及 absent-phase 扩展契约 | 从来源完整的简单体系建立 family/root-aware 跨相数面输运、ghost 同步与事务式相变，接入时间推进 | 非零通量、多单元、MPI 分区及相出现/消失情况下，独立导数核对、组分/能量守恒与最终 history/time commit 均通过；未知物性仍明确拒绝 | **受阻项明确**：当前 restart 禁止 faces/ghost；Sample-6 manufactured 数据不计作物理完成证据 |
| **M3 可运行研究工具** | M1 的目标算例已验证；采用 SW92 相变流动时还需 M2 | 正式算例配置、初始化、边界/井配置、时间推进、结果输出与持久化 checkpoint/restart | 新算例无需修改内核；输入/模型/版本/求解配置可追溯；中断恢复与连续计算在预先声明误差内一致；失败状态不写入 accepted checkpoint | 待建立完整配置驱动链；现有库能力可复用 |
| **M4 工程规模与效率** | M3 提供可重复的完整算例；小规模性能基线可提前建立 | 分析 EOS/flash、装配、通信、线性求解与 I/O 成本，按证据改善结构和资源使用 | 同精度、同工作量报告总耗时、峰值内存、Newton/KSP 迭代、缩步/失败次数及强/弱扩展结果；记录硬件/软件与重复测量波动，优化后保持科学判据 | 待建立端到端与规模基线；局部加速不代替此项 |
| **M5 场景扩展与产品** | M1 的科学证据；按目标依赖 M2/M3/M4；所需参数、数据与许可齐备 | 按实际场景选择多井/schedule、复杂网格、CPA、原油 lumping 或新增物理；接入流动任务、场结果展示和发行 | 每项扩展单独验收；模型/参数/导数/有效域有来源及独立验证；产品实现配置—计算—分析，签名/公证、许可、发布/更新和部署条件分别满足 | 待按场景选定范围，不同时承诺全部扩展 |

针对超临界水研究，先用来源完整的简单体系闭合流动，再决定原油重组分/lumping 的范围；组分参数、输运/热量物性及实验依据是进入下一体系的前置条件。PR76 高级混合规则、CPA 物理三相与 closure、H 相分类、盐库存、电解质/反应、固相/水合物、有限速率传质以及 MPFA/非正交处理，均由实际算例和独立证据驱动，不从现有接口或固定 molality 能力推定已支持。

矩阵随能力变更维护：核对代码与注册测试后更新对应行及必要的基线说明；“已验证”须关联具体提交、运行或可复现参考，注明体系、工况、平台和未覆盖部分。重构前后一致只能支持行为保持，不能单独证明原实现科学正确。里程碑完成须逐条关闭判据；当前代码、测试、文档和历史设计均允许被新证据纠正。按一个可验收问题组织 PR、按小切片推进，保留唯一文档入口和既有测试 ownership；优先跟踪物理/守恒误差、失败率、时间、内存和维护成本，代码行数只作辅助指标。

项目许可证仍由负责人决定。参考成熟软件不构成代码或数据复用授权；资料不足时明确阻塞范围，保留已有可验证能力。
