# MPMC_HNU

面向多相、多组分计算的模块化高性能计算平台，采用可移植的 C++20 计算后端与独立 React/TypeScript 前端。计算核心可以作为库无界面运行。

当前已有独立 AD、PR76/SW92/CPA PT 闪蒸后端、局部灵敏度与部分 physics closure，以及已合并 PR #118 建立的可变相数 fully implicit 非等温多组分流动链：自然变量 1/2/3 相、组分/能量守恒、TPFA Darcy/重力、相渗/毛细压契约、PETSc SNES/GMRES/ASM、相变重启、自适应时间步和 Peaceman 单井控制；同时已有 PT 服务、Web/Electron 前端和 Android 工程应用。安装包候选和工程预览仍不等于正式发行。具体能力与证据范围见下表，测试结果以对应提交的 [GitHub Actions](https://github.com/lihuihnu/MPMC_HNU/actions) 为准。

## 阅读入口与文档分工

本页只维护项目概览、依赖边界和待完成方向；接口、公式、参数、构建命令及验证详情在所属模块维护，历史审计保留在专题文档中。新增进展应更新对应专题与必要的概览条目，避免在多个章节重复粘贴整段状态。

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

| 层次 | 已实现 | 关键边界 |
| --- | --- | --- |
| AD | `Dual<T,N>`、常见初等函数、固定维数 `value_and_jacobian`、分块 `value_and_jacobian_runtime<K>` | 仅依赖标准库；运行期驱动复用固定宽度 Dual。未实现反向模式、通用高阶或稀疏传播。 |
| PR76 | 有序参数、纯/混合/PT 物性、有限 TPD、汽液 PT、maximum-three-phase 编排与 continuation、统一 PT backend | 两相有独立高精度二元/三元参考；max3 除 synthetic structural regression 外，还有 Li–Firoozabadi 六组分文献工程状态的独立数值三相回归，以及 Heringer 2026 来源一致的物理 3→2 边界回归。有限搜索不是全局证明，这些证据也不建立通用相形态分类或宽工况实验精度。 |
| SW92 | corrected-original 物性、fixed-family stability/VLE、Whitson dual-model observables、Xu-style max2、Profile-C authoritative 1/2/3-phase PT | 各算法 profile 保持独立。Profile-C 有 Sample-6 与 n-butane/H2O 实验三相线验证；固定 NaCl molality 不等于盐库存守恒，H 相仍为 `nonaqueous_unclassified`。 |
| CPA | explicit-site-pair 参数、缔合、PT 物性、minimum-Gibbs stability、VLE、max3 与统一 PT backend | 甲醇(2B)/水(4C) 两相已有可追溯物理回归；max3 仍为 synthetic structural validation，缺兼容的物理 VLLE oracle。 |
| 灵敏度与 physics | PR76 内部两相、SW92 Profile-C 固定相集合的局部隐式导数及 closure；两者的 model-neutral component inventory/local Jacobian | 导数绑定已接受的光滑相/根/family 分支；PR76 standalone 单相 closure、CPA flash sensitivity/closure 尚未实现。Inventory 单位为 `mol/m³ fluid`，不是孔隙体积累积项。 |
| Mesh / discretization foundation | 1D/2D/3D structured/unstructured topology/geometry、location-aware fields、DoF/partition、Gmsh 4.1 ASCII、VTU ASCII、声明的 GRDECL 子集、DMPlex/PetscSection/PetscSF 适配；TPFA admissibility/half/static transmissibility、owned connection schedule 与 MPIAIJ symbolic preallocation | `mesh` 不依赖 PETSc/MPI 或具体流动物理；TPFA 几何/传导与符号结构属于 `discretization`。数值守恒装配、mobility/gravity/Darcy flux、时间推进和 PETSc 求解由上层 `flow` / `flow_discretization` / PETSc bridge 消费这些基础能力。详见 [Mesh 模块](modules/mesh/README.md)。 |
| Flow / PETSc | 1/2/3 相 variable-cardinality 自然变量；Backward Euler 组分与能量守恒；TPFA 相势/重力、上风 mobility、组分与焓通量、导热；完整 residual/Jacobian 与 ragged global numbering；PETSc `SNESNEWTONLS + BT -> GMRES + restricted ASM`；post-SNES phase scan、1↔2↔3 rebuild/rebind、自适应 retry/cutback/growth、accepted physical-time/history commit；PR76 真实 transient 与来源完整三相 sour-gas 回归 | TPFA 只在现有 admissibility gate 允许的面上使用；production phase-property / transition 物理证据目前以 PR76 路径最完整。未因此自动获得 MPFA、扩散/弥散、反应、地质力学、裂缝、AMR 或任意 EOS/构成模型的生产覆盖。详见 [Flow 模块](modules/flow/README.md)。 |
| Well | Cartesian/diagonal-K Peaceman WI、fixed-BHP 单/多 completion owner-only source、variable-cardinality rebind；单井 fixed-total-molar-rate 的全局 BHP unknown 与解析四块 Jacobian；minimum-BHP 切换、accepted control state、带滞回的 BHP→rate reactivation，以及 control-arbitration × phase-transition transactional restart | 当前是单井控制基线；不包含 multi-well network / priority、wellbore pressure-drop/heat-loss、surface/phase-rate control、maximum-BHP、completion schedule 或设施模型。详见 [Well 模块](modules/well/README.md)。 |

三个预置 PT 后端均通过 [统一能力与结果契约](modules/flash/pt_flash_backend.md) 暴露 1/2/3 相能力、版本与有序组分身份。预置 backend inventory 本身仍是冻结快照；组分增减、替换、重排不能原地修改该快照。Electron 的免登录 **PR76 Expert** 工作台另行支持用完整的新参数快照创建不可变运行时 PR76 模型，可新增、删除、重排组分并显式给出 `Tc`、`Pc`、偏心因子、摩尔质量、完整 `kij` 与 solver settings，再由同一 C++ 相稳定/最多三相内核求解。这个能力不是通用物性数据库；内建 PR76 甲烷/乙烷/丙烷、SW92 CO₂/淡水、CPA 甲醇/水快照仍只覆盖声明的窄文献体系。

### PR76 验证边界与历史证据

严格 PR76 当前生产路径继续使用 classical vdW one-fluid mixing 与显式常数对称 `kij`。其 pure/mixing/root/`ln(phi)`/TPD/two-phase/max-three-phase/boundary/continuation/publication 链路已经在已合并的 [PR #112](https://github.com/lihuihnu/MPMC_HNU/pull/112) 做系统审计；在该审计覆盖的数学与数值契约内，没有保留的 production implementation defect。完整公式、接口和当前测试入口由 [热力学模块](modules/thermodynamics/README.md) 与 [闪蒸模块](modules/flash/README.md) 维护，根 README 不再复制逐 PR 的历史过程。

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
- **网格与求解链：** `core` 保持可选最小基础层；`numerics` 提供领域无关算法；`mesh` 负责 topology/geometry/field/DoF/I/O/partition，`discretization` 负责 TPFA admissibility/transmissibility 与 PETSc symbolic structure；`flow`/`flow_discretization` 在其上拥有组分/能量守恒、Darcy/热通量、variable-cardinality 自然变量和时间推进，PETSc bridge 负责 MPIAIJ 数值装配与 SNES/KSP/PC 求解。基础模块公共领域接口仍不得泄漏 PETSc/MPI 类型。
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

根 [CMakeLists.txt](CMakeLists.txt) 和 [CMakePresets.json](CMakePresets.json) 当前只配置 AD；热力学、闪蒸、physics、服务和产品采用各自独立入口。请从上述模块文档选择构建命令，避免把根目录构建误认为整个项目的验证。

开发流程以 [AGENTS.md](AGENTS.md) 为唯一规则入口：先审计，做可回退增量，再按受影响依赖选择测试。正式自动化中，普通 Linux 编译、单元测试和常规矩阵默认使用 GitHub 官方 `ubuntu-24.04`，Windows/macOS 使用对应官方托管 runner；当前私有 `mpmc_hnu` 白名单仅覆盖 `flow_discretization_petsc.yml` 的 PETSc/MPI 2-rank 集成/求解 Gate，以及 `cpa_performance_audit.yml` 的经审计长时 paired performance audit。workflow 在私有 runner 上仍必须显式安装/核验固定依赖，不依赖持久机器的偶然环境。纯文档核对内容、相对引用和差异范围，代码变更运行必要增量及受影响下游。解析/数值回归、synthetic fixture、实验验证与性能基准分别报告，不把已有工作流当作已经通过的证据。

| 后续方向 | 所需证据或前置条件 |
| --- | --- |
| 补齐模型能力 | PR76 若引入高级混合规则，先建立公式/导数/参数来源/有效范围契约并补独立宽温区验证；不移植测试 helper 或用受阻三元目标反标。CPA 继续补兼容物理三相参考与灵敏度/closure，PR76 standalone 单相 closure 按独立增量推进。 |
| 扩展流动与求解覆盖 | PR #118 已建立 fully implicit 1/2/3 相自然变量守恒、时间推进与 PETSc SNES/KSP/PC 主链；后续新增 MPFA/非正交修正、扩散/弥散、更多 EOS/构成模型的 production bridge、multi-well/schedule、checkpoint/restart 或场尺度优化时，分别冻结契约并以守恒、收敛和独立物理证据验收。 |
| 细化含水体系 | H 相物理分类、盐库存、电解质/反应、固相/水合物及有限速率传质均需新的模型与证据，不由当前 fixed-molality PT 能力自动获得。 |
| 完成应用发行与结果能力 | 产品签名/公证、项目许可、发布/更新策略、真实部署，以及可视化和结果导出分别验收。工程 APK/MSI/PKG/DEB 不能自动升级为正式发行。 |

项目许可证仍由负责人决定。参考成熟软件不构成代码或数据复用授权；资料不足时明确阻塞范围，保留已有可验证能力。
