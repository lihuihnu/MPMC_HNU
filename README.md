# MPMC_HNU

面向多相、多组分计算的模块化高性能计算平台，采用可移植的 C++ 计算后端与独立 Web 前端。

> **当前状态：AD、PR76 汽液 PT 闪蒸/隐式灵敏度、SW92 corrected-original fixed-family stability / one-family VLE / Whitson dual-model observables、Xu-style asymmetric Gate 3A + Gate 3B.1/3B.2/3B.3 maximum-two-phase、Profile-C fixed-molality authoritative 1/2/3-phase PT publication / fixed-phase-set implicit sensitivity、PR76/SW92 physics thermodynamic closure，以及 model-neutral PT component inventory / local Jacobian。** 已提供独立的 C++20 `mpmc::ad::Dual<T, N>`、常见初等函数、`value_and_jacobian`、分块 `value_and_jacobian_runtime<K>`、独立增量测试入口与官方 runner 工作流；功能边界和使用方法见 [AD 模块说明](modules/ad/README.md)。已新增[有序组分与 PR76 数据契约](modules/thermodynamics/README.md)，并提供经原文核验的[纯组分 a(T)、b 数值核](modules/thermodynamics/pr76_pure.md)及温度导数增量测试；已增加[运行期经典混合参数](modules/thermodynamics/pr76_mixture.md)，区分完整/约化组成导数并验证顺序与变维数；已增加 [PT 候选相性质核](modules/thermodynamics/pr76_phase.md)，给定 p、T、相组成求可用 Z 根与 ln(phi)，含局部隐式导数和退化诊断；已实现 [TPD 相稳定性搜索](modules/flash/README.md) 与 [汽液 PT 相分裂基线](modules/flash/pt_split.md)，给定 p、T、总体摩尔组成 z，联合检查物料守恒与逸度平衡，并对候选两相作共同切平面复核；包含边界回归、TPD 回溯停滞修复和[未确定原因摘要](modules/flash/diagnostic_summary.md)，并增加[模型无关 PT phase-set 表示契约](modules/flash/phase_set.md)。已增加[内部两相收敛解隐式灵敏度](modules/flash/pt_sensitivity.md)，并由 [physics thermodynamic-closure consumption contract](modules/physics/README.md) 将已接受状态映射为相摩尔分率、组成、Z、摩尔密度及可用的局部线性化，同时把 residual 可用性与 Newton 基点可用性分开。SW92 已依据 Søreide–Whitson 1992 原文及作者勘误实现 [corrected-original thermodynamics kernel](modules/thermodynamics/sw92.md)，并接入[同一 phase-family 的 PT stability](modules/flash/sw92_stability.md)和[one-family PT VLE](modules/flash/sw92_family_vle.md)；在此基础上已实现 [Whitson dual-model observable orchestration](modules/flash/sw92_dual_model.md)，从同一个 ordered snapshot 分别执行完整 AQ 与 NA family run，并仅在已审计的二元情形按 water-rich / water-poor 提取兼容性 observable。该双模型 Profile **不是** AQ/NA 联合热力学 phase set：两次 run 的 phase fraction、material balance、common tangent 与 Gibbs 判据保持独立，不得拼接成共同守恒结果。另一路 [`SW92-equilibrium/xu-asymmetric-gibbs/v1`](modules/flash/sw92_xu_asymmetric_audit.md) 已按 Gate 分阶段实现：Gate 3A [asymmetric stability foundation](modules/flash/sw92_asymmetric_stability.md) 建立 lower-envelope feed reference/common tangent 和 AQ/NA finite TPD；Gate 3B.1 [family-aware fixed-pair primitive](modules/flash/sw92_asymmetric_fixed_pair.md) 对显式 `(F0,F1)`（允许 AQ+AQ / AQ+NA / NA+NA）求共同 EOS-component material balance/common reduced chemical potentials 并逐相检查 lower-envelope family assignment；Gate 3B.2 [witness orchestration](modules/flash/sw92_asymmetric_orchestration.md) 从完整 family-tagged negative-witness 集合生成 family-neutral plan、mandatory same-family alternatives、执行 fixed-pair、做 slot-swap/equivalent-candidate 去重并进行 guarded Gibbs selection；Gate 3B.3 [maximum-two-phase final acceptance](modules/flash/sw92_asymmetric_max2.md) 进一步重算 acceptance-critical pair evidence、检查 pair-vs-lower-feed Gibbs、对 selected-pair common tangent 同时执行最终 AQ/NA finite stability，并仅在全部门槛通过时由 family-aware `accepted_phase_set()` 发布一相或两相结果。第三条 Profile-C `SW92-equilibrium/phase-assigned-aq-na-joint/v1` 已完成 C1 `W(AQ)+H(NA)` joint candidate、C2a1 H-side NA additional-phase witness、C2a2 morphology/family/root/topology separation、C2b.1 unordered `W(AQ)+H0(NA)+H1(NA)` 三相 candidate、C2b.2 W-present H-multiplicity review、相消失邻接 topology fresh re-solve，以及 [authoritative `PtPhaseSetResult` publication adapter](modules/flash/sw92_profile_c_phase_set.md)。生产入口 `solve_sw92_profile_c_pt_phase_set(...)` 只运行一次 boundary-aware Profile-C PT driver，再纯投影已拥有的收敛 phase property/provenance；在声明的 finite-search/topology contract 下可 authoritative 发布 1/2/3 相。对一个已接受 authoritative phase set，现已提供 [fixed-phase-set implicit sensitivity](modules/flash/sw92_profile_c_sensitivity.md)：SW92 pure/mixing/phase-property arithmetic 支持 AD，局部 IFT 在固定 phase count、AQ/NA family、selected root 与 representation slots 下给出 `beta/x/Z/c` 对 `(p,T,z_reduced)` 的 Jacobian；[SW92 physics closure](modules/physics/README.md) 将该 Jacobian 与 primal 原子地发布到同一 variable-cardinality physics snapshot。在这两类 physics closure 之上，[PT component inventory](modules/physics/component_inventory.md) 已统一消费 fixed-VLE 与 variable-cardinality payload，以 `v_bar=sum(beta/c)` 构造每总流体体积的 `c_mix` 与 ordered component molar inventory `a_i`，并解析组合 `d beta/dq`、`d x/dq`、`d c/dq` 发布 local Jacobian；0D fixed-fluid-volume isothermal residual 仅作为 integration regression，用 fresh PR76/SW92 re-solve 交叉核验，不是 production discretization API。物理证据包括 Mortezazadeh–Rasaei Sample-6 多组分三相回归以及 Reamer et al. n-butane/H2O 实验三相共存线。**这仍不是数学全局稳定性证明：`global_stability_proven=false` 保持；H0/H1 继续为 `nonaqueous_unclassified`，不以 Z/root 冒充 LV/LL morphology；当前 NaCl molality 是给定模型参数，不是 salt-inventory conservation；fixed-phase-set Jacobian 不跨相出现/消失、family/root switching、临界/近重根或 topology 切换外推；component inventory 是 total-fluid-volume quantity，不是 pore-volume accumulation，也不引入 porosity/saturation。** CPA、PR76 standalone 单相 physics closure、pore-volume accumulation、production conservation residual/flux、网格、离散、time stepping、global Jacobian 与全局求解器尚未实现。已有软件与模型数值回归不等于普适实验验证或性能达标；测试结果以具体提交的 GitHub Actions 日志为准。开发约束见 [AGENTS.md](AGENTS.md)。

## 1. 项目目标与基本原则

以通用、高效的自动微分（Automatic Differentiation，AD）模块为起点，逐步建立热力学、两相与三相闪蒸、网格、流动、多相多组分物理模型、数值算法、可视化和结果导出能力。

- **可移植、可组合：** 计算核心不依赖前端、网络服务、特定操作系统或特定 HPC 供应商；支持无界面、单进程使用，再按需求增加并行后端。
- **科学正确、证据可追溯：** 公式、参数、验证数据和性能结论必须有来源；区分模型假设、数值近似与实验事实，不伪造数据或测试结果。
- **审计在前、小步实现：** 每个增量先审查全局依赖和局部算法，再编写、验证与提交。允许有证据、有迁移方案地推翻并重构既有设计，不把早期选择视为永久约束。

## 2. 初始技术基线

下表区分当前 AD 基线与后续设计选择；未说明已实现的项均为规划，不代表已安装或通过跨平台验证。具体依赖版本在首次引入时核查官方支持情况并锁定。

| 层次 | 初始选择 | 边界与理由 |
| --- | --- | --- |
| 计算核心 | C++20，标准库优先 | AD 算术和初等函数已实现，且只依赖标准库；已增加 PR76 纯组分、经典混合系数、PT 候选相性质核、C++20/double 汽液 PT 闪蒸、内部两相收敛解灵敏度、SW92 corrected-original thermodynamics、fixed-family stability / one-family VLE、Whitson dual-model observable orchestration、Xu-style Gate 3A + Gate 3B.1/3B.2/3B.3 maximum-two-phase，以及 Profile-C C1/C2a1/C2a2/C2b.1/C2b.2、boundary-aware disappearance-neighbor re-solve、authoritative 1/2/3-phase publication、fixed-phase-set implicit sensitivity、physics closure 与 model-neutral PT component inventory/local Jacobian。PR76 两相和 SW92 Profile-C 1/2/3 相均已有最小 thermodynamic-closure consumption；两种 generic closure 均可进一步映射到每总流体体积 component inventory。Profile-C H morphology、salt inventory、跨 phase/topology/root/family 边界导数、pore-volume accumulation、production conservation residual/flux、网格/离散与全局求解仍待独立 gate。其余计算能力按路线开发，不将 GPU、MPI 或专有指令集作为基础依赖。 |
| 构建与测试入口 | CMake 3.21+、CMake Presets、CTest | 已提供独立 `mpmc::ad`、`mpmc::thermodynamics`、`mpmc::flash`、opt-in `mpmc::flash_sensitivity` 与 `mpmc::physics` 目标及对应增量测试入口；PR76 汽液 PT、灵敏度、physics closure、SW92 thermodynamics、fixed-family stability / VLE、dual-model compatibility observables、asymmetric Gate 3A/3B.1/3B.2/3B.3、Profile-C C1/C2a1/C2b.1/C2b.2、boundary re-solve、authoritative phase-set publication、fixed-phase-set sensitivity、SW92 physics closure 与 PT component inventory 均有独立或增量验证。Profile-C publication/sensitivity/physics 已覆盖 1/2/3 相、physical Sample-6、component permutation、H-slot symmetry、AD/IFT derivative checks、boundary/unavailable semantics 与 provenance guard；component-inventory suite 另覆盖 fixed-VLE/variable-cardinality consumption、phase-slot invariance、derivative-unavailable propagation，以及 PR76/SW92 0D fixed-volume residual fresh-re-solve Jacobian cross-check。[E1] |
| 前端 | React + TypeScript + Vite 单页应用 | 前后端独立开发与部署；不为计算平台默认引入 SSR 或另一套服务端业务逻辑。[E2] |
| 服务通信 | Protocol Buffers + gRPC；浏览器通过 gRPC-Web 适配层接入 | 契约先行、消息版本化；传输对象与计算核心类型分离。[E3][E4] |
| 大体量结果传输 | 独立 HTTP 二进制分块通道 | 网格、场变量不默认转为 JSON 数组或 Base64；支持有界缓冲、按需读取和校验。具体格式在结果模块开发时确定。 |
| 科学可视化 | 优先评估 vtk.js | 仅作为前端可视化适配层，不进入计算核心；大数据采用按需加载与分辨率分级。[E5] |
| 并行与求解扩展 | 可替换执行策略与 `SolverBackend` | 先提供串行基线；OpenMP、MPI、PETSc 或 GPU 后端须经独立需求、许可与收益审计后引入。 |

**可移植性目标：** Linux、Windows、macOS，覆盖 GCC、Clang、MSVC；处理器架构和编译器版本以实际验证矩阵为准，未测试的平台不得标记为支持。不依赖绝对路径、大小写不敏感文件系统或仅在单一平台可用的系统调用。通用发行构建不得默认开启 `-march=native` 或不受控的 `fast-math`。

**通信约束：** 浏览器不能被当作原生 gRPC 客户端。当前官方 gRPC-Web 实现支持 unary 和特定模式的服务端流，不支持客户端流或双向流；后续实现须复核所选版本的能力，不承诺不存在的双向通信。[E4] 首期服务契约以任务提交、查询、取消、结果索引和进度事件为目标；进度流采用受支持的适配模式，大数组走独立二进制通道。出现真实双向交互需求时，再审计 WebSocket 等补充方案，不预先建设多套通信栈。

任务协议应明确请求 ID、任务 ID、版本、超时、取消、重试幂等性、背压、错误码和能力查询。二进制数据须携带类型、维度、布局、字节序、单位和校验信息；跨机器不得直接传输 C++ 对象内存。服务端必须独立校验输入，设置资源限制，并在远程部署时提供认证、授权和加密传输。单次 AD 运算或网格单元计算不得变成远程调用。

## 3. 模块架构与依赖方向

下表描述目标职责；**当前 AD 已有计算代码，`modules/thermodynamics` 已有通用数据/来源基础、完整 PR76 基线以及 SW92 corrected-original double/AD-capable phase-family PT 候选物性核；`modules/flash` 已有有限多初值相稳定性搜索、PR76 汽液 PT 相分裂、最终共同切平面复核、模型无关 phase-set 表示、未确定诊断、已接受内部两相收敛解隐式灵敏度，以及 SW92 fixed-family stability / one-family VLE / Whitson dual-model observable orchestration / Xu-style Gate 3A + Gate 3B.1/3B.2/3B.3 max2 / Profile-C authoritative 1/2/3-phase PT publication与 fixed-phase-set implicit sensitivity；`modules/physics` 已有 PR76 fixed-VLE 与 SW92 Profile-C variable-cardinality thermodynamic-closure consumption contract，以及同时消费两种 generic closure 的 model-neutral PT component inventory/local Jacobian；其余流动/离散模块尚未创建**。当前 SW92 flash 路径保留显式 AQ/NA family 与 cubic-root branch 分离：fixed-family stability / VLE 是单模型内部自洽 primitive；Whitson dual-model Profile 保留 AQ/NA 两个完整 run 并提取兼容性 observable；Xu-style 路径则使用 lower-envelope family model，Gate 3A 做 common-tangent stability，Gate 3B.1 做显式 family-pair joint equations，Gate 3B.2 做 family-neutral candidate orchestration，Gate 3B.3 做 pair-vs-feed Gibbs 与 selected-pair final AQ/NA stability 后才允许 family-aware `accepted_phase_set()`。Profile-C 固定 physical `W->AQ`、`H->NA`：C1 求 W+H candidate，C2a1 在 C1 common tangent 上搜索 additional NA witness，C2b.1 求 symmetric/unordered W+H0+H1 joint candidate，C2b.2 做 W-present H-multiplicity final review，boundary-aware adapter 对相消失必须 fresh re-solve 邻接 topology；只有完整 finite topology/provenance chain 闭合后，`solve_sw92_profile_c_pt_phase_set(...)` 才把已拥有结果纯投影到 generic `PtPhaseSetResult`。C2a2 仍规定 H morphology 不由 family 或 cubic root 决定，因此 authoritative phase set 只发布 `aqueous` 或 `nonaqueous_unclassified`，不发布未经验证的 `LV/LL`。Xu `max2` 与 Profile-C accepted finite-search result 都不声明数学全局稳定或盐守恒。当前 physics 可消费已接受 PR76 两相状态，以及 SW92 Profile-C authoritative 1/2/3-phase 状态；在对应 fixed smooth branch 的 flash sensitivity 可用时，同一 snapshot 还拥有 `beta/x/Z/c` 对 reduced PT-feed coordinates 的局部线性化，并以 `can_seed_newton()` 与 residual availability 分离。`build_pt_component_inventory(...)` 再将 `beta/x/c` phase-resolved state 映射为 `c_mix` 和 ordered component molar inventory `a_i`，并在源 derivative 可用时发布同坐标 local Jacobian；该量按 **total fluid volume** 定义，不是 pore-volume accumulation，也不改变全局 Newton admissibility。该能力仍不等于 production conservation residual、flux、饱和度、网格、离散或全局求解器。只在对应增量需要时创建文件与构建目标，不预生成空模块、占位实现或插件框架。

| 模块 | 职责 |
| --- | --- |
| `core` | 单位、基础数据约定、错误与诊断；保持最小依赖。 |
| `ad` | 通用标量 AD、导数传播与导数访问；仅依赖标准库，不依赖 `core` 或其他项目模块。 |
| `numerics` | 通用求根、线性代数接口、非线性算法、容差与收敛诊断。 |
| `thermodynamics` | 组分与参数、热力学模型、物性及其导数、适用范围与模型能力。 |
| `flash` | 相稳定性、相数判定、两相及三相平衡与结果诊断。 |
| `mesh` | 网格数据结构、创建、拓扑与几何有效性检查。 |
| `physics` | 守恒方程、本构关系，以及水独立或水参与互溶的物理模型。 |
| `discretization` | 在物理模型与网格上进行离散、残差和 Jacobian 组装。 |
| `solvers` | 时间推进、线性/非线性求解协调与可选 `SolverBackend` 适配。 |
| `io` | 网格导入、结果与检查点读写、格式适配、来源元数据。 |
| `runtime` / `api` | C++ 应用 API、能力目录、参数校验、Compute Worker、任务与通信适配。 |
| `frontend` | 建模交互、任务控制、结果展示与可视化，不重复实现科学模型。 |

依赖必须无环。`ad` 是仅依赖 C++ 标准库的独立模块，不依赖 `core` 或其他项目模块；通用数值模块可以按需要使用 AD，但不得反向依赖热力学。热力学可以使用 AD 和通用数值算法，但不得调用闪蒸；闪蒸通过热力学能力接口求解平衡。物理与离散层使用上述能力，应用层负责协调；网格和格式适配不应污染 AD 或热力学核心。

计算库必须可以被 C++ 程序直接调用，而不启动服务或浏览器。第三方求解器、网格库和通信库通过适配层接入；公共领域接口不得暴露 PETSc/MPI、Protobuf、HTTP 或 UI 类型。模块化首先通过清晰的 C++ 接口和构建目标实现，不等于微服务化，也不承诺不同编译器之间的 C++ 二进制 ABI 兼容。

## 4. 第一阶段：通用自动微分

先建立可独立构建、可独立验证的 AD 模块，再将其用于物性与平衡残差。核心函数按标量类型编写，使普通浮点数与 AD 标量能够复用同一数学表达式；不得在计算中无意提取普通数值而丢弃导数。

首个代码增量已提供固定维数一阶前向 `Dual<T, N>`，只处理值/导数语义、种子与基本算术。现已在独立数学头文件中增加常见初等函数；逐函数定义域、异常与增量测试入口见 [初等函数契约](modules/ad/math.md)。使用与功能边界见 [AD 模块说明](modules/ad/README.md)。已增加固定维数 [函数值与 Jacobian 提取接口](modules/ad/differentiate.md)；后续按实际需求扩展其他数学函数和求导策略。反向模式（VJP）、稀疏传播、高阶导数与表达式优化作为后续能力，每项均需需求与代价审计，不在第一次实现中一次性承诺完成。[E6]

明确可微区间、定义域、分支、别名、零/负参数、溢出、对象生命周期和线程安全语义。AD 并不使不可微函数变得可微，也不保证模型本身正确。对 EOS 根、CPA 缔合位点方程和闪蒸等隐式问题，应区分残差导数、迭代过程导数与收敛解导数；在相应阶段审计隐式求导及内外层容差，不能仅因数值求解收敛就宣称导数正确。

验证以可推导的解析函数为基础，辅以有适用条件的差分步长扫描、复步法或独立 AD 实现交叉核验。有限差分不得充当生产 AD 实现，也不得是唯一正确性证据。性能比较必须固定问题、精度、构建选项和环境，同时报告时间、内存与误差；目前没有任何性能达标结论。

运行期输入/输出维数已由独立的 [分块 Jacobian 驱动](modules/ad/runtime_differentiate.md) 提供：复用 `Dual<T,K>`，每次传播固定 K 个方向，保留原固定接口；这不等于实现动态标量梯度、变组分业务管理或实时前后端链路。

## 5. 热力学模型与科学边界

当前已实现的[PR76 数据契约](modules/thermodynamics/README.md)按稳定组分 ID 绑定参数并形成运行期有序快照；公共来源、单位与身份层可由其他热力学模型复用，不将 PR 的常数对称 kij 规则强加给它们。PR76 已增加纯组分、经典混合 a/b 与 [PT 候选相性质](modules/thermodynamics/pr76_phase.md)，支持 Z、ln(phi) 和简单根分支的局部导数；根数不等于相数。在此基础上已实现 [TPD 相稳定性搜索](modules/flash/README.md)、[PR76 汽液 PT 基线](modules/flash/pt_split.md)与[已接受内部两相收敛解隐式灵敏度](modules/flash/pt_sensitivity.md)，联合验收守恒、逸度与最终相集合的有限稳定性搜索结果，并为固定相集合/根分支提供局部一阶导数；不提供全局稳定性或跨相边界光滑性认证，也不等于已实现三相或含水联合闪蒸。首个 [physics thermodynamic closure](modules/physics/README.md) 消费上述 PR76 能力并增加 Z/相摩尔密度及其链式导数，不新增经验热力学模型。

[SW92/corrected-original thermodynamics kernel](modules/thermodynamics/sw92.md) 已依据 Søreide–Whitson 1992 原论文及其作者勘误实现：水/brine alpha、AQ water-pair correlations、显式 NA water-pair数据/H2S Eq.(17)、family-specific non-water BIP、经典 PR mixing、Z 根与 ln(phi) 候选物性。`SwPhaseFamily::{aqueous,nonaqueous}` 与 cubic root index 明确分离；NaCl 在当前 profile 中是外部 fixed molality，不作为守恒的 Na+/Cl- EOS 组分。该热力学核保留 double primal 路径，同时 pure/mixing/phase-property arithmetic 已支持 fixed-phase-set sensitivity 所需的 AD 数值类型；生产 flash 的离散 phase/family/root/topology 选择本身不通过 AD 求导。它已接入 [fixed-family PT stability](modules/flash/sw92_stability.md) 与 [one-family PT VLE](modules/flash/sw92_family_vle.md)：一次 run 固定同一 family/molality，stability 在该 family 内选择 mechanically admissible minimum-Gibbs root，split 继续把 lowest/highest-Z 仅作为请求的数值 candidate role，并复用 generic TPD / Rachford-Rice / material-balance / fugacity / final common-tangent 逻辑。[equilibrium-algorithm contract](modules/flash/sw92_equilibrium_algorithm.md) 区分模型身份与算法身份；[Whitson dual-model observable Profile](modules/flash/sw92_dual_model.md) 从同一 ordered snapshot 独立执行 AQ 与 NA 两个完整 run，只在二元 compatibility scope 内提取 water-rich AQ 与 water-poor NA observable，并显式输出 `cross_model_equilibrium_ratio`，不构成共同 material balance。经 [Xu-style asymmetric-Gibbs 科学审计](modules/flash/sw92_xu_asymmetric_audit.md) 和 [Gate 3B joint-split 设计审计](modules/flash/sw92_gate3b_joint_split_audit.md)，另一条独立算法路径已实现 Gate 3A [asymmetric stability](modules/flash/sw92_asymmetric_stability.md)、Gate 3B.1 [fixed pair](modules/flash/sw92_asymmetric_fixed_pair.md)、Gate 3B.2 [candidate orchestration](modules/flash/sw92_asymmetric_orchestration.md) 与 Gate 3B.3 [max2 final acceptance](modules/flash/sw92_asymmetric_max2.md)。Gate 3B.3 从 retained family-aware phases 重算共同化学势/物料衡算/Gibbs/common tangent 关键证据，比较 pair 与 lower-feed Gibbs，并把两个 candidate compositions 都加入 AQ/NA final starts；只有两边 final finite search 均 `no_instability_found` 才允许 `accepted_phase_set()` 发布。该 Xu max2 能力仍不包含 global stability certification、三相或 salt inventory；现有 SW92 导数/physics handoff 专属于下述 Profile-C authoritative fixed-phase-set contract，不能反向解释为 Xu max2 的 derivative contract。

第三条独立算法路径 `SW92-equilibrium/phase-assigned-aq-na-joint/v1` 保留 SW92 physical phase-assignment 语义而不使用 Xu lower envelope。Profile-C C1 [`sw92_phase_assigned_joint.hpp`](modules/flash/include/mpmc/flash/sw92_phase_assigned_joint.hpp) 实现固定 `W(AQ)+H(NA)` candidate；C2a1 [`sw92_phase_assigned_h_side_witness.hpp`](modules/flash/include/mpmc/flash/sw92_phase_assigned_h_side_witness.hpp) 从重新核验的 C1 common tangent 做 NA-family finite additional-phase witness search；C2a2 [hydrocarbon-role audit](modules/flash/sw92_phase_assigned_c2a2_hydrocarbon_role_audit.md) 确认 `physical hydrocarbon morphology`、`SwPhaseFamily::nonaqueous`、cubic root/Z 与 phase instance 必须分离；C2b.1 [`sw92_phase_assigned_three_phase.hpp`](modules/flash/include/mpmc/flash/sw92_phase_assigned_three_phase.hpp) 求 unordered `W(AQ)+H0(NA)+H1(NA)` candidate，使用二维 generalized Rachford-Rice 与两组 W-H reduced chemical-potential residual；C2b.2 对 W-present H multiplicity 做最终有限 review。相消失不能因小 fraction 直接删相：[`sw92_phase_assigned_boundary.hpp`](modules/flash/include/mpmc/flash/sw92_phase_assigned_boundary.hpp) 对请求降相的邻接 topology fresh re-solve，并重做其专属 stability/topology review。最终 [`sw92_profile_c_phase_set.hpp`](modules/flash/include/mpmc/flash/sw92_profile_c_phase_set.hpp) 以该 boundary-aware driver 为唯一生产求解来源，通过 [`solve_sw92_profile_c_pt_phase_set(...)`](modules/flash/sw92_profile_c_phase_set.md) 将已有收敛 phase fraction/composition/activity/Z 与完整 provenance **纯投影**为通用 `PtPhaseSetResult`；它不重新调用 EOS、不重求根、不重跑 TPD/守恒/平衡，也不修改容差。authoritative `accepted` 只表示当前 fixed-molality Profile-C finite-search/topology contract 已闭合，可发布 1/2/3 相；`global_stability_proven=false` 仍固定。H0/H1 无论槽位如何都只为 `nonaqueous_unclassified`/NA，relative Z 仅用于 deterministic representation，不是 LV/LL 分类器。对该 authoritative phase set，可通过 [`differentiate_sw92_profile_c_phase_set(...)`](modules/flash/sw92_profile_c_sensitivity.md) 获取 `q=(p,T,z_0,...,z_{N-2})` 下的 fixed-phase-set implicit Jacobian：输出 `d beta/dq`、`d x/dq`、`d Z/dq`、`d c/dq`，以 AD 构造局部平衡方程 `A=∂F/∂u`、`B=∂F/∂q` 并解 `du/dq=-A^{-1}B`。内部使用 log-ratio/softmax chart 保持 simplex/positivity 并改善 trace-component conditioning；不对 topology、TPD、SSI、family/root selection、相出现/消失或 H morphology 求导。随后 [`build_sw92_profile_c_thermodynamic_closure(...)`](modules/physics/README.md) 先重现 authoritative primal，再原子消费该 sensitivity；导数成功时 physics snapshot 同时提供 primal 与 variable-cardinality linearization，derivative-only failure 则保留 primal 但 `can_seed_newton()==false`，`solution_not_accepted` 一致性失败则撤销整个 snapshot。物理验证已包括 Mortezazadeh–Rasaei 2017 Sample-6 多组分 W+H0+H1 interior 状态和 Reamer et al. 1944 n-butane/H2O 实验三相线；这些证据验证声明范围内的模型/算法，不扩大为普适全局证明。

在上述 thermodynamic closure 之上，[`build_pt_component_inventory(...)`](modules/physics/component_inventory.md) 以 phase-resolved 形式计算 `v_bar=sum(beta_alpha/c_alpha)`、`c_mix=1/v_bar`、`s_i=sum(beta_alpha*x_alpha,i)` 与 `a_i=c_mix*s_i`，并要求 `s_i≈z_i` 及等价 feed-form identity。若源 local derivative 可用，则解析组合 `d beta/dq`、`d x/dq`、`d c/dq` 得到 `d c_mix/dq` 与 `d a_i/dq`，并逐列检查 differentiated phase/feed identity；源 derivative unavailable 时仍可保留 valid inventory primal，但不伪造 Jacobian。该层同时支持现有 PR76 fixed-VLE 与 SW92 variable-cardinality closure；量纲是 `mol/m^3 fluid`，不包含 pore volume、porosity 或 saturation。0D fixed-fluid-volume isothermal residual `R_i=V*a_i-N_i_target` 只作为 integration test，Jacobian 已用 fresh symmetric PR76/SW92 re-solves 交叉验证；它不是 production conservation/discretization API。

模型选择通过能力目录与配置完成，不把 EOS 公式硬编码进闪蒸算法。每个模型应报告模型 ID/版本、可用组分、允许相态、所需参数、混合规则、导数能力、适用范围及资料来源；不支持的组合显式拒绝，不静默回退到另一模型。

| 模型 | 当前状态与实施约束 |
| --- | --- |
| PR（Peng–Robinson） | PR76 已作为首个完整基线实现至汽液 PT、内部两相灵敏度和最小 physics closure；其 accepted two-phase closure 现也可进入 model-neutral PT component inventory/local Jacobian。仍须明确版本、alpha、混合规则、BIP 与任何体积修正，不能把普通 PR 等同于已验证的含水互溶模型。[S1] |
| SW（Søreide–Whitson） | 已实现 `SW92/corrected-original/PR76-base/NaCl-molality` double/AD-capable thermodynamics、fixed-family stability / one-family VLE、`SW92-equilibrium/whitson-dual-model-observables/v1` compatibility orchestration、`SW92-equilibrium/xu-asymmetric-gibbs/v1` Gate 3A + Gate 3B.1/3B.2/3B.3 maximum-two-phase，以及 Profile-C C1/C2a1/C2a2/C2b.1/C2b.2、fresh disappearance-neighbor routing、authoritative generic `PtPhaseSetResult` publication、fixed-phase-set implicit sensitivity 与 variable-cardinality physics thermodynamic closure；authoritative closure 现也可进入 model-neutral PT component inventory/local Jacobian。Profile-C 可在当前 finite-search/topology contract 下发布 1/2/3 相，并在固定 accepted phase set / family / root / representation slots 的光滑局部域内发布 `beta/x/Z/c` 对 reduced PT-feed coordinates 的 Jacobian；已由 Sample-6、n-butane/H2O 实验三相线、component permutation、H-slot symmetry、AD/IFT checks、physics handoff 与 component-inventory regressions 验证。但 `global_stability_proven=false`，H morphology 仍是 `nonaqueous_unclassified`，fixed molality 不代表 salt-inventory conservation，Jacobian 不跨离散相/根/family/topology 边界外推。[S2][S6] |
| CPA（Cubic-Plus-Association） | 尚未实现；实现时须明确立方项版本、缔合位点方案、缔合参数、交叉缔合及混合规则；缔合子问题须同时满足残差与导数要求。[S3] |

原始论文、官方物性资料或有许可的数据集应在实现前核对。仅有摘要、搜索结果或二手转述时，不声称已核验完整公式和参数表。缺少全文、勘误、组分数据、交互参数或商业软件访问时，应明确缺口并向项目负责人索取，不猜测、不伪造、不绕过许可。

参数必须记录组分身份、来源、单位、温压范围、版本及转换过程。未知参数不得静默设为零；经批准采用零相互作用参数等假设时，也必须记录假设与适用边界。首期不默认包含电解质反应、固相、水合物或非平衡传质，这些需单独建模与验证。

## 6. 两相、三相闪蒸与水的处理

### 6.1 当前可用：PR76 汽液 PT 与 SW92 fixed-family / dual-model compatibility / Xu-style max2 / Profile-C authoritative 1/2/3-phase PT

现有 `solve_pr76_pt_vle` 输入压力 p（Pa）、温度 T（K）和总体摩尔组成 z，依次执行进料 TPD、失稳试探组成的守恒初始化、汽液相分裂和最终共同切平面复核。使用 CMake 目标 `mpmc::flash`，无需前端或网络服务；完整接口、选项、数值门槛与构建入口见 [汽液 PT 契约](modules/flash/pt_split.md)。对最终接受的内部两相状态，可通过 opt-in `mpmc::flash_sensitivity` 和 [隐式灵敏度契约](modules/flash/pt_sensitivity.md) 获取 `beta/x/y/logK` 对 `(p,T,z_reduced)` 的局部一阶导数；该导数不跨泡露点、相消失、临界/根切换或相集合变化外推。

以下为调用片段：调用方先建立有来源、单位和模型约定记录的 `PrParameterSet parameters`；提供有限且正的 `p_pa`、`t_kelvin`，以及有限、非负且与参数快照顺序和数量一致的摩尔组成容器 `z`（如 `std::vector<double>`）。参数由调用方提供，接口不内置流体数据库；输入只接受组成和在契约规定的舍入范围内偏离 1，不将负组成裁剪为零或填补缺失组分。输入／参数错误仍可能抛出异常，不能假定所有失败都封装在返回状态中。

```cpp
#include <mpmc/flash/pr76_split.hpp>

const auto phase =
    mpmc::thermodynamics::Pr76Phase<double>::from_parameters(parameters);
mpmc::flash::Pr76VleEvaluator evaluator(phase);
const auto result =
    mpmc::flash::solve_pr76_pt_vle(p_pa, t_kelvin, z, evaluator);
const auto& solution = result.solution;
// 先按 solution.status 分流；candidate() 存在不等于最终稳定性复核通过。
```

| `solution.status` | 调用方应采用的解释 |
| --- | --- |
| `single_phase_no_instability_found` | 初始有限搜索未检出失稳；单相候选组成取记录的 `solution.initial_stability.feed`，不强行创建第二相。 |
| `two_phase_no_instability_found` | 两相通过联合方程、相区分及 Gibbs 检查，最终共同切平面有限搜索未检出进一步失稳；可按本基线的两相结果使用。 |
| `phase_set_unstable` | 两相方程已收敛，但最终复核找到进一步失稳证据；保留候选和证据，不能作为已接受的两相平衡结果。 |
| `indeterminate` | 初始 TPD、分裂或最终复核未完成可靠判定；保留阶段和终止诊断，不转换为单相或两相成功。 |

仅在 `two_phase_no_instability_found` 下按已接受两相处理；将非空 `solution.candidate()` 指针记为 `pair` 后，汽相摩尔分率为 `pair->fractions.vapor_fraction`，液相、汽相摩尔组成分别为 `pair->fractions.liquid` 和 `pair->fractions.vapor`。其他状态也可能保留已收敛候选，须结合 `solution.final_stability` 判读。`equations_converged()` 只表示存在被选中的方程收敛候选，不等于稳定性通过；即使底层 `attempt.point` 满足方程，若小相低于数值接受门槛，也可能没有顶层候选。

所有结果的 `global_stability_proven` 均为 `false`；“未检出”不是全局认证。`solution.diagnostic` 的分阶段摘要、`attempts` 与初始／最终 TPD 的 trial 记录用法见 [未确定原因摘要](modules/flash/diagnostic_summary.md)。机器分流应读取已有枚举、候选与配额字段，不匹配整条诊断字符串；`phase_disappearance` 不等于物理上严格不存在该相。结果拥有数据，候选指针随结果对象存活；evaluator 仅顺序复用，不能并发共用同一实例。

SW92 当前提供三条不同用途的平衡路线。`solve_sw92_pt_family_vle` 只处理一个显式 `SwPhaseFamily`，结果算法身份为 `SW92-equilibrium/fixed-family-vle-primitive/v1`。[`solve_sw92_whitson_dual_model_observables`](modules/flash/sw92_dual_model.md) 的算法身份为 `SW92-equilibrium/whitson-dual-model-observables/v1`，从同一个 snapshot 分别执行 AQ 和 NA 两次完整 one-family run；自动 target-phase extraction 当前仅限二元 water+gas compatibility case，任何 `cross_model_equilibrium_ratio` 都是两个独立模型结果之间的 observable，不能解释为共同 Rachford-Rice `K` 或共同 phase fraction。Xu-style 路线的完整 maximum-two-phase 入口为 [`solve_sw92_xu_asymmetric_max2`](modules/flash/sw92_asymmetric_max2.md)，算法身份 `SW92-equilibrium/xu-asymmetric-gibbs/max2-logK-SSI-RR/common-tangent/v1`：它执行 Gate 3A 初始 lower-envelope stability、Gate 3B.1 fixed-family-pair equations、Gate 3B.2 family-neutral witness/candidate orchestration，并在 Gate 3B.3 对 selected pair 重算共同化学势/物料衡算/Gibbs/common tangent、检查 pair-vs-feed Gibbs，把两个 candidate compositions 同时加入 AQ 与 NA final search；仅 `single_phase_no_instability_found` 或 `two_phase_no_instability_found` 时 `accepted_phase_set()` 非空。`phase_set_unstable`、`pair_gibbs_above_feed` 或 `indeterminate` 可以保留 candidate/diagnostics，但不得作为 accepted phase set 使用。该结果最多两相且 stability 为 finite multistart，仍然不是全局证明。

Profile-C 的 production publication 入口为 [`solve_sw92_profile_c_pt_phase_set(...)`](modules/flash/sw92_profile_c_phase_set.md)。它固定 physical-role/family mapping：`W->AQ`、每个 `H->NA`；内部执行 no-W topology、C1、C2a1、必要的 C2b.1/C2b.2 以及 disappearance-neighbor fresh re-solve。任何降相都必须由邻接 topology 重新求解并重新通过其专属有限 stability/topology review，不能因 small fraction 或 route status 直接删相。只有 topology/provenance 完整闭合后才把已有 phase property 投影到通用 `PtPhaseSetResult`，支持 accepted 1/2/3 phase。模型无关 payload 不增加 aqueous/liquid/vapor 标签；SW92 wrapper 另存 phase family 与 Profile-C role，H0/H1 始终为 `nonaqueous_unclassified`。`accepted` 是当前 finite-search contract 下的 authoritative publication，不是 mathematical global-stability proof。对已接受结果，`differentiate_sw92_profile_c_phase_set(...)` 提供 fixed phase set / fixed family / fixed root / fixed slot 的 `beta/x/Z/c` 对 `q=(p,T,z_reduced)` 局部 Jacobian；`build_sw92_profile_c_thermodynamic_closure(...)` 将该 Jacobian 与 authoritative primal 原子发布，只有 linearization status 为 available 且 payload 与 primal/provenance/shape 一致时 `can_seed_newton()==true`。phase-boundary、ill-conditioning、unsupported support 或 derivative property/arithmetic failure 均不伪造 Jacobian；`solution_not_accepted` 一致性失败不保留可用 primal。

对于 PR76 fixed-VLE closure 或 SW92 Profile-C variable-cardinality closure，可进一步调用 [`build_pt_component_inventory(...)`](modules/physics/component_inventory.md)。该 model-neutral 层按 phase-resolved `beta/x/c` 构造 `c_mix` 与 ordered component inventory `a_i [mol/m^3 fluid]`，并在 source linearization available 时发布 `d c_mix/dq` 与 `d a_i/dq`；它同时用 feed form 和 differentiated material-balance identity 做一致性 guard。source derivative unavailable 时仍可保留 valid inventory primal，但不补零、旧 Jacobian 或 hidden finite difference。此处的 inventory 不是 pore-volume accumulation，也没有 porosity/saturation/rock model。

[泡点／露点邻域回归](modules/flash/boundary_regression.md) 和[极近露点门槛审计](modules/flash/dew_limit_audit.md) 记录有限采样范围、独立高精度二元参考与相消失边界；原边界文档中的 TPD 停滞是历史记录，后续修复依据与回归见 [TPD 模块说明](modules/flash/README.md)。这些是软件和模型数值验证，不是全局稳定性证书。SW92 Profile-C 的物理三相验证还包括 [Sample-6 authoritative completion audit](modules/flash/sw92_authoritative_three_phase_audit.md) 与 [n-butane/H2O 实验三相线](modules/flash/sw92_three_phase_experimental_validation.md)。当前不提供 CPA flash、生产泡点／露点 API、PR76 standalone 单相 physics closure、validated LV/LL morphology classifier、salt-inventory conservation、pore-volume accumulation、production conservation residual/flux、网格/离散/time stepping/global Jacobian/流动求解，或跨泡露点/相消失/临界/family/root/topology switching 的闪蒸导数。

### 6.2 后续目标：多相与水处理约束

以下是后续扩展的目标约束。当前 Profile-C 已能在 **fixed-molality、有限搜索/拓扑契约** 下 authoritative 发布一相、两相或三相 PT phase set，并在固定 accepted phase set 的光滑局部域内提供 implicit Jacobian、physics thermodynamic closure 与 total-fluid-volume component inventory/local Jacobian；这不等于已实现盐库存守恒、数学全局稳定性认证、validated LV/LL morphology、pore-volume accumulation、production conservation residual/flux/离散或 flow coupling。

首期以给定温度、压力、总组成的非反应 `T-p-z` 平衡闪蒸为基础，默认相间同温同压；毛细压导致的相压差等扩展需另行定义。算法必须允许最终返回稳定单相、两相或三相，而不是强制返回指定数量的非零相。

两相包含需要验证的气液平衡（VLE）和液液平衡（LLE）；三相目标为气相、烃富液相和水富液相共存。水是化学组分，水富液相是相态，二者不得混为一谈；相标签也不能仅由 EOS 根的排列硬编码。

| 水处理模式 | 模型约定 |
| --- | --- |
| `immiscible_water` | 水作为独立相处理，明确禁止水与其他组分跨该相交换；烃体系平衡与独立水相的物性、质量守恒分别处理。不得称其为全组分互溶平衡。 |
| `mutual_solubility` | 水可进入模型允许的气相和烃富相，其他允许组分也可进入水富相；在同一组分守恒与一致的平衡条件下联合求解，而非闪蒸后再拼接固定水量。 |

水处理模式、可迁移组分和相集合必须显式配置并随结果保存。模型不具备所选模式所需能力时直接报错。这里的相间传质首先指平衡分配；有限速率传质是后续流动物理问题，不由平衡闪蒸自动实现。

对实际存在的相，验收至少检查以下关系及其数值容差：[S4][S5]

```text
sum(alpha, beta_alpha) = 1
sum(i, x_i_alpha) = 1
z_i = sum(alpha, beta_alpha * x_i_alpha)
f_i_alpha = f_i_beta  （适用于允许跨相迁移且实际参与平衡的组分）
beta_alpha >= 0，x_i_alpha >= 0
```

其中 `beta` 为摩尔相分率，`x` 与 `z` 为摩尔组成。迹量或零含量组分、消失相不能靠无条件取对数或硬裁剪处理。除守恒、组成与逸度/化学势一致性外，还要进行相稳定性判定，例如有明确搜索范围与容差的切平面距离（TPD）分析；求得一个方程根不等于证明全局稳定。[S4]

三相不得用两个互不一致的两相结果拼凑。SW92 Whitson dual-model observable Profile 仍只是两个独立 model pass，不得冒充共同守恒结果。Xu-style Gate 3B.3 已闭合 maximum-two-phase 的共同 EOS-component material balance、common chemical potentials、lower-envelope family selection、candidate Gibbs 选择与最终 AQ/NA common-tangent finite stability，因此可以发布**在声明的 max2/finite-search 边界内**的 family-aware 一相或两相 phase set；但这仍不证明不存在第三相，不守恒盐库存，也不自动赋予 water-rich/hydrocarbon-rich 物理标签。Profile-C 已进一步闭合 unordered `W(AQ)+H0(NA)+H1(NA)` joint equations、W-present final H-multiplicity review、disappearance-neighbor fresh re-solve 与 publication provenance guard，因此 `solve_sw92_profile_c_pt_phase_set(...)` 可以发布当前 fixed-molality contract 内的 authoritative 1/2/3-phase set，并可在固定 accepted phase set 上提供局部 implicit Jacobian 和 physics closure。这里的 authoritative 含义严格限定为“声明的 finite-search/topology contract 已接受”；H0/H1 仍不被命名为 liquid/vapor，`global_stability_proven=false`，NaCl 仍是 prescribed molality 而非 conserved salt inventory。失败必须返回结构化诊断，区分输入/参数错误、模型越界、未收敛、稳定性未确认与成功状态。近平衡退化点、临界区和相切换处的导数有效性须明确，不把固定相分支导数解释为跨相边界的光滑导数。

## 7. 递进路线与验收关口

| 阶段 | 交付范围 | 进入下一阶段前的证据 |
| --- | --- | --- |
| M0：文档引导 | 本 README 与根级 AGENTS（已完成初始化） | 初始化时仅文档；后续各阶段单独提供实现和验证证据。 |
| M1：AD 最小基础 | 已提供最小构建与测试入口、固定维数值/导数语义、基础算术 | 解析导数、边界行为和必要编译验证；每项分开小步提交。 |
| M2：AD 可用能力 | 已提供常见初等函数及逐函数测试入口；继续按需补齐函数和导数访问 | 独立交叉验证、导数保真与必要资源检查。 |
| M3：PR 与两相基线 | 已提供有序组分/PR76 参数契约、纯/混合系数、PT 候选相性质、有限 TPD 搜索、汽液 PT 相分裂、模型无关 phase-set 表示、内部两相收敛解灵敏度及最小 physics thermodynamic closure；已集成边界回归、TPD 停滞修复、露点门槛审计及诊断摘要 | 已有解析与独立高精度二元/三元参考的守恒、逸度、边界、未确定状态、隐式导数与 closure 链式导数回归；进一步扩展前仍需适用体系证据，不作全局、跨相边界或普适实验验证声明。 |
| M4：含水模型与互溶 | 已实现 SW92 corrected-original double/AD-capable thermodynamics、fixed-family stability、one-family VLE、Whitson dual-model compatibility observables、Xu-style Gate 3A + Gate 3B.1/3B.2/3B.3 maximum-two-phase joint path，以及 Profile-C C1/C2a1/C2a2/C2b.1/C2b.2 与 boundary-aware topology closure | 已有可追溯 SW92 参数、二元高精度 thermodynamics/stability/VLE/dual-model observable、Xu max2 回归、Profile-C AQ/NA joint anchors、additional-H witness、三相 equations、final finite review 与 phase-disappearance fresh-neighbor evidence。 |
| M5：三相闪蒸与局部 closure/inventory | **当前 fixed-molality Profile-C authoritative 1/2/3-phase PT publication、fixed-phase-set implicit sensitivity、variable-cardinality physics thermodynamic closure，以及同时支持 PR76 fixed-VLE/SW92 phase-set 的 model-neutral PT component inventory/local Jacobian 已实现**；后续扩展聚焦 validated morphology、salt inventory、pore-volume accumulation 与 production conservation residual/flux/Jacobian assembly | 已有独立 Decimal(80) structural 三相 reference、Mortezazadeh–Rasaei Sample-6 真实 8 组分三相物理锚点、Reamer et al. n-butane/H2O 实验三相线、fresh disappearance-neighbor routing、publication provenance guard、component permutation、H-slot symmetry、AD/IFT sensitivity checks、physics atomic-handoff regression、phase-slot/component permutation inventory checks，以及 PR76/SW92 0D fixed-volume fresh-resolve Jacobian cross-check；GCC Debug+ASan/UBSan、Clang Release、MSVC Release 对 inventory 与受影响 closure 下游验证。authoritative 仍不等于 global proof，H morphology 仍 unresolved，inventory 仍非 pore-volume accumulation。 |
| M6：网格与输入 | 创建、导入、拓扑与几何校验，逐一支持格式 | 小型已知网格、非法输入、单位与索引检查。 |
| M7：流动与数值计算 | 水独立/互溶物理模型、离散、时间推进、求解器，逐项开发 | 守恒、制造解或解析解、收敛性与相关耦合回归。 |
| M8：应用与结果 | 任务 API、前端、可视化、结果导出，按可信计算能力接入 | 契约与失败路径、必要读写往返、可视化数据一致性。 |

这是依赖关系与验收路线，不是一次性开发任务或时间承诺。服务接口与前端可在科学核心接口稳定后以独立增量接入；导出适配也可在验证需求出现时提前建设最小范围。不得为展示界面而填入冒充计算结果的数据。

## 8. 验证与协作

所有正式自动化测试使用 **GitHub 官方托管 runner（GitHub-hosted runner）**，不以 self-hosted runner 作为项目测试要求。[E7] 仅运行可能影响正确性或稳定性的必要增量验证，但测试影响范围按依赖传播判定，不能只按修改文件名决定。

纯文档变更仅做必要的内容、格式、引用和差异审查；AD、共享数值核心、公共接口、构建或编译选项变更必须扩展到受影响的已有下游。跨平台风险变更运行必要平台组合；平台未经实际测试不得宣称通过。共享托管 runner 上的计时不能单独证明 HPC 性能收益。

文档初始化阶段没有创建代码或 CI。当前 AD、热力学数据契约、PR76 纯/混合系数、PT 候选相性质核、TPD 稳定性搜索、汽液 PT 基线、模型无关 phase-set 表示、内部两相收敛解隐式灵敏度、SW92 corrected-original thermodynamics、fixed-family stability / one-family VLE、Whitson dual-model observables、Xu-style Gate 3A + Gate 3B.1/3B.2/3B.3 max2、Profile-C C1/C2a1/C2b.1/C2b.2、boundary-aware disappearance-neighbor re-solve、authoritative phase-set publication、fixed-phase-set sensitivity、PR76 physics closure、SW92 Profile-C physics closure 与 model-neutral PT component inventory/local Jacobian 均有增量测试入口与官方 runner 工作流；PR76 相分裂集成覆盖边界、停滞、露点门槛和诊断回归，并保留独立高精度二元/三元参考，灵敏度和 closure 另有独立 Decimal 高精度导数交叉核验。SW92 thermodynamics、fixed-family stability、one-family VLE、dual-model Profile、Gate 3A 与 Gate 3B.1 分别保留独立 Decimal(80) 数值锚点/回归；Gate 3B.2 验证 family-neutral witness plan、mandatory same-family alternatives、attempt-limit withholding、slot-swap dedup、distinct-candidate Gibbs near-tie 与真实 CO2/H2O full-plan selection；Gate 3B.3 独立 Decimal(80) 继续重算 lower-feed/pair Gibbs、selected-pair common tangent 与 prescribed AQ/NA TPD，并由 C++ regression 覆盖 accepted two-phase/single-phase、mandatory final starts、final negative/indeterminate rejection、pair-above-feed rejection、component permutation、资源配额与 publication guard。Profile-C C1 独立 Decimal(80) 直接解 phase-assigned AQ+NA cross-family equations；C2a1 独立 reference 重算 C1 tangent 上 retained-H TPD；C2b.1 独立 Decimal(80) reference 验证 generalized three-phase equations；Mortezazadeh–Rasaei Sample-6 回归进一步验证真实 8 组分 W+H0+H1 interior 状态、共同化学势、相分率和 runtime component permutation；Reamer et al. n-butane/H2O 回归独立验证实验三相压力与共存组成；boundary suite 验证 phase disappearance 必须 fresh re-solve 邻接 topology；authoritative publication suite 验证 1/2/3 phase projection、Sample-6 phase property/activity、incipient-W routing、错误 provenance、component permutation、H-slot symmetry 与 public-header self-containment；Profile-C sensitivity suite 以 AD 构造 fixed-set local equations 并结合 fresh re-solves、Sample-6、component permutation/H-slot symmetry 与 derivative boundary guards 验证 `beta/x/Z/c` Jacobian；SW92 physics closure suite 则验证 1/2/3-phase atomic primal/linearization handoff、unavailable failure policy 与 generic reduced-feed v2 contract；component-inventory suite 验证 fixed-VLE/variable-cardinality phase-resolved `c_mix/a_i`、phase-slot/component permutation、source-derivative unavailable 传播、malformed payload guard，并以 test-only 0D fixed-fluid-volume residual 对 PR76/SW92 fresh symmetric re-solves 交叉检查 local Jacobian。上述 accepted finite-search results 均不构成数学全局稳定性证明。CPA、PR76 standalone 单相 physics closure、validated H morphology、salt-inventory conservation、pore-volume accumulation、production conservation residual/flux/global Jacobian assembly、网格、离散与全局求解器尚未实现。具体执行结果查阅对应 PR/提交的 Actions 日志，不能将配置了工作流视为测试已通过。完整审计、提交与验证规则以 [AGENTS.md](AGENTS.md) 为准。

## 9. 待确认事项与参考资料

**SW 已由用户确认为 Søreide–Whitson 模型。** 当前已对 `SW92/corrected-original/PR76-base/NaCl-molality` thermodynamics kernel 核验并实现所需原文公式与作者勘误；fixed-family root/Gibbs stability、one-family PT VLE、`SW92-equilibrium/whitson-dual-model-observables/v1`、`SW92-equilibrium/xu-asymmetric-gibbs/v1` Gate 3A/3B.1/3B.2/3B.3 max2，以及 Profile-C C1/C2a1/C2a2/C2b.1/C2b.2/boundary-aware/publication 均已单独审计、实现并保留各自算法身份与验证边界。Profile-C authoritative publication 之后还已增加 fixed-phase-set AD/IFT sensitivity 与 SW92 variable-cardinality physics thermodynamic closure；二者均严格绑定 fixed accepted phase set、AQ/NA family、selected roots 与 representation slots，不把离散 topology selection 或 morphology classifier 当作可微映射。其上层 model-neutral PT component inventory 已统一消费 PR76 fixed-VLE 与 SW92 phase-set closure，并仅把已有 `beta/x/c` 及其 local derivative 链式组合为 total-fluid-volume inventory/Jacobian；它不是新 EOS/flash model，也不是 pore-volume conservation discretization。Profile-C C2a2 确认任意 H phase 的 physical liquid/vapor morphology 不能从 NA family、absolute Z 或 cubic-root index 自动推出；authoritative phase-set 因而仍将 H0/H1 表示为 unordered/representation-canonicalized `nonaqueous_unclassified` NA phase instances，而不是 W+L+V 或 W+L1+L2 的 morphology 声明。Whitson dual-model Profile 采用当前工程资料支持的“两次独立 AQ/NA model pass”编排模式，但仓库没有引入当前 refresh 项目的新拟合相关式、参数、salinity 方法或第三方源代码；这些若未来采用必须作为新的模型/数据 profile 重新核验。Xu-style max2 可以在同一个 EOS-component inventory/common chemical-potential 问题内发布 family-aware 一相或两相有限搜索结果；Profile-C 现在也可以在自己的 fixed-role/fixed-molality finite topology contract 内发布 authoritative 1/2/3 phase `PtPhaseSetResult`，并为该 accepted state 提供局部 sensitivity/physics closure/inventory。两者的 `global_stability_proven=false`，固定 molality 仍不代表 salt-inventory conservation；不能把 Whitson dual-model observables、Xu joint result 与 Profile-C phase-assigned result 相互替代，也不能把 Profile-C fixed-set Jacobian 套用到 Xu route 或离散 phase/topology transitions。首次使用其他模型前仍需取得对应原文、必要参数与合法验证资料；资料暂缺只阻塞依赖该资料的实现，不阻塞无关工作。项目许可尚未确定，不能因参考成熟开源软件就自动采用其代码或许可证。

工程参考用于借鉴边界、接口和验证方法，不代表这些软件已成为依赖。原始科学文献列为后续核验入口，不代表所有模型、所有公式与所有数据表均已获取；对已实现的 SW92 corrected-original profile，其核验范围和来源见 [SW92 thermodynamics 文档](modules/thermodynamics/sw92.md)，固定 family stability、一模型 VLE、双模型 observable、asymmetric stability、Gate 3B.1 fixed-pair、Gate 3B.2 orchestration 与 Gate 3B.3 max2 分别见 [SW92 stability](modules/flash/sw92_stability.md)、[SW92 one-family VLE](modules/flash/sw92_family_vle.md)、[SW92 dual-model Profile](modules/flash/sw92_dual_model.md)、[SW92 asymmetric stability](modules/flash/sw92_asymmetric_stability.md)、[SW92 Gate 3B.1 fixed pair](modules/flash/sw92_asymmetric_fixed_pair.md)、[SW92 Gate 3B.2 orchestration](modules/flash/sw92_asymmetric_orchestration.md) 与 [SW92 Gate 3B.3 max2](modules/flash/sw92_asymmetric_max2.md)。Profile-C 相关边界与实现见 [C1 joint header](modules/flash/include/mpmc/flash/sw92_phase_assigned_joint.hpp)、[H-side NA witness header](modules/flash/include/mpmc/flash/sw92_phase_assigned_h_side_witness.hpp)、[C2a2 morphology audit](modules/flash/sw92_phase_assigned_c2a2_hydrocarbon_role_audit.md)、[C2b.1 unordered three-phase candidate](modules/flash/sw92_phase_assigned_three_phase.md)、[boundary fresh re-solve](modules/flash/sw92_phase_assigned_boundary.md)、[authoritative completion audit](modules/flash/sw92_authoritative_three_phase_audit.md)、[experimental three-phase validation](modules/flash/sw92_three_phase_experimental_validation.md)、[Profile-C authoritative phase-set publication](modules/flash/sw92_profile_c_phase_set.md)、[Profile-C implicit sensitivity](modules/flash/sw92_profile_c_sensitivity.md)、[physics thermodynamic closure](modules/physics/README.md) 与 [PT component inventory/local Jacobian](modules/physics/component_inventory.md)。

- [E1] [CMake Presets 官方文档](https://cmake.org/cmake/help/latest/manual/cmake-presets.7.html)。
- [E2] [React：从头构建应用](https://react.dev/learn/build-a-react-app-from-scratch)；[Vite 官方指南](https://vite.dev/guide/)。
- [E3] [gRPC 核心概念](https://grpc.io/docs/what-is-grpc/core-concepts/)；[Protocol Buffers 语言指南](https://protobuf.dev/programming-guides/proto3/)。
- [E4] [gRPC-Web 官方实现与流支持说明](https://github.com/grpc/grpc-web#streaming-support)。
- [E5] [vtk.js 官方文档](https://kitware.github.io/vtk-js/docs/)。
- [E6] [CppAD 官方文档](https://coin-or.github.io/CppAD/)，用于 AD 模式、接口与独立核验方法参考，不预先绑定实现路线。
- [E7] [GitHub-hosted runners 官方文档](https://docs.github.com/actions/reference/runners/github-hosted-runners)。
- [E8] [C++ Core Guidelines](https://isocpp.github.io/CppCoreGuidelines/CppCoreGuidelines)；[OPM 模块说明](https://opm-project.org/?page_id=274)；[DuMux 官方文档](https://dumux.org/docs/doxygen/master/)。借鉴规范、分层与科学软件工程，复用代码前另行审计许可。
- [S1] Peng, D.-Y.; Robinson, D. B. (1976). *A New Two-Constant Equation of State*. DOI: [10.1021/i160057a011](https://doi.org/10.1021/i160057a011)。
- [S2] Søreide, I.; Whitson, C. H. (1992). *Peng-Robinson predictions for hydrocarbons, CO2, N2, and H2S with pure water and NaCl brine*. DOI: [10.1016/0378-3812(92)85105-H](https://doi.org/10.1016/0378-3812(92)85105-H)。对应当前已实现的 SW92 corrected-original thermodynamics profile；具体已核验公式、勘误、参数语义与实现边界见 [SW92 thermodynamics 文档](modules/thermodynamics/sw92.md)。
- [S3] Kontogeorgis, G. M.; Voutsas, E. C.; Yakoumis, I. V.; Tassios, D. P. (1996). *An Equation of State for Associating Fluids*. DOI: [10.1021/ie9600203](https://doi.org/10.1021/ie9600203)。
- [S4] Michelsen, M. L. (1982). *The isothermal flash problem. Part I. Stability*. DOI: [10.1016/0378-3812(82)85001-2](https://doi.org/10.1016/0378-3812(82)85001-2)。
- [S5] Michelsen, M. L. (1982). *The isothermal flash problem. Part II. Phase-split calculation*. DOI: [10.1016/0378-3812(82)85002-4](https://doi.org/10.1016/0378-3812(82)85002-4)。
- [S6] Xu, G.; Haynes, W. D.; Stadtherr, M. A. (2005). *Reliable Phase Stability Analysis for Asymmetric Models*. DOI: [10.1016/j.fluid.2005.06.016](https://doi.org/10.1016/j.fluid.2005.06.016)。当前 Xu-style 路径采用其 asymmetric lower-envelope/common-tangent formulation 构建 Gate 3A 与 maximum-two-phase Gate 3B；实现使用有限 multistart，不继承论文 interval-analysis 的全局保证。
