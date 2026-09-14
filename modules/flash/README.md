# PT 闪蒸模块

`mpmc::flash` 提供有限 TPD 相稳定性搜索、两相/三相平衡、相集合发布和模型无关 PT backend。PR76、SW92 Profile-C、CPA 均已接入统一后端；固定分支局部导数由 opt-in `mpmc::flash_sensitivity` 提供。热力学公式与参数由 [thermodynamics](../thermodynamics/README.md) 负责。

本页是导航；每个专题维护自己的接口、数值门槛、来源和验证证据。原 README 中的完整 TPD 契约、首轮验证和 PR #12 回溯修复记录已保留在 [pt_stability.md](pt_stability.md)，原路径中的历史代码注释仍可经本页到达该文档。

## 公共契约与 PR76

| 内容 | 入口 |
| --- | --- |
| 统一 configured backend、能力发现与冻结的 v1 handoff | [PT flash backend](pt_flash_backend.md) |
| 通用相集合、相变边界与诊断 | [Phase set](phase_set.md)、[phase transition](pt_phase_transition.md)、[diagnostic summary](diagnostic_summary.md) |
| 有限 TPD 搜索、PR76 适配与回溯修复历史 | [PT stability](pt_stability.md) |
| 汽液平衡、接受条件与隐式导数 | [PT split](pt_split.md)、[PT sensitivity](pt_sensitivity.md) |
| PR76 maximum-three-phase 与连续状态路径 | [Three-phase baseline](pr76_three_phase.md)、[PT continuation](pr76_pt_continuation.md) |
| 二元边界、露点与独立三元参考 | [Boundary regression](boundary_regression.md)、[dew-limit audit](dew_limit_audit.md)、[ternary regression](ternary_regression.md) |

## SW92：按算法 profile 查阅

| 路径 | 用途与入口 |
| --- | --- |
| 算法身份与科学范围 | [Equilibrium algorithm](sw92_equilibrium_algorithm.md)、[profile separation audit](sw92_equilibrium_profile_separation_audit.md) |
| Fixed-family primitive | [Stability](sw92_stability.md)、[one-family VLE](sw92_family_vle.md)；单次运行固定一个 AQ/NA family。 |
| Whitson dual-model observables | [Dual-model profile](sw92_dual_model.md)；两次独立模型运行的兼容性 observable，不是联合守恒相集合。 |
| Xu-style asymmetric max2 | [Asymmetric stability](sw92_asymmetric_stability.md)、[fixed pair](sw92_asymmetric_fixed_pair.md)、[orchestration](sw92_asymmetric_orchestration.md)、[max2 acceptance](sw92_asymmetric_max2.md)；使用 lower-envelope family model。 |
| Profile-C 生产入口 | [Authoritative 1/2/3-phase publication](sw92_profile_c_phase_set.md)、[fixed-phase-set sensitivity](sw92_profile_c_sensitivity.md)；固定 `W->AQ`、`H->NA` 与 prescribed molality。 |
| Profile-C 求解与边界 | [Joint W+H](sw92_phase_assigned_joint.md)、[additional-H witness](sw92_phase_assigned_h_side_witness.md)、[no-W topology](sw92_phase_assigned_no_w.md)、[unordered W+H0+H1](sw92_phase_assigned_three_phase.md)、[three-phase closure](sw92_phase_assigned_three_phase_closure.md)、[PT driver](sw92_phase_assigned_pt.md)、[boundary fresh re-solve](sw92_phase_assigned_boundary.md) |
| Xu max2 验证 | [Independent numerical validation](sw92_asymmetric_max2_validation.md)、[physical ternary regression](sw92_asymmetric_max2_physical_ternary.md) |
| Profile-C 三相证据 | [Authoritative completion / Sample-6 audit](sw92_authoritative_three_phase_audit.md)、[experimental three-phase line](sw92_three_phase_experimental_validation.md) |

### 保留的设计与历史审计

这些文档记录各增量当时的前置条件、否决方案和证据。阶段性的“下一步/未实现”应结合上表中的后续实现阅读，不把不同 profile 的结论互相替代。

| 审计主题 | 入口 |
| --- | --- |
| Xu formulation 与 joint split | [Xu asymmetric audit](sw92_xu_asymmetric_audit.md)、[Gate 3B design](sw92_gate3b_joint_split_audit.md) |
| Profile-C stability/topology | [C2 stability audit](sw92_phase_assigned_c2_stability_audit.md)、[C2a topology](sw92_phase_assigned_c2a_topology_audit.md) |
| H morphology 与后续拓扑审查 | [C2a2 hydrocarbon role](sw92_phase_assigned_c2a2_hydrocarbon_role_audit.md)、[post-C2b.1 topology/morphology](sw92_phase_assigned_post_c2b1_topology_morphology_audit.md) |

## CPA

| 内容 | 入口 |
| --- | --- |
| Minimum-Gibbs PT stability | [CPA stability](cpa_stability.md) |
| VLE、max3、最终复核与降相重新求解 | [CPA flash](cpa_flash.md) |
| 甲醇/水两相文献参数与实验验证 | [CPA physical validation](cpa_physical_validation.md) |

CPA 与 PR76 的 max3 synthetic structural fixtures 不等于真实流体的物理三相验证。CPA 仍需完整兼容的 associating VLLE 参数与独立三相 oracle。

## 共用使用边界与验证

所有当前有限搜索都保持 `global_stability_proven=false`。候选方程收敛不等于相集合已接受；`indeterminate` 不转换为成功。根侧、AQ/NA family 或数值槽位不自动赋予液/汽形态；固定相集合 Jacobian 不跨相出现/消失、根/family/topology 切换外推。具体状态和阈值在相应契约中维护。

模块构建目标见 [CMakeLists.txt](CMakeLists.txt) 和 [sensitivity/CMakeLists.txt](sensitivity/CMakeLists.txt)；各专题指向 [tests/flash](../../tests/flash/) 下的独立工程与参考脚本。正式测试按 [AGENTS.md](../../AGENTS.md) 选择必要增量及受影响下游，不因更新文档启动无关数值测试。

上层消费入口为 [physics closure](../physics/README.md)、[component inventory](../physics/component_inventory.md) 和 [PtService](../runtime/README.md)。它们各自维护消费契约，不在本页重复复制整段实现状态。
