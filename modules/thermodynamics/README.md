# 热力学模块

`mpmc::thermodynamics` 提供有序组分/来源基础，以及 PR76、SW92、CPA 的参数、物性和各自声明的导数能力。它不依赖闪蒸或传输层；相稳定性与相分裂由上层 [flash](../flash/README.md) 实现。

## 参数与物性入口

| 范围 | 文档 | 主要边界 |
| --- | --- | --- |
| 公共组分/来源、PR76 参数 | [有序组分与 PR76 参数契约](pr_parameters.md) | 按稳定 ID 建立不可半更新的有序快照，记录来源、单位、版本、适用范围；缺失参数与合法零值分离。 |
| PR76 纯组分 | [a(T)、b 数值核](pr76_pure.md) | PR76 原始 alpha 与声明的 SI 常数，含温度导数。 |
| PR76 混合 | [经典混合规则](pr76_mixture.md) | 显式常数对称 kij；完整/约化组成坐标与运行期变维数。 |
| PR76 PT 物性 | [候选相性质](pr76_phase.md) | Z 根、ln(phi)、简单根局部导数与退化诊断；根数不等于相数。 |
| SW92 | [Corrected-original profile](sw92.md) | Søreide–Whitson 原文及作者勘误；AQ/NA family 与根分支分离，NaCl 为给定 molality。 |
| CPA 参数/缔合 | [CPA 基线](cpa.md) | 显式立方参数、位点与交叉缔合记录；不自动继承 PR 参数。 |
| CPA PT 物性 | [密度根与相性质](cpa_pt_phase.md) | 有限密度根扫描、机械可接受性与物性失败语义。 |

各专题维护模型身份、公式、来源、精度及验证范围。数值输入通过校验、位于声明区间，不等于物理模型已获实验验证；不同 EOS 的参数与 profile 不可只改名称后互换。

## 构建与下游

本模块目标及配置见 [CMakeLists.txt](CMakeLists.txt)。独立测试按实际改动选择：

| 契约/内核 | 独立 CMake 测试目录 |
| --- | --- |
| 公共/PR76 数据契约 | [contracts](../../tests/thermodynamics/contracts/) |
| PR76 纯组分 / 混合 / PT | [pr76](../../tests/thermodynamics/pr76/)、[pr76_mixture](../../tests/thermodynamics/pr76_mixture/)、[pr76_pt](../../tests/thermodynamics/pr76_pt/) |
| SW92 | [sw92](../../tests/thermodynamics/sw92/) |
| CPA 缔合 / PT | [cpa_baseline](../../tests/thermodynamics/cpa_baseline/)、[cpa_pt_phase](../../tests/thermodynamics/cpa_pt_phase/) |

具体命令、参考生成器及必要下游范围在对应专题中。公共物性或导数修改须覆盖受影响的 stability、flash 和 physics 回归；规则见 [AGENTS.md](../../AGENTS.md)。原参数契约文档中的阶段性“尚未实现”和资料缺口保留为历史记录，当前 SW92/CPA 状态以各自专题为准。

## Selected-phase fugacity contract for compositional flow

`selected_phase_fugacity.hpp` is the thermodynamics-owned bridge used by later natural-variable flow code. It does not select phases or roots. Callers must supply an already chosen branch:

- PR76: selected algebraic root index;
- SW92: selected family, NaCl molality and algebraic root index;
- CPA: selected PT density-root index.

The common output is ordered `ln(phi_i)` at the exact supplied `p,T,x`. PR76 and SW92 advertise `scalar_generic_first_order` because their selected-root phase kernels already preserve the caller scalar type and reject unreliable root derivatives. CPA currently advertises `value_only`: its PT density-root/association path is scientifically valid for fugacity values, but the repository does not yet contain the association-state and density-root IFT derivatives required for a complete `p,T,x` Jacobian. The contract therefore forbids presenting CPA as differentiable until that missing derivative layer is implemented and independently validated.

CPA `near_multiple`/tangent root topology is rejected by the selected-phase façade, and an out-of-range selected root is an error; the bridge never silently changes root identity.
