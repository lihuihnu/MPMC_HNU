# PR76 三组分汽液 PT 端到端回归

## 1. 写入前审计与范围

基线为 `main@230b8c523146e5eba918ba0c8cf7c4fe532f530c`。本增量只验证现有运行期组分链路在一个非平凡三组分 PR76 汽液 PT 问题上的端到端行为，不修改生产数值实现、公共状态枚举或默认容差。

审计链路为：

1. `PrParameterSet::create` 以运行期 `order.size()` 建立有序快照；所有被选异组分对必须显式提供，缺失不能自动解释为 `kij=0`。
2. `Pr76Mixture` 和 `Pr76Phase` 的存储、混合双重求和、组成长度检查及工作区均按运行期组分数处理，不含二元专用分支。
3. `test_pt_stability` 按活动组分数生成 feed、uniform 和逐活动顶点初值；有限多初值仍不构成全局稳定性证明。
4. `solve_pr76_pt_vle` 要求进料长度与参数快照一致，并把有序组分 ID、数据集版本和 PR76 约定带入结果。改变组分集合或顺序时重新建立快照/模型/求值器；不把一个固定模型原地热切换为另一维数。

因此本轮不需要修改生产头。新增三元参数/独立参考/端到端回归，并把它接入现有 PT phase split 测试工程与官方 runner。

## 2. 三元参数与可追溯性

体系为甲烷／乙烷／丙烷。模型保持项目既有
`PR76/classical-vdw/constant-kij/no-translation` 和
`PR76/printed-coefficients/R-SI-2019/exact-sqrt2/PT-v1` 约定。

纯组分参数取 Deiters 与 Bell 的公开论文 *Calculation of phase envelopes of fluid mixtures through parametric marching*，DOI `10.1002/aic.16730`，Table 1：

| 组分 | Tc / K | Pc / MPa | omega |
| --- | ---: | ---: | ---: |
| methane | 190.555 | 4.595 | 0.0 |
| ethane | 305.4 | 4.88 | 0.099 |
| propane | 369.825 | 4.248 | 0.15308 |

压力在输入记录中显式执行 `MPa * 1e6 -> Pa`。本回归**不使用**该论文为其甲烷／丙烷二元算例列出的拟合 `k12`。

三对交互参数分别为 methane/ethane、methane/propane、ethane/propane，均以有来源的**显式零记录**提供。依据 Peng 与 Robinson 原始论文 *A New Two-Constant Equation of State*，DOI `10.1021/i160057a011`，期刊 p.62 对 methane/ethane/propane 三元算例的说明：该例没有使用 interaction coefficients。这里把“未使用”落实为本项目常数对称 `kij` 契约中的三个显式零，而不是依赖缺失参数的隐式默认值。

这是一个组合的、可追溯的软件回归数据集：纯组分数值与零交互系数的依据来自两份不同文献，并不声称复现任一论文的全部实验条件、参数拟合或最终数值结果。

## 3. 独立参考

参考程序 `reference_ternary_decimal.py` 只使用 Python 标准库 `Decimal(80)`。固定：

- `T = 220 K`
- `p = 4,000,000 Pa`
- 液相 `x_CH4 = 0.5`

未知量取 `x_C2H6`、`y_CH4`、`y_C2H6`，直接求三组分

```text
ln(x_i) + ln(phi_i^L) = ln(y_i) + ln(phi_i^V), i=1..3
```

参考侧独立实现 PR76 纯组分参数、经典双重混合求和、Z 三次式和逸度公式，并在三维相组成空间用 Newton 法求解。数值 Jacobian只服务独立参考 Newton；不调用生产 C++、TPD、多初值搜索、logK-SSI 或 Rachford-Rice，也不从生产输出反推期望值。

收敛 tie-line 为：

```text
x = (0.5,
     0.3334151832635363826660932077339112222787,
     0.1665848167364636173339067922660887777213)

y = (0.8943594827590593429076380789028395017818,
     0.09378930148285197934244636668815656792870,
     0.01185121575808867774991555440900393028947)

ZL = 0.1251092543567705217556002938357282057135
ZV = 0.6883465259522418692089462285643535773927
```

Decimal 参考要求三组分化学势平衡残差 `<1e-60`。随后才用精确 Decimal 杠杆关系构造 `betaV=0.25/0.50/0.75` 三个内部进料状态，并把 36 个参考量写入只读头文件。脚本每次 CI 重新生成并逐项核对头文件到 `1e-35 * max(1, |value|)`。

这些状态是为数值回归选取的模型状态，不是实验数据，也不是相包络或全局稳定性认证。

## 4. C++ 端到端验收

新增四项独立 CTest：

- `ternary_parameters_runtime`：检查三个显式零 `kij` 及其来源；删除任一 pair 必须得到 `missing_parameter`，不能隐式补零；同一工作区顺序经历 `1 -> 3 -> 2 -> 3` 组分模型求值，验证重建快照后的运行期维数链路。
- `ternary_flash_reference`：默认初始化对三个进料执行完整 `solve_pr76_pt_vle`，必须完成初始失稳证据、相分裂和最终共同切平面复核；相组成、相分率、`ZL/ZV` 与独立参考比较。
- `ternary_permutations`：对中间进料执行全部 6 种组分排列；每种都重建有序参数快照，按组分 ID 映射后必须得到相同物理解。
- `ternary_final_indeterminate`：只把最终稳定性物性调用预算降到 1；已收敛两相候选必须保留，但顶层仍返回 `indeterminate`，验证既有“方程收敛不等于最终相集合已接受”的语义。

独立参考比较沿用现有二元 PR 回归的 `2e-9 * (1+|expected|)` 预算，不放宽生产停止条件。生产验收仍严格使用：

| 条件 | 既有值 |
| --- | ---: |
| 最大 `|ln(fL/fV)|` | `1e-11` |
| 逐组分绝对守恒 | `1e-12` |
| 正进料组分相对守恒 | `1e-10` |
| 最小相分率 | `1e-10` |
| TPD 驻点残差 | `1e-8` |
| 负 TPD 证据 | `1e-10` 加既有算术 guard |

测试还用新的 PR76 相性质调用在返回组成上重新计算逸度残差；只给不同运算路径 `32*epsilon(double)` 的算术余量，保存的生产残差仍必须严格小于原 `1e-11`。

## 5. 明确不变项与验证范围

原 `split_test.cpp`、全部二元 `golden` 字面值、`reference_split_decimal.py`、边界/露点参考、生产 `pt_split.hpp`/`pt_stability.hpp`/PR76/RR 头、状态枚举和默认门槛均不修改。

CI 只扩展既有 PT phase split 工程：编译新三元测试目标、运行 `^flash[.]pt_split[.]`，并新增 Decimal(80) 三元参考复核。原二元、边界、露点和诊断检查仍按原工作流执行；不因测试增量机械运行无关 AD 或独立热力学套件。正式 C++ 结果以本增量 PR 的 GitHub 官方 runner 日志为准，不能用参考脚本本身替代 C++ 验证。
