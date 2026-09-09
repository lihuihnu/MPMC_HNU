# 汽液 PT 相分裂基线

在 [TPD 稳定性搜索](README.md) 上新增 C++20/double 相分裂：进料稳定性 → 负 TPD 试探组成 → 守恒初值 → Rachford–Rice / 对数 K 逐次代入 → 最终共同切平面复核。不是全局认证闪蒸器；相角色为明确请求的流体候选分支，不是由根数推断实际相数。

## 1. 审阅与写入前审计

本增量基于 PR #9 的 `3c09b8ac879d04b1bdbde0536ef5e93f02fed271`（tree `57efcfedb4978b34d8fbd45aeed714b49e027220`），主线仍为 `c48ac3a`。审阅意见已以 COMMENT 写入 PR #9，review `5150222616`；原作者复核不是第三方独立审稿，未自我 APPROVE、未合并或强推。

复核既有契约、根/逸度接口、失败记录、测试与依赖后，未发现阻止复用其有限局部搜索契约的必改生产缺陷。明确两个后续风险：`lowest_sampled` 不是可靠负证据接口，应读取 `negative_tpd` 的 trial.point；原入口仅支持进料参考，不能直接用于最终两相共同切平面。PR #9 正文将 Clang/MSVC job ID 交换，审阅评论已更正；实际三个作业仍全部 success。

全局取舍：不复制稳定性优化器，不把 PR 公式放入通用相分裂；只扩展原搜索的参考构造，新增独立 RR、相分裂与 PR 适配头。保留原 AD/thermodynamics 源码和旧测试。与原样保留接口但重复搜索、或一次引入全 Newton/Hessian/三相框架相比，共用驱动＋受保护逐次代入是本轮最小可验证范围。共享稳定性入口变化需要复跑原22项稳定性测试；新27项测试及受影响平台必须在官方runner运行。

科学核验：已查看 [R1] 期刊211–212页公式及 [R2] 作者稿10–11页参数、表4–5渲染页。原 Rachford–Rice 全文镜像本轮访问失败，公式依据已读的 [R1] 式(11)–(14)，不声称复现未读原文。下述守恒初值、指数缩放、误差判据与回溯是本实现独立推导，不复制第三方代码或引入新依赖。参数夹具不重新拟合，不以同一生产实现生成期望值。

## 2. 接口与对象语义

公共头：`rachford_rice.hpp`、`pt_split.hpp`、`pr76_split.hpp`，仍使用 `mpmc::flash` CMake目标。可独立 `add_subdirectory(modules/flash)`；根级 AD 构建入口不改。核心无网络、前端或 AD 依赖，不计算迭代过程导数或平衡解灵敏度。

```cpp
#include <mpmc/flash/pr76_split.hpp>
const auto phase = mpmc::thermodynamics::Pr76Phase<double>::from_parameters(parameters);
mpmc::flash::Pr76VleEvaluator evaluator(phase);
const auto result = mpmc::flash::solve_pr76_pt_vle(p_pa, t_kelvin, z, evaluator);
const auto& solution = result.solution;
if (const auto* pair = solution.candidate()) {
    // pair 是方程收敛的候选解；必须再读取 solution.status 和 final_stability。
    const double beta_v = pair->fractions.vapor_fraction;
    const auto& x = pair->fractions.liquid;
    const auto& y = pair->fractions.vapor;
    (void)beta_v; (void)x; (void)y;
}
```

p 单位 Pa，T 单位 K，beta 为摩尔汽相分率，x/y/z 为有序摩尔组成。正温压、有限性、组成和与支持集要求沿用稳定性契约。`solve_pt_vle` 明确记录并使用进料舍入范围内的归一化；不接收任意未归一化输入，不更改负数或添加缺失组分。

通用 `solve_pt_vle(p,T,z,stability_provider,phase_provider,options,initial_starts,final_starts)` 接受两个同步可重复回调，二者必须是同一模型、组分顺序和标准参考态。稳定性回调枚举所有允许分支的最低Gibbs值；相分裂回调接受 `PtPhaseRole` 并返回候选分支的 ln(phi)、Z。类型/数组校验无法验证回调的热力学一致性，此责任属于模型提供者。

`iterate_pt_split` 是给定 logK 初值的低层迭代，不包含初始或最终稳定性，不能把其 converged 单独解释为闪蒸成功。其结果保留最后接受的迭代点；只有 status=converged 同时通过方程、守恒与相区分要求。完整结果 `candidate()` 仅暴露被选中的方程收敛点，返回指针随结果对象存活；临时对象调用被禁止。结果拥有所有向量；适配器拥有参数快照，工作区只顺序复用，不能并发共用。

## 3. 失稳组成与守恒初始化

先用旧入口 `test_pt_stability` 对进料运行有限搜索。未确定时停止分裂并报告；未检出时返回单相候选及其有限搜索结论，不强行创建第二相。仅使用越过原阈值/guard的 `negative_tpd` 试探点 w。

归一化 w 本身不是平衡相，也不能直接把 K=w/z 当成有内部汽相分率的初值。选

```text
a = min(0.1, 0.5*min_active(z_i/w_i))
u_i = (z_i-a*w_i)/(1-a)
```

数学上 `z=a*w+(1-a)*u`，且互补相保留正活跃组成。u 只允许舍入范围内的归一化；不能表示时拒绝该初值。尝试 `(x=u,y=w)` 与 `(x=w,y=u)` 两种角色，再取 logK=log(y)-log(x)，不预先猜测失稳点属于气相还是液相。默认最多16个分裂尝试，所有尝试共享20000次单相物性调用预算；记录来源 trial 索引与角色，不能伪装为独立物理资料。受限未尝试部分以 attempt_limit_reached 记录。

## 4. Rachford–Rice 与迭代

[R1] 式(11)–(14)给出

```text
K_i = y_i/x_i = phiL_i/phiV_i
F(beta) = sum z_i*(K_i-1)/(1-beta+beta*K_i) = 0
x_i = z_i/(1-beta+beta*K_i), y_i = K_i*x_i
```

F 在非平凡区间单调递减。只在 `0<beta<1` 求根；端点符号用 log-sum-exp 的 `log(sum z*K)` 和 `log(sum z/K)` 判断，64*eps内视为未分辨，**不能据初始K的无内部根宣布稳定单相**。K=1返回degenerate。

内部计算用 `exp(-abs(logK))`、`expm1` 分别重排正/负logK的分子分母，不形成可能溢出的exp(+logK)。二分最多192次；保留RR残差，并检查原始x/y和在64*eps内。随后明确归一化供EOS使用，重新计算返回组成上的逐组分守恒。正trace下溢为0时报告不能表示，不填下限；无负闪蒸、根裁剪或越界beta修补。

外层残差为 `r_i=log(x_i/y_i)+lnphiL_i-lnphiV_i`，在全部正进料组分上取最大绝对值。方向为logK的逐次代入修正，经 `s=max(1,max|r|/2)` 限制实际对数步，再逐次减半尝试。接受条件为残差达到原容差，或满足 `norm_new < norm_old*(1-1e-4*alpha/s)`。这是固定点的残差保护，不是声称具有已证明全局收敛性的Newton方法；保护失败保留诊断。

默认最大512个外层更新、每步24次回溯；计数包含失败与拒绝的单相物性调用。内存分配/逻辑/输入域错误不吞掉；数值性质失败保留结构化问题，可在回溯中缩步恢复。默认工作区与结果规模受组分、稳定性初值/存储、尝试数和评估次数限制；这些不是精确字节配额，服务层仍需限制参数元数据。

## 5. 联合验收与最终相集合

候选方程收敛要求同时满足：

| 条件 | 默认阈值与含义 |
| --- | --- |
| 逸度 | `max_active abs(log(fL/fV)) <= 1e-11`，不加摩尔分数权重。 |
| 守恒绝对值 | `max_i abs((1-beta)*x_i+beta*y_i-z_i) <= 1e-12`。 |
| 守恒相对值 | 对每个正z_i用其自身作分母，最大相对误差<=1e-10；无绝对下限掩盖trace。 |
| 相分率 | beta与1-beta均>1e-10；更小视为相消失未分辨，不硬设零再称成功。 |
| 相区分 | max abs(logK)>1e-7，且 `ZV-ZL>1e-8*max(ZV,ZL)`，排除同组成/同密度伪两相。 |

PR适配器从机械可用根选择请求的低Z/高Z候选，用原phase核计算性质；病态或未分辨根直接诊断。两种请求在单根区可能落到同一分支，不能仅因调用两次就算两相。密度排序不能证明任意混合物必为VLE，本基线只求给定的汽液候选分支；LLE、临界区、纯物质饱和时相分率不可确定等情形可以返回未确定。

方程收敛后选择已收敛候选中Gibbs最低的点，计算相对进料的Gibbs变化，并用

```text
d_i = midpoint(log(x_i)+lnphiL_i, log(y_i)+lnphiV_i)
D_set(w) = sum w_i*[log(w_i)+lnphi_i(w)-d_i]
```

对整个允许试探相集合重新搜索，不再使用进料化学势。新入口 `test_pt_stability_against` 共用旧驱动，以 `lnphi_equivalent=d-log(z)` 进行代数参考转换；等价参考不是真实进料相。结果存 `imposed_log_activity`，而 `reference` 留空；不虚构参考物性调用。原接口数值表达式、下降逻辑与容差保持不变。

复核初值至少包括两个已求得相、原默认组成搜索和调用者的额外初值。默认有效负值阈值为用户指定的最终TPD阈值加 `0.5*max|log(fL/fV)|`，明确记录为 common_reference_allowance 及 final_stability.options.tpd_tolerance；它表示两个数值参考之间的差异，**不是到精确平衡解的严格误差界**，也不改变相分裂的1e-11逸度标准。仍保留原算术guard，其亦不是EOS总误差界。

| 完整结果 | 严格含义 |
| --- | --- |
| single_phase_no_instability_found | 初始有限搜索未检出失稳；不提供全局单相证明。 |
| two_phase_no_instability_found | 不同的两相候选满足联合方程标准、Gibbs未显著高于进料，且共同切平面有限搜索未检出进一步失稳。 |
| phase_set_unstable | 两相方程收敛，但共同切平面发现更低Gibbs试探分支；保留候选和证据，拒绝两相成功，不自动拼接第三相。 |
| indeterminate | 初始未确定、无可用分裂、退化、预算耗尽，或方程已收敛但最终复核未完成。后者仍保留candidate、equations_converged=true。 |

所有结果的 global_stability_proven 恒为false。equations_converged 只描述已选候选方程，不替代status；仅有零逸度残差不足以证明全局平衡。Gibbs比较保护量来自明确算术幅值，不是严格区间认证。

## 6. 验证与边界

新增27项CTest覆盖：解析固定K的RR、零/迹量、极端logK、边界/退化、非法输入/预算、原参考与外部参考等价、原进料失稳但共同切平面未失稳的反例、守恒种子、精确与非精确SSI初值、数值/编程异常、单/两相候选、初始/最终未确定、第三候选分支失稳、相消失、PR两体系6个进料、顺序/所有权/根失败及公共头。因共用稳定性驱动改变，原22项亦属于本轮必要回归。

解析夹具用常数理想分支 `phiL=(2,1/4), phiV=(1,1)`，z=(.6,.4)的独立解析解为x_A=3/7、y_A=6/7、beta=.4。附加常数分支lnphi=(-1,2)在该两相共同切平面下产生负值，用于验证“方程收敛不等于相集合稳定”；这些是制造数学模型，不是物理流体。

PR参数仅来自[R2]，没有将其TPD驻点当作闪蒸平衡组成。`reference_split_decimal.py` 采用80位十进制、原Z三次式、直接PR逸度、在两个二元相组成中求解的独立Newton法，再按杠杆规则求beta；它不调用生产SSI或RR。参考中的有限差分仅用于独立Newton的导数，不是生产AD，也不是唯一正确性依据。脚本复核EOS残差与log逸度差，输出六状态各五个数值（部分相组成/Z因同温压二元系而重复）；与C++字面值核对相对1e-35，C++数值回归使用 `2e-9*(1+|expected|)`，联合守恒/逸度仍执行生产原阈值。普通解析比较使用 `2e-12*(1+|expected|)`。

此版本的文献测试是模型数值回归，不是实验验证。未实现认证全局搜索、联合三相、反应、电解质、固相、水合物、水互溶或水独立模式、SW/CPA；未验证macOS、其他架构、并发竞态或HPC性能。原有用户组分来源/适用范围检查仍由参数快照及物性模型执行，不自动认为普通PR适用于含水互溶。

```sh
cmake -S tests/flash/pt_split -B build/pt-split -DCMAKE_BUILD_TYPE=Debug
cmake --build build/pt-split --target mpmc_pt_split_tests --config Debug
ctest --test-dir build/pt-split -C Debug -R '^flash[.]pt_split[.]' --verbose --no-tests=error
python tests/flash/pt_split/reference_split_decimal.py
```

正式结果以具体提交的官方runner日志为准。新增workflow独立执行相分裂；原稳定性workflow仅将关联头路径收窄，并允许依赖PR的基分支。纯相分裂修改不再误触发原稳定性全套；共享稳定性变化仍运行两套。不会重跑未受影响的旧AD/thermo测试。两个PR保持开放，新PR以PR #9为依赖；未经授权不合并主线。

## 7. 原始依据

[R1] Gernert, J.; Jäger, A.; Span, R. (2014), Calculation of phase equilibria for multi-component mixtures using highly accurate Helmholtz energy equations of state, Fluid Phase Equilibria 375, 209–218. DOI `10.1016/j.fluid.2014.05.012`。本轮核验式(11)–(14)、(19)–(22)、(27)–(30)。[论文](https://www.sintef.no/globalassets/project/impacts/phase-equilibria-gernert-et-al-fpe-2014.pdf)。没有声称复现该论文全部初始化和Newton方案。

[R2] Hua, J. Z.; Brennecke, J. F.; Stadtherr, M. A., Reliable Computation of Phase Stability Using Interval Analysis: Cubic Equation of State Models，1997年10月修订作者稿，第4.3/4.4节、表4/5；甲烷参数源于前面第4.1节。[作者稿](https://academicweb.nd.edu/~markst/zm97a.pdf)。原文正文/表4最后一个N2进料.65/.60的差异保留于PR #9审计，本次平衡参考使用明确失稳的.18/.30/.44；CO2使用.20/.30/.43。未复制其区间算法、实验数据或全文。
