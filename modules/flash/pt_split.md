# 汽液 PT 相分裂基线

在 [TPD 稳定性搜索](README.md) 上新增 C++20/double 相分裂：进料稳定性 → 负 TPD 试探组成 → 守恒初值 → Rachford–Rice / 对数 K 逐次代入 → 最终共同切平面复核。不是全局认证闪蒸器；相角色为明确请求的流体候选分支，不是由根数推断实际相数。

## 1. 审阅与写入前审计

本增量基于 PR #9 的 `3c09b8ac879d04b1bdbde0536ef5e93f02fed271`（tree `57efcfedb4978b34d8fbd45aeed714b49e027220`），主线仍为 `c48ac3a`。审阅意见先以 COMMENT 写入 PR #9，review `5150222616`，再编写代码。原作者复核不是第三方独立审稿，未自我 APPROVE、未合并或强推。本增量位于依赖 #9 的 PR #10。

复核既有契约、根/逸度接口、失败记录、测试与依赖后，未发现阻止复用其有限局部搜索契约的必改生产缺陷。明确两个后续风险：`lowest_sampled` 不是可靠负证据接口，应读取 `negative_tpd` 的 trial.point；原入口仅支持进料参考，不能直接用于最终两相共同切平面。PR #9 正文将 Clang/MSVC job ID 交换，审阅评论已更正；实际三个作业仍全部 success。

全局取舍：不复制稳定性优化器，不把 PR 公式放入通用相分裂；只扩展原搜索的参考构造，新增独立 RR、相分裂与 PR 适配头。保留原 AD/thermodynamics 源码和旧测试。共用驱动＋受保护逐次代入避免本轮引入全 Newton/Hessian/三相框架。共享稳定性入口变化需要复跑原22项稳定性测试。首轮验证发现原外层回溯目标选择有缺陷，修正依据与真实失败见第7节。

科学核验：已查看 [R1] 期刊211–212页公式及 [R2] 作者稿第4.3/4.4节参数、表4–5渲染页。原 Rachford–Rice 全文镜像本轮访问失败，公式依据已读的 [R1] 式(11)–(14)，不声称复现未读原文。守恒初值、指数缩放、Gibbs方向导数与数值保护由本实现推导，不复制第三方代码或引入新依赖。参数夹具不重新拟合，不以同一生产实现生成期望值。

## 2. 接口与对象语义

公共头：`rachford_rice.hpp`、`pt_split.hpp`、`pr76_split.hpp`，仍使用 `mpmc::flash` CMake目标。可独立 `add_subdirectory(modules/flash)`；根级 AD 构建入口不改。核心无网络、前端或 AD 依赖，不计算迭代过程导数或平衡解灵敏度。

```cpp
#include <mpmc/flash/pr76_split.hpp>
const auto phase = mpmc::thermodynamics::Pr76Phase<double>::from_parameters(parameters);
mpmc::flash::Pr76VleEvaluator evaluator(phase);
const auto result = mpmc::flash::solve_pr76_pt_vle(p_pa, t_kelvin, z, evaluator);
const auto& solution = result.solution;
if (const auto* pair = solution.candidate()) {
    // 方程收敛的候选解；仍须读取 solution.status 和 final_stability。
    const double beta_v = pair->fractions.vapor_fraction;
    const auto& x = pair->fractions.liquid;
    const auto& y = pair->fractions.vapor;
    (void)beta_v; (void)x; (void)y;
}
```

p 单位 Pa，T 单位 K，beta 为摩尔汽相分率，x/y/z 为有序摩尔组成。正温压、有限性、组成和与支持集要求沿用稳定性契约。完整驱动记录并使用进料舍入范围内的归一化；不接收任意未归一化输入，不更改负数或添加缺失组分。

通用 `solve_pt_vle(p,T,z,stability_provider,phase_provider,options,initial_starts,final_starts)` 接受两个同步可重复回调，二者必须是同一模型、组分顺序和标准参考态。稳定性回调返回全部允许分支中的最低Gibbs值；相分裂回调接受 `PtPhaseRole` 并返回候选分支的 ln(phi)、Z。类型/数组校验无法验证回调的热力学一致性，此责任属于模型提供者。

`iterate_pt_split` 是给定 logK 初值的低层迭代，不包含初始或最终稳定性，不能把其 converged 单独解释为闪蒸成功。其结果保留最后接受的迭代点；只有 status=converged 同时通过方程、守恒与相区分要求。完整结果 `candidate()` 仅暴露被选中的方程收敛点，返回指针随结果对象存活；临时对象调用被禁止。公开索引被调用者改成越界、point清空或attempt状态改为未收敛时，访问器返回nullptr，不进行无效解引用；这不是对调用者任意改写数值后的重新科学验证。结果拥有所有向量；适配器拥有参数快照，工作区只顺序复用，不能并发共用。

## 3. 失稳组成与守恒初始化

先用旧入口对进料运行有限稳定性搜索。未确定时停止分裂并报告；未检出时返回单相候选及其有限搜索结论，不强行创建第二相。仅使用越过原阈值/guard的 `negative_tpd` 试探点 w。

归一化 w 本身不是平衡相，直接取 K=w/z 也不保证存在内部汽相分率。选

```text
a = min(0.1, 0.5*min_active(z_i/w_i))
u_i = (z_i-a*w_i)/(1-a)
```

数学上 `z=a*w+(1-a)*u`，且互补相保留正活跃组成。u 只允许舍入范围内的归一化；不能表示时拒绝初值。尝试 `(x=u,y=w)` 与 `(x=w,y=u)` 两种角色，再取 logK=log(y)-log(x)，不预先猜测失稳点属于气相还是液相。默认最多16个分裂尝试，共享20000次单相物性调用预算；记录来源 trial 索引与角色。受限未尝试部分以 attempt_limit_reached 记录。

## 4. Rachford–Rice 与 Gibbs 下降保护

[R1] 式(11)–(14)给出

```text
K_i = y_i/x_i = phiL_i/phiV_i
F(beta) = sum z_i*(K_i-1)/(1-beta+beta*K_i) = 0
x_i = z_i/(1-beta+beta*K_i), y_i = K_i*x_i
```

F 在非平凡区间单调递减。只在 `0<beta<1` 求根；端点符号用 log-sum-exp 的 `log(sum z*K)` 和 `log(sum z/K)` 判断，64*eps内视为未分辨。初始K无内部根不能宣布稳定单相，K=1返回degenerate。

内部计算用 `exp(-abs(logK))`、`expm1` 重排正/负logK的分子分母，不形成exp(+logK)。二分最多192次；检查RR残差以及原始x/y和在64*eps内。随后明确归一化供EOS使用，并重新检查返回组成上的逐组分守恒。正trace下溢为0时报不能表示，不填下限；无负闪蒸或beta裁剪。

外层残差为 `r_i=log(x_i/y_i)+lnphiL_i-lnphiV_i`，在全部正进料组分上取最大绝对值。方向为logK逐次代入修正 `r/s`，其中 `s=max(1,max|r|/2)`。该方向不保证最大残差每步下降，故不能将其作为全程单调回溯目标。

对归一化、守恒的内部RR解，定义 `c_i=x_i*y_i/z_i`、`H=sum((y_i-x_i)^2/z_i)>0`。将汽相组分摩尔数记为 `v_i=beta*y_i`，RR隐式微分与守恒给出

```text
dbeta = sum(c_i*dlogK_i)/H
dv_i  = c_i*[dbeta+beta*(1-beta)*dlogK_i]
dG    = -sum(r_i*dv_i)
G'(0) = -[beta*(1-beta)*sum(c_i*r_i^2) + (sum(c_i*r_i))^2/H]/s
```

这里G为每摩尔进料的无量纲总Gibbs，公共组分参考项因守恒相消。最后一步由同温同压Gibbs–Duhem关系给出，不需要EOS完整Jacobian或数值差分。非退化可表示时是下降方向；H、方向或斜率不可表示时返回诊断，不宣称全局收敛。

每次alpha从1开始减半，最多24次尝试。预测下降 `pred=-alpha*G'(0)` 大于新旧两点Gibbs算术guard之和时，要求 `G_new <= G_old-1e-4*pred`（gibbs_armijo可配置）。预测下降不可分辨时，要求 `G_new<=G_old+guard_sum`，且残差达到原逸度标准，或 `norm_old-norm_new>1e-4*(alpha/s)*norm_old`。两个分支互斥，不能靠舍入后的G相等接受空步或残差增大的步。接受进展不是最终收敛；全套联合验收标准保持。

默认最多512次外层更新；计数包含失败与拒绝的单相物性调用，另记录可分辨Gibbs下降步数及暂时增残差的接受步数。输入域/逻辑/内存分配异常不吞掉；数值物性失败可在回溯中缩步恢复。资源上限不是精确字节配额，服务层仍须限制参数元数据。

## 5. 联合验收与最终相集合

| 条件 | 默认阈值与含义 |
| --- | --- |
| 逸度 | `max_active abs(log(fL/fV)) <= 1e-11`，不加摩尔分数权重。 |
| 守恒绝对值 | `max_i abs((1-beta)*x_i+beta*y_i-z_i) <= 1e-12`。 |
| 守恒相对值 | 每个正z_i用自身作分母，最大相对误差<=1e-10；无绝对下限掩盖trace。 |
| 相分率 | beta与1-beta均>1e-10；更小视为相消失未分辨，不裁剪为零再称成功。 |
| 相区分 | max abs(logK)>1e-7，且 `ZV-ZL>1e-8*max(ZV,ZL)`，排除同组成/同密度伪两相。 |

PR适配器从机械可用根选择请求的低Z/高Z候选，用原phase核计算性质；病态或未分辨根直接诊断。两种请求在单根区可能落到同一分支，不能因调用两次就算两相。密度排序不能证明任意混合物必为VLE，本基线只求给定的汽液候选分支；LLE、临界区、纯物质饱和时相分率不可确定等情形可以返回未确定。

从已满足上述要求的候选中选择Gibbs最低者，计算相对进料的Gibbs变化，再采用

```text
d_i = midpoint(log(x_i)+lnphiL_i, log(y_i)+lnphiV_i)
D_set(w) = sum w_i*[log(w_i)+lnphi_i(w)-d_i]
```

对整个允许试探相集合重新搜索，不再使用进料化学势。新入口 `test_pt_stability_against` 共用旧驱动，通过 `lnphi_equivalent=d-log(z)` 进行代数参考转换；等价参考不是真实进料相。结果存 `imposed_log_activity`，`reference` 留空，不虚构参考物性调用。原入口的表达式、下降逻辑与容差保持不变。

复核初值包括两个已求得相、默认组成搜索及调用者额外初值。有效负值阈值为用户最终TPD阈值加 `0.5*max|log(fL/fV)|`，记录为 common_reference_allowance 和 final_stability.options.tpd_tolerance；这是两个数值参考间的差异，不是到精确平衡解的严格误差界，也不改变1e-11逸度标准。算术guard亦不是EOS总误差界。

| 完整结果 | 严格含义 |
| --- | --- |
| single_phase_no_instability_found | 初始有限搜索未检出失稳，不提供全局单相证明。 |
| two_phase_no_instability_found | 不同两相候选满足联合方程标准、Gibbs未显著高于进料，且共同切平面有限搜索未检出进一步失稳。 |
| phase_set_unstable | 两相方程收敛，但共同切平面发现更低Gibbs试探分支；保留候选和证据，拒绝两相成功，不自动拼接第三相。 |
| indeterminate | 初始未确定、无可用分裂、退化、预算耗尽，或方程已收敛但最终复核未完成；后者仍保留candidate和equations_converged=true。 |

全部global_stability_proven=false。equations_converged只描述所选候选方程，不替代status。Gibbs比较和回溯guard来自明确算术幅值，不是区间认证。

## 6. 必要验证与科学边界

新增31项CTest：原27项覆盖解析RR、零/trace/极端logK、边界与预算、外部参考等价/共同切平面反例、守恒种子/SSI、异常隔离、单/两相、初始/最终未确定、第三分支失稳、相消失、PR两体系6个进料、顺序/所有权/根失败及公共头；修正增加解析Gibbs斜率、舍入接受分支、PR暂时增残差路径、公开结果访问器防护4项。共享稳定性驱动改变，原22项亦属必要回归。

理想解析夹具使用 `phiL=(2,1/4), phiV=(1,1)`，z=(.6,.4)的解为x_A=3/7、y_A=6/7、beta=.4。附加lnphi=(-1,2)分支用于验证“方程收敛但相集合进一步失稳”。Gibbs斜率测试以有理数组成及解析ln2恒等式分别验证对角项和beta耦合项，舍入接受测试明确拒绝空步、增残差和G上升。这些是制造数学模型，不是实验流体。

PR参数来自[R2]，未将其TPD驻点当成平衡组成。独立 `reference_split_decimal.py` 采用Decimal(80)、原Z三次式、直接PR逸度、二元相组成空间Newton及杠杆规则，不调用生产SSI或RR。有限差分仅用于参考Newton的导数，不是生产AD，也不是唯一正确性依据。脚本检查EOS/log逸度残差，生成六状态各五值（同温压二元相组成/Z有重复），字面值核对相对1e-35。C++ PR值回归用 `2e-9*(1+|expected|)`；联合守恒/逸度仍执行原生产阈值，普通解析比较用 `2e-12*(1+|expected|)`。

正式编译/执行只在GitHub官方GCC Debug+ASan/UBSan、Clang Release、MSVC Release配置进行，以具体提交日志为准。两个workflow按必要头/构建/测试路径选择，原稳定性工作流收窄头路径，支持依赖PR的基分支；GitHub按整个PR差异触发时仍可能在后续提交重跑已通过稳定性集成，不声称按单个提交实现缓存跳过。没有运行未受影响的旧AD或thermo套件。

```sh
cmake -S tests/flash/pt_split -B build/pt-split -DCMAKE_BUILD_TYPE=Debug
cmake --build build/pt-split --target mpmc_pt_split_tests --config Debug
ctest --test-dir build/pt-split -C Debug -R '^flash[.]pt_split[.]' --verbose --no-tests=error
python tests/flash/pt_split/reference_split_decimal.py
```

文献测试是模型数值回归，不是实验验证。没有认证全局搜索、三相求解、反应、电解质、固相、水合物、水互溶/独立水模式、SW/CPA、解灵敏度、macOS/其他架构、并发竞态或HPC性能结论。原有参数来源/适用范围检查仍保留，不自动认为普通PR适用于含水互溶。

## 7. 首轮官方失败与修正记录

首轮相分裂run `34317176470`：GCC/Clang编译通过，各22/27；MSVC遇C4702编译失败，未执行测试。共享稳定性run `34317176472` 三个作业均success。修正前审计意见已写入PR #10，comment `5596683428`。

生产缺陷是强制最大逸度残差单调下降。独立binary64算术检查N2进料(.3,.7)、试探(.5,.5)的守恒初值：一步残差约.03699765→.07523118，总G约-1.19811566→-1.20034184，确实是Gibbs下降却被旧规则拒绝。由第4节的RR隐式微分和Gibbs–Duhem推导修正步长接受，不放宽1e-11逸度、1e-12绝对/1e-10相对守恒、相区分或分率标准；原六个进料、预算与30个参考值不变。独立Python诊断不冒充生产C++正式测试。

另外两项属于测试契约/可移植性问题：K=1的trace在返回的显式舍入归一化后不应要求逐位等于输入，修正为零严格为零、正trace保持正、返回相对误差<=128eps、乘回已记录原始相组成和后的相对误差<=8eps；没有以绝对误差下限隐藏trace丢失，原守恒断言不变。MSVC将两个必抛异常的静态lambda模板实例识别为不可返回，导致调用后的验证被报C4702；改用测试层std::function动态故障回调，仍注入同样异常，不改生产编译选项或关闭告警。

同时补充可变结果对象的索引/空point/状态检查，防止调用者改写公开记录后在candidate()解引用无效存储。新增回归保留原输入和期望，修正通过与否以随后官方日志为准。

## 8. 原始依据

[R1] Gernert, J.; Jäger, A.; Span, R. (2014), Calculation of phase equilibria for multi-component mixtures using highly accurate Helmholtz energy equations of state, Fluid Phase Equilibria 375, 209–218. DOI `10.1016/j.fluid.2014.05.012`。核验式(11)–(14)、(19)–(22)、(27)–(30)。[论文](https://www.sintef.no/globalassets/project/impacts/phase-equilibria-gernert-et-al-fpe-2014.pdf)。未声称复现全文算法；本轮Gibbs回溯方向导数单独推导。

[R2] Hua, J. Z.; Brennecke, J. F.; Stadtherr, M. A., Reliable Computation of Phase Stability Using Interval Analysis: Cubic Equation of State Models，1997年10月修订作者稿，第4.3/4.4节、表4/5；甲烷参数见第4.1节。[作者稿](https://academicweb.nd.edu/~markst/zm97a.pdf)。正文/表4末个N2进料.65/.60的差异保留于PR #9审计，本次平衡参考用明确的.18/.30/.44；CO2用.20/.30/.43。未复制区间算法、实验数据或全文。
