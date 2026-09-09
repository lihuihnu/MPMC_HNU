# PT 闪蒸基础：相稳定性搜索

本增量提供给定 `p,T,z` 的非反应、同温同压、流体候选相 **TPD 相稳定性搜索**，以及 PR76 适配器。它不计算平衡相组成、相分率或最终相数，不是完整闪蒸器。返回的组成是试探点／失稳证据，不是闪蒸结果。

## 1. 编写前交叉审计

基线为 `main@c48ac3a3534e93e005beb2db29fdd293676cfdcc`，tree `d2639183af265feee11c676251f982d11de54af9`。已读取根 AGENTS、README、参数／混合／PT／根接口、相关构建、历史及开放事项；当时没有开放 Issue/PR，父目录没有附加 AGENTS。以下是单一实现者从多方面与多份原始研究交叉核验，不冒称有外部独立审稿人。

| 方面 | 风险与取舍 |
| --- | --- |
| 热力学 | 相同进料组成有多个根时，参考态必须比较 Gibbs 能；负 TPD 不要求搜索先收敛。机械稳定不等于组成稳定。 |
| 算法 | 多初值不能证明找到全部极小点；初值收敛到平凡点或鞍点不能证明全局稳定。输出不提供 `stable` 成功枚举。 |
| 数值 | 零进料、迹量、根退化、非光滑根切换、累计舍入与迭代预算分别处理；不能靠小组成更新量退出。 |
| 架构 | 新建 `flash` 调用既有热力学，不反向依赖、不复制 PR 公式。保留 PR76/AD 的源码、接口、编译设置及原测试不变。 |
| 验证 | 理想混合／正则溶液解析判据、独立 Gibbs 表达式、AD 切向导数、原始文献数值问题与高精度独立参考相互交叉。 |

**来源限制与路线修订：** Michelsen (1982), DOI `10.1016/0378-3812(82)85001-2` 的全文镜像本轮未成功读取；不声称核验或逐式复现原版。已直接核验 Gernert 等 [R1] 期刊第 211–212 页式 (21)–(24)，并用 Hua 等 [R2] 第 3–4 页的切平面条件／参考根选择及第 10–11、21–22 页的 PR 问题交叉验证。本实现是带下降保护的归一化逐次代入搜索，不是原版 Michelsen 算法，也不是 Hua 的区间全局算法。

比较方案：只有 TPD 单点函数不足以完成搜索；直接固定点迭代没有下降保护；一次引入 Hessian、Newton/BFGS、通用优化器或区间证明会扩大依赖与审计范围。因此先选可推导下降方向的对数组成步＋回溯，保留用户补充初值与清晰失败诊断。**生产搜索只支持 double**，不实现搜索过程 AD 或稳定性解灵敏度。

## 2. 接口与结果契约

`pt_stability.hpp` 仅依赖标准库。`pr76_stability.hpp` 依赖既有热力学头，不依赖 AD。CMake 目标 `mpmc::flash` 引入 `mpmc::thermodynamics`；根 AD 构建入口不变。前端、网络、MPI、PETSc、SW/CPA 或注册框架都不是前置条件。

```cpp
#include <mpmc/flash/pr76_stability.hpp>
// parameters 为调用方已建立、来源可追溯的 PrParameterSet；不内置参数数据库。
const auto phase = mpmc::thermodynamics::Pr76Phase<double>::from_parameters(parameters);
mpmc::flash::Pr76StabilityEvaluator evaluator(phase);
const auto result = mpmc::flash::test_pr76_pt_stability(p_pa, t_kelvin, z, evaluator);
// result.search.status；不要把 lowest_sampled->composition 当作平衡组成。
```

通用入口 `test_pt_stability(p,T,z,provider,options={},extra_starts={})` 接受可重复、同步的 `provider(p,T,w)->StabilityPhase`。提供者须保持组分顺序、共同化学势参考态、Gibbs–Duhem 一致性，并选取该组成下允许流体分支中最低 Gibbs 能分支。接口检查数组维数与有限性不能证明这些热力学前提成立；不能将任意不一致的 ln(phi) 回调称为物理模型。

| 状态 | 严格含义 |
| --- | --- |
| `unstable` | 至少一个可行试探点的负 TPD 越过声明阈值与舍入保护量；保留该证据。即使其他初值失败也不能抹掉已有证据。 |
| `no_instability_found` | 本次全部要求初值达到驻点残差要求，但未找到越阈值的负值。**不是全局稳定性证明，也不保证驻点是极小值。** |
| `indeterminate` | 无足够负值证据，且至少一个搜索未完成／非光滑／物性失败／预算耗尽。调用方不得把它当单相成功。 |

`global_stability_proven` 固定为 false；`lowest_sampled` 仅为最低已采样值，不是已求得的全局最小值。每条 trial 保留原始初值、终点、TPD、未加权驻点残差、分支诊断、迭代数、物性调用数、回溯次数、拒绝的物性调用数、最近物性问题和终止原因。成功恢复后仍可能保留历史拒绝问题，终止状态与历史问题不可混淆。

结果拥有自身存储，并记录 p/T、原进料和、所用归一化进料与 options。PR 包装另存 dataset_id/revision、有序组分 ID、模型 profile、PT 数值约定和根选项；原参数来源由拥有快照的 evaluator.model().parameters() 查询。evaluator 可顺序复用，禁止并发共用；独立实例无共享可变状态。回调不得保留临时 span。未进行并发竞态验证或性能基准。

## 3. 输入、零组分与搜索范围

p 为 Pa，T 为 K，均有限且正；组成数组非空、有限、位于 [0,1]，其和只接受 `64*epsilon(double)` 内的舍入偏差。驱动层**显式**除以已记录的进料和，并同样对搜索初值／对数更新归一化；不改变原 PT 核“仅校验、不归一化”的契约。不接受任意未归一化数据，不裁剪负数。单点 `tangent_plane_distance` 只校验已归一化输入，不替调用者修改它。

活跃集合为 `z_i>0`，不设工程 cutoff。`z_i=0` 的组分在所有试探相中保持严格零（封闭非反应物料约束）；不能在该支持集之外添加组分。搜索初值在每个活跃组分上必须严格正，否则拒绝。单点 TPD 允许活跃试探分数为零，用 `lim(w log w)=0` 计算值，但该边界点的驻点残差为无穷，不能用内部光滑收敛条件接受。

默认初值集合：进料、活跃组分均匀分配，以及每个活跃组分的 `0.1*z+0.9*e_i` 富集点。单活跃组分只用进料；其他完全重复的初值保留而不伪称独立信息。调用者可补充初值，或关闭默认初值后提供非空集合。此策略是明确的启发式，不保证覆盖所有吸引域；不使用未核验 Wilson 修正或隐藏随机数。

对数更新使用移位 log-sum-exp。若活跃分量在 exp／归一化或初值生成时丢失为零，拒绝该步或报告不能表示，**不填 epsilon 下限**。资源上限先于初值分配和回调检查：默认最多 256 组分、1024 初值、262144 个初值元素及 100000 次物性调用。实际内存还包括终点／残差、模型快照与回调存储，不把元素上限说成精确字节配额。服务层仍须限制输入文件和来源字符串。

## 4. 数学与停止判据

在活跃支持集上，定义

```text
q_i = log(w_i/z_i) + ln(phi_i(w)) - ln(phi_i(z))
D   = sum_i w_i*q_i
r_i = q_i - D                    （精确 sum(w)=1 时）
```

Gernert 等式 (21)–(22) 给出归一化 TPD。固定 p/T、光滑且热力学一致的分支满足 `sum_i w_i*d ln(phi_i)=0`，因此在单纯形切向上 `dD=sum_i q_i*dw_i`。这给出一阶梯度，**不是有限差分，也不需要组装 EOS 完整 Jacobian**。测试以已有 AD 对整个 TPD 表达式求切向导数交叉核验该消元。

令 `s=max(1,max|r_i|/log_step_limit)`，`delta_i=-r_i/s`，采用

```text
w_i(alpha) = exp(log(w_i)+alpha*delta_i) / sum_j exp(log(w_j)+alpha*delta_j)
D'(0) = -sum_i w_i*r_i^2/s <= 0
```

若无步长缩放且 alpha=1，该表达式归一化后即 [R1] 式 (24) 的逐次代入更新；回溯与截断**步长**是本实现独立推导的保护，不截断组成或化学势。浮点运算中心使用 `D/sum(w)`，避免归一化尾数使残差均值有微小偏差。

默认 `max|r_i|<=1e-8` 为驻点收敛条件，绝不使用 `max|w_i*r_i|` 或仅使用组成步幅。否则最多 512 步，每步最多 32 次尝试，从 alpha=1 每次减半，以 `D_new<=D+1e-4*alpha*D'(0)` 接受。

在预测下降已不大于两点舍入保护量之和时，**改用残差进展规则替代 Armijo**，允许 `D_new<=D+guard_sum` 且驻点残差严格下降至少 10% 的步，避免正常舍入导致停滞；该规则**不是收敛判据**，最终仍必须通过原残差阈值。物性失败可拒绝后缩步恢复。预算统计包含参考态、拒绝步与失败调用；资源耗尽不伪造收敛。

采用补偿求和，近似相同组成用 log1p 计算比值。负证据要求

```text
D < -tpd_tolerance - guard
默认 tpd_tolerance = 1e-10
 guard = 256*epsilon(double) *
         [1 + sum_i w_i*(|log(w_i)|+|log(z_i)|+|lnphi_i(w)|+|lnphi_i(z)|)]
```

guard 是公开的算术保护尺度，**不是 EOS 总误差的严格上界，不是区间认证**。不能以任意稀释 trace 或抬高容差隐藏不收敛。任何已求得的可靠负点立即结束该条搜索，不要求该证据同时是驻点；其他要求的初值继续执行并保留结果。没有相分裂迭代或跨相导数。

## 5. PR76 适配与失败隔离

从既有 `roots_full` 枚举 `Z>B` 根，只将 `slope_sign>0`（等价于 dp/dv<0）的根用于流体 Gibbs 比较；中间机械不稳定根不被称为稳定相。以补偿的 `sum_i w_i*(lnphi_i(candidate)-lnphi_i(best))` 比较，选择最低者，不固定选最大／最小根。

根集合未分辨、根迭代失败、非有限比较或所选根 `derivative_valid=false` 时，保守地返回物性诊断。虽然本搜索不调用根 AD，仍使用既有病态根标志限制不可靠的搜索区域。两个允许根 Gibbs 差落在 `256*eps*[1+sum w*(|phi_a|+|phi_b|)]` 以内时，所选值仍取实际较小者，但 smooth=false，不能按光滑驻点接受。该 tie 阈值也不是严格认证。

只捕获明确的 `StabilityPropertyError` 为搜索失败。PR 适配器将数值根／range 异常映射为该类型；温压数据越界、维数／参数错误仍抛给调用者。回调维数错误为 logic_error；内存分配失败及其他程序错误不吞掉。返回 reference 失败时不发布半套参考逸度。

## 6. 必要验证与科学边界

独立 CTest 工程目前包含 22 个条目：理想 TPD／搜索、正则溶液 TPD／搜索、平凡鞍点反例、零／迹量、非法输入、预算、故障分类、回溯恢复、舍入区下降保护、下溢保护、非光滑、所有权恢复、最低 Gibbs 根、PR 的 AD 梯度交叉核验、PR 根失败、24 种排列与变组分数、两个文献问题、独立高精度锚点、独立公共头。正式执行结果以具体提交的 GitHub Actions 日志为准，源文件存在不等于通过。

普通固定数值比较用 `2e-13 + 2e-11*|expected|` 无量纲预算；完整 TPD 的 AD 与 Gibbs–Duhem 切向导数比较用 `5e-12 + 2e-10*|expected|`，覆盖独立根、逸度和消元路径的舍入。收敛／负证据检查使用生产声明阈值，不用这些参考比较预算替代它们。参考脚本以 Decimal(80) 原 Z 三次式 Newton 修正、完整双重混合和直接逸度式生成 3 个固定 TPD 锚点，文本核对相对误差 1e-35；不导入生产实现，不用生产输出拟合期望值。

文献数值回归来自 [R2] 的 PR 问题：氮气／乙烷 270K、76bar；CO2／甲烷 220K、60.8bar。纯组分参数和显式 kij 均带文献来源与 bar→Pa 换算记录，仅位于测试夹具。各选五个进料检查文献稳定情况对应“未检出”、文献失稳情况对应负证据；另对两个公开试探组成直接核验负 TPD。**[R2] 正文给氮气／乙烷最后进料 .65，表 4 给 .60；本测试明确采用表 4 的 .60。** 文献表格有效位数有限，不拿它们作高精度期望值；保留 PR76 印刷系数、精确根式和本项目现代 SI R，未拟合论文未知的历史实现常数。

这些是模型间数值回归和代数／软件验证，**不是实验数据验证，也不是移植区间算法的全局可靠性保证**。没有认证全局搜索、反应、电解质、固相、水合物、含水互溶、SW/CPA、相分率、macOS／其他架构、性能或并发竞态结论。化学组分／相态／水处理能力须由上层有依据地选择，不因本接口接受数组就视为有效。

新增 workflow 仅构建／执行稳定性集成；使用已有官方 Ubuntu GCC Debug＋ASan/UBSan、Clang Release、Windows MSVC Release 配置及已审计 checkout 固定提交。未来 thermo/AD 头或构建改变时触发本下游；本次不改它们，所以不重跑无关旧套件。文档更新不另开全量测试，未修改许可证、保护规则或设置。

```sh
cmake -S tests/flash/stability -B build/flash-stability -DCMAKE_BUILD_TYPE=Debug
cmake --build build/flash-stability --target mpmc_flash_stability_tests --config Debug
ctest --test-dir build/flash-stability -C Debug -R '^flash[.]stability[.]' --verbose --no-tests=error
python tests/flash/stability/reference_decimal.py
```

### 首轮官方验证发现的回溯缺陷及修正

首轮 run `34307917973` 在 GCC/Clang 通过全部 21 项；MSVC 编译通过、20/21 项通过，失败项为 `backtracking_recovery`，两组文献回归均通过。复核发现原逻辑把舍入区残差进展规则与普通 Armijo 用 OR 连接，因此后者仍可凭舍入后的目标值相等或假下降接受残差增大的步，绕过保护。这是生产步长接受逻辑缺陷，不把它归咎于测试容差。

独立 binary64 算术检查给出一个可说明该缺陷的正则溶液例：interaction=-10，z=(0.5,0.5)，w0=0.5+91*2^-36；半步会将未加权残差由约 1.5891e-8 增至 3.1781e-8，但目标值的相消残留可从约 3.8579e-17 变成约 -1.2219e-17，原 Armijo 仍会接受。精确 Gibbs 函数在此附近严格凸，残差放大对应错误的过大步长。这是独立数学诊断，不冒称取得了 MSVC 内部迭代轨迹。

修正为互斥分支：下降可分辨时用 Armijo；下降不可分辨时仅允许目标值在舍入保护范围内且残差下降至少 10% 的步。保留原收敛／负值阈值、预算、原失败输入和全部原测试，新增 `roundoff_descent` 检查不可接受空步／增残差步、正常下降和四个近驻点凸模型搜索。回溯恢复测试增加状态与残差日志。残差进展要求在非常平坦区域可能拒绝缓慢步，届时返回 indeterminate，不能为追求收敛报告放宽阈值。

## 7. 来源

- [R1] Gernert, J.; Jäger, A.; Span, R. (2014), *Calculation of phase equilibria for multi-component mixtures using highly accurate Helmholtz energy equations of state*, Fluid Phase Equilibria 375, 209–218. DOI: `10.1016/j.fluid.2014.05.012`. [SINTEF 提供的论文](https://www.sintef.no/globalassets/project/impacts/phase-equilibria-gernert-et-al-fpe-2014.pdf)，本轮核验 PDF 第 3–4 页式 (21)–(24)。不把其正试探结果的启发式稳定假设当成数学证明。
- [R2] Hua, J. Z.; Brennecke, J. F.; Stadtherr, M. A., *Reliable Computation of Phase Stability Using Interval Analysis: Cubic Equation of State Models*, 作者稿 April 1997, revised October 1997. [作者主页 PDF](https://academicweb.nd.edu/~markst/zm97a.pdf)，核验 PDF 第 5–6 页理论与参考根、12–13 页参数、23–24 页表 4–5；本文数值回归采用该版本，不声称已核验最终出版版或全部勘误。不复制其区间算法或源码。

**下一步建议：** 审阅本增量后，实现仅限已获明确失稳证据体系的汽液相分裂基线，复用 PR 物性，检查 Rachford–Rice 守恒与逸度残差，并对最终相集合复核稳定性／保留未确定状态；不自动把本次“未检出”升级为认证单相。该后续不属于本次实现授权。
