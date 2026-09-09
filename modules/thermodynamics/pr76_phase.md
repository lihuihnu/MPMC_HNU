# PR76：PT 候选相性质核

本增量计算 `p,T,w -> Z,ln(phi_i)`，用于后续 PT 闪蒸的试探相求值。**它不是闪蒸器，不决定稳定相、相数、相分率或各相平衡组成。** 三个 EOS 根不是三个共存相；负斜率条件也不能替代组成空间的相稳定性搜索。普通浮点路径仅依赖标准库和已有热力学头，AD 由调用者包含 `<mpmc/ad/math.hpp>`。

## 1. 来源、数值约定与审计取舍

Peng & Robinson (1976), *A New Two-Constant Equation of State*, DOI [10.1021/i160057a011](https://doi.org/10.1021/i160057a011)。本轮对照[作者提供的原文 PDF](https://www.researchgate.net/profile/Ding-Yu-Peng-2/publication/231293953_New_Two-Constant_Equation_of_State/links/5b50a161aca27217ffa634b5/New-Two-Constant-Equation-of-State.pdf)第 2 页／期刊第 60 页的清晰页面图像核验式 (4)–(8)、(19)–(22)。此前用户提供 PDF 的摘要保留在 [纯组分说明](pr76_pure.md)；本轮不把在线页面核验称作重新校验该文件摘要，不重新分发 PDF、截图或全文。

保留 `PR76/classical-vdw/constant-kij/no-translation`，沿用原纯组分印刷系数和现代 SI `R=8.31446261815324 J/(mol K)`。**原文式 (19) 的对数项印为 2.414、0.414，本核明确采用由式 (4) 分母分解得到的精确 1±sqrt(2)**，保证逸度与所用 EOS 代数一致；这不是声称原文印有更多有效数字。实现约定为 `pr76_pt_convention = PR76/printed-coefficients/R-SI-2019/exact-sqrt2/PT-v1`，不改变参数数据集的身份。

使用 [IDAES 的通用立方 EOS 官方推导](https://idaes-pse.readthedocs.io/en/stable/explanations/components/property_package/general/eos/cubic.html)交叉核对精确根式和逐组分混合项；不复制其代码或参数，不引入第三方库。实现独立编写。

原 `Pr76Mixture::evaluate_full/evaluate_reduced`、返回类型与 a/b 累加顺序不变。新增私有、仅供 `Pr76Phase` 使用的编译期路径，在同一次纯组分／上三角混合遍历中计算 `s_i=sum_j w_j*a_ij`；不复制混合公式，不为获取行和而组装完整 Jacobian，不暴露容易错配的可变参数视图。该变化须运行原 15 项混合核回归。旧 AD、数据契约、纯组分源码及其测试／构建入口不修改。

专用三次方程仅用于 PR76，留在 `thermodynamics`，本次不创建空 `numerics`、`flash` 或动态插件框架。相性质不调用闪蒸，后续闪蒸使用相性质而不重复 EOS 公式。

## 2. 输入、坐标与结果

`Pr76Phase<T>::from_parameters(parameters, limits={})` 拥有一个已校验、有序的参数快照，T 为 float/double/long double。参数固定，不对 Tc/Pc/omega/kij 求导。重排、增减或替换组分须创建新的快照／模型；数值数组顺序始终与模型一致，同长度的错误 ID 顺序不能由裸数组识别。

| 接口 | 契约 |
| --- | --- |
| `roots_full(p,T,w,workspace,options={})` | 普通浮点数，n 个完整摩尔分数；返回所有已分辨的 `Z>B` 根及状态。 |
| `roots_reduced(p,T,w_independent,workspace,options={})` | 普通浮点数，前 n−1 个分数；最后一个为 1−sum；n=1 时为空。 |
| `evaluate_full(p,T,w,root_index,workspace,options={})` | 普通浮点或 `Dual<T,K>`；在本次状态重新求根并求所选根的 ln(phi)。 |
| `evaluate_reduced(...)` | 同上，使用约化组成。 |

p 单位 Pa、T 单位 K，均须有限且 >0；所有 AD 种子须有限。已声明的温压边界均按闭区间检查，未知仍为未知。完整／约化组成校验及导数语义继承 [混合核](pr76_mixture.md)：完整输入允许无量纲和误差 `64*epsilon(T)`，不归一化、不裁剪；完整种子独立传播，约化种子通过最后分数传播负的种子和。非约束完整偏导只是指定表达式的形式导数；约束偏导为 full_i−full_last，不能混称。

`root_index` 是本次可用根的递增 Z 下标，不是固定的相标识。任何状态变化后都须重新检查拓扑；不承诺跨根消失、出现或合并连续跟踪旧下标。中间根可供诊断，不自动标为物理稳定相。

结果 `Pr76PhaseValues<Number,T>` 包含 `z`、拥有存储的 `ln_phi` 向量、普通浮点诊断 `root_set` 与所选下标。返回的数值和导数均须有限；异常不发布部分结果。`root_set.count` 仅在 status=success 时有结果含义。每根保留 Z、独立计算的 X=Z−B、缩放多项式残差、原有理 EOS 残差、斜率比、导数有效标志、斜率符号和迭代次数。这些是数值／机械诊断，不是相稳定性结论。

## 3. 根计算与精度边界

由原文式 (5)–(7)：

```text
A = a*p/(R*T)^2; B = b*p/(R*T)
F(Z) = Z^3 + (B-1)Z^2 + (A-3B^2-2B)Z + B^3+B^2-A*B
```

用 `X=Z-B>0` 排除不可用根，再用 `s=max(1,B,sqrt(abs(A)))`、`y=X/s` 缩放：

```text
H(y) = y^3 + c2*y^2 + c1*y + c0
c2 = (4B-1)/s
c1 = (A+2B^2-4B)/s^2
c0 = -2B^2/s^3
```

按 H' 的实驻点分割正半轴，在 Cauchy 根界内对异号的单调区间二分。稳定求 H' 的二次根，使用 FMA 求多项式；不以直接 Cardano 在多根／小根处可能发生的抵消作为生产基线。保存 X 避免重新用 Z−B 相减丢失自由体积。

| 数值规则 | 含义 |
| --- | --- |
| 驻点处 `abs(H)<=64*eps*sum(abs(terms))` | 浮点符号不足以可靠区分拓扑，返回 near_multiple，不猜测重数或给出假根。 |
| 二分区间宽度 `<=4*eps*max(y,min_normal)`，或相邻端点／精确零 | 停止迭代后仍必须通过下面的残差检查。 |
| 多项式缩放残差 `<=128*eps` | 向后误差判据，不是临界区根的前向误差保证。 |
| `abs(1/X-A/(Z^2+2BZ-B^2)-1)/(1+abs(1/X)+abs(A/(Z^2+2BZ-B^2)))<=512*eps` | 回代原始有理 EOS，防止仅有三次式小残差。 |
| 斜率比 `abs(H')/(3y^2+2abs(c2)y+abs(c1)) > 8*sqrt(eps)` | 本实现允许局部导数的保守数值门槛，不是普遍条件数上界。 |

所有上述量无量纲，eps 为所选 T 的机器精度；不是为了拟合实验而设置的科学容差。`Pr76RootOptions.max_iterations` 默认每区间 2048，只能改变工作量上限，不能放宽精度。达到上限返回 iteration_limit。缩放常数、根或非零端点无法表示，或回代不合格，返回 unrepresentable。不承诺全部浮点动态范围可用；极低 B 导致 c0 下溢时拒绝整个拓扑，而不是只返回一个未全面检查的气体样根。

允许有限的零／负 A 作为已有混合核的代数延拓，不裁剪为正，也不声称这种参数或状态已物理验证。多组分仍继承各纯组分 a_i>0 的主平方根限制。行和需要比加权总 a 更强的可表示性，无法表示 s_i 或其导数时明确失败。

## 4. 稳定逸度表达式与局部隐式求导

令 `s_i=sum_j w_j*a_ij`，`r_i=b_i/b`。式 (19) 写为

```text
ln(phi_i) = r_i*(Z-1) - ln(Z-B)
            - A/(2*sqrt(2)*B)*(2*s_i/a-r_i)
              * ln((Z+(1+sqrt(2))*B)/(Z+(1-sqrt(2))*B))
```

生产路径先符号约去 a 和 B：设 `d=X+(2-sqrt(2))*B`，`u=2*sqrt(2)*B/d`，`L=log1p(u)/u`，则吸引项为

```text
[(2*s_i-a*r_i)*p/(R*T)^2] * L/d
```

因此 a=0 不需要除零或 epsilon 正则化。各 `ln(phi_i)` 不使用 `log(w_i)`，允许零／迹量组分，不删除它们的加入导数。

当 |u|<=1/8 时，L 用交错 Taylor 多项式，阶数 `N=digits(T)/3+2`；省略的值项不超过 `|u|^(N+1)/[(N+2)(1-|u|)]`，斜率项不超过 `|u|^N/(1-|u|)^2`，均低于 eps。其余区域用 `log1p(u)/u`，避免在 u≈0 时做带抵消的导数商法。分段是有误差界的数值求值策略，不是修改 EOS。

在 |Z−1|<1/4 且 B<1/8 的气体样区域，使用 EOS 恒等式

```text
Z-1 = B/X - (A/Z)/(1+2B/Z-(B/Z)^2)
ln(X) = log1p((Z-1)-B)
```

保留 Z 已舍入为 1 时仍非零的 O(p) 逸度项；其余区域使用保存的 X 的 log。两个分支表示同一数学函数，均保留 AD 依赖。

根先用普通浮点数收敛，再在固定简单根上应用隐函数定理：`dy/dtheta=-H_theta/H_y`。缩放 s 在该局部线性化中固定，`dZ=s*dy+dB`。实现让 AD 求 H 对输入的偏导，只移除 H 的原值以保持收敛 Z 不变；不对二分分支／迭代次数求导，不使用有限差分。该局部修正只支持已有一阶 Dual，不声称高阶或嵌套 AD 正确。

若所选根 `derivative_valid=false`，普通浮点仍可读取通过残差检查的原值，但 AD 求值抛出 `ill_conditioned_derivative`，包括所有种子为零的 Dual。导数只在固定拓扑、固定根的局部光滑域有效，不是跨相边界的全局导数。

## 5. 错误、所有权与使用

根求解返回 `success/near_multiple/iteration_limit/unrepresentable`。相求值将后三者转换为带 `code()` 的 `Pr76PhaseError`；导数拒绝也用此类型。非法根下标为 out_of_range，非法维数／求根选项为 invalid_argument，非法温压、组成和种子为 domain_error，无法表示的 A/B、行和、逸度或导数为 range_error；底层参数／数学／分配异常保留传播。

模型可拷贝／移动构造、不可赋值，移出后的模型不可继续计算。参数引用仅从左值取得。工作区不可拷贝／移动，仅保存每次覆盖的临时量，不跨 Jacobian 批次缓存 AD 值；同一个工作区不能并发使用。失败后内容未指定，但下一次调用重新覆盖，已经返回的结果不受影响。

准备模型为 O(n²) 存储，普通相求值混合项 O(n²)，额外行和与结果 O(n)。根迭代有界；调用返回向量会分配内存，工作区容量不足也会分配。未进行性能测量，不称为已经达到 HPC 性能目标。

```cpp
#include <mpmc/thermodynamics/pr76_phase.hpp>
// parameters is an audited PrParameterSet; w follows its ordered component IDs.
namespace th = mpmc::thermodynamics;
const auto model = th::Pr76Phase<double>::from_parameters(parameters);
th::Pr76PhaseWorkspace<double> work;
const auto candidates = model.roots_full(p_pa, t_k, w, work);
if (candidates.status != th::Pr76RootStatus::success) {
    // Report the diagnostic; do not use count or fabricate a phase.
    throw std::runtime_error("Unresolved PR76 roots");
}
const auto candidate = model.evaluate_full(p_pa, t_k, w, candidates.count - 1, work);
// Largest Z candidate only, NOT a proven stable vapor or a completed flash.
```

这段代码依赖调用方提供合法参数、组成和温压，不包含真实物性数据；完整可执行人工构造见测试。

## 6. 必要验证与下一步

新增独立 `tests/thermodynamics/pr76_pt`，18 个 CTest：可因式分解的精确根、拓扑／重根、尺度、相性质参考值、完整导数、约化导数、Gibbs–Duhem、纯组分退化、低压极限、log1p 比值、零／迹量组成、4 组分 24 种排列、运行期 1→4→2→3→1、输入域、所有权／失败恢复、零／负吸引项、病态导数拒绝、独立公共头。其中前 16 项各测试 float/double/long double；病态导数构造及公共头为 double。排列和标量类型属于子情形，不重复计入 CTest 数。

独立判据包括原 Z 三次式的 Cardano／三角法（仅在远离退化的参考点使用）、原始参数的完整双重混合求和、原 Z 方程隐式导数与直接式 (19) 的手工解析导数；不复用生产 shifted solver 或生产预计算系数来构造期望值。另由 Python 标准库 Decimal(80)、原 Z Newton 和直接对数式再生成 12 个数值与 48 个约化导数锚点（直接对约化坐标应用解析链式法则），脚本仅读取测试中的锚点核对，不读取生产实现、不修改参考值。全部参数明确为 synthetic_test，不是实际流体或实验数据。

普通非相消数值比较采用 `4096*eps(T)` 相对容差、零绝对容差；可精确表示的结构零单独检查，Gibbs–Duhem 的零值采用 `8192*eps*sum(abs(weighted_terms))` 缩放。此预算覆盖三次求根、对数、混合和独立导数路径的累计舍入，不替代生产的严格回代判据。低压测试由解析第一维里系数核验非零逸度及导数，不用固定绝对误差吞掉小量。数值锚点文本与 Decimal 计算以 1e-32 相对误差核对。

首次官方运行 `34302464615` 暴露了测试判据问题：在 `w=(1,0,0)`，直接微分式 (19) 后，现存组分的每个完整组成导数严格为零。因为 `dg*e+g*de=0`，其余 Z/B 导数项由 EOS 恒等式消去；当第二组分为 eps 时，该导数仅为 O(eps)，而被相减的各项仍为 O(1)。初始测试错误地用独立参考计算的舍入残留（约 4e-20）当相对误差尺度，与 float 的约 1.5e-8 舍入残差比较。修正为**相消残差判据**：逐项误差除以独立解析式六个加项的绝对值之和，系数仍为 `4096*eps(T)`，无任意绝对底限；纯组分顶点还独立核验解析参考的零恒等式。所有原输入、导数条目及溶质非零加入导数的逐项相对检查保留，生产代码、其容差及其他测试容差未改。此处不声称可以给接近零的导数提供高相对精度；报告的是有明确尺度的绝对舍入误差。MSVC 同轮还发现测试常量 kij 遮蔽模板内局部名，已仅将测试常量改名为 fixture_kij，不抑制告警。

第二次官方运行 `34303031144` 中 GCC/Clang 全部通过；MSVC 编译成功但一项约化导数参考比较失败。MSVC 的 [long double 与 double 使用相同表示](https://learn.microsoft.com/en-us/cpp/cpp/fundamental-types-cpp)，原参考用该精度完成直接式解析导数后再相减，不能假定其有额外精度。对中间根、组分 1、独立 w0，原生产值为 −0.004773184757624893，原参考为 −0.004773184757618676；独立 Decimal(80) 原 Z 隐式导数与直接式 (19) 链式法则给出 −0.004773184757622566848374166371093646，生产值相对误差约 4.87e−13，满足原 `4096*eps(double)` 约 9.09e−13 的判据。故将该固定状态的全部 48 个约化导数期望值改为由独立 80 位解析路径再生成的锚点，避免按平台精度生成不同质量的参考；保留全部原输入、根、输出及导数检查，**原相对容差和零绝对容差不变，生产实现不变**。脚本在 CI 中逐项核验 12 个原数值锚点及 48 个导数锚点，共 60 项；不从被测程序输出拟合期望值。

```sh
cmake -S tests/thermodynamics/pr76_pt -B build/thermo-pt -DCMAKE_BUILD_TYPE=Debug
cmake --build build/thermo-pt --target mpmc_thermodynamics_pr76_pt_tests --config Debug
ctest --test-dir build/thermo-pt -C Debug -R '^thermo[.]pr76_pt[.]' --verbose --no-tests=error
python tests/thermodynamics/pr76_pt/reference_decimal.py
```

正式编译和执行仅用 GitHub 官方托管 runner。当前增量运行 PT 与受影响的原混合核套件，以及工作流的依赖选择回归；不配置／构建／运行无关旧 AD、契约或纯组分套件。保留现有 GCC Debug + ASan/UBSan、Clang Release、MSVC Release 环境与告警。未来修改 pure 或 AD 时必须纳入 PT 下游，公共参数变化纳入全部相关热力学套件。结果以具体提交的 Actions 日志为准；配置和测试源码不等于测试已通过。

本轮不验证真实流体精度、相平衡、SW/CPA、水互溶、macOS／其他架构、并发竞态或性能。不修改项目许可、分支保护或其他仓库设置。

**下一步：实现供 PT 闪蒸使用的相稳定性测试。** 先核验 Michelsen TPD 原文、相分支及搜索／失败契约，复用本核的 ln(phi) 和局部导数，验证稳定、失稳及未确定情况；有限初值的局部搜索不能冒充全局证明。随后才实现相分裂与守恒／平衡复核。用户已将主线改为 PT 闪蒸，本建议替代早期混合核文档中的独立 p(T,v,w) 建议，不自动开发后续模块。
