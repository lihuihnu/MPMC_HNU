# PR76 纯组分 a(T)、b：原文核验与最小数值核

本增量只计算单个组分的吸引参数 `a(T)` 和共体积 `b`，不计算混合物参数、压力、压缩因子、逸度、EOS 根或闪蒸。输入来自已有 `PrParameterSet`；公共契约、全部 AD 实现与旧测试保持不变。新头 `pr76_pure.hpp` 只依赖数据契约和标准库，普通浮点路径不要求 AD；需要 AD 时由调用方包含 `mpmc/ad/math.hpp`。

## 1. 本次实际核验的原始材料

用户提供 `New_Two-Constant_Equation_of_State.pdf`，6 页，对应期刊页 59–64：Peng, D.-Y.; Robinson, D. B. (1976), *A New Two-Constant Equation of State*, **Ind. Eng. Chem., Fundam. 15(1), 59–64**, DOI `10.1021/i160057a011`。

```text
SHA256: 0ea7cb1a8c6206f9e929d7d7730399bb6eecf310259cd7146d9eaadb4b4a6036
```

已对照 PDF 渲染图而非仅靠文本提取核验：第 2 页（期刊第 60 页）的式 (9)、(10)、(12)、(13)、(17)、(18)，以及第 6 页的 R、T、v 等符号定义。首尾页包含相邻文章片段，不归入本模型。只记录书目信息、公式定位与文件摘要；未将 PDF/截图/论文全文或参数表提交仓库，上传并不被解释为获准公开再分发。

| 原文式号 | 核验内容 |
| --- | --- |
| (9) | `a_c = a(T_c) = 0.45724 R^2 T_c^2 / P_c` |
| (10) | `b(T_c) = 0.07780 R T_c / P_c` |
| (12) | `a(T) = a(T_c) alpha(T_r, omega)` |
| (13) | `b(T) = b(T_c)` |
| (17) | `alpha^(1/2) = 1 + kappa (1 - T_r^(1/2))` |
| (18) | `kappa = 0.37464 + 1.54226 omega - 0.26992 omega^2` |

保留原文印刷十进制系数；不替换为由临界条件重新求得的更多位系数，不使用 PR78 的高偏心因子分段修正、PRSV、体积平移或 SW 水 alpha。`pr76_profile` 仍为原数据配置；新增数值约定标识 `pr76_pure_convention = PR76/printed-coefficients/R-SI-2019` 与数据集版本分别记录。

原文只把 R 定义为气体常数，未给出其数值。本项目另行选择现代 SI 的 `R = 8.31446261815324 J/(mol K)`，来自 `N_A k_B`。[BIPM 的 SI 定义常数](https://www.bipm.org/en/measurement-units/si-defining-constants) 给出精确的 `N_A=6.02214076e23 mol^-1` 与 `k_B=1.380649e-23 J/K`；也核对了 [NIST 常数表](https://physics.nist.gov/cuu/Constants/Table/allascii.txt)。这个十进制 SI 值是项目的外部约定，不冒充论文采用的历史常数；C++ 将它舍入到选定的基础浮点类型。

## 2. 接口与数据来源

```cpp
#include <mpmc/thermodynamics/pr76_pure.hpp>

// parameters 是来自合法资料并已通过契约校验的 PrParameterSet。
const auto pure = mpmc::thermodynamics::Pr76Pure<double>::from_parameters(parameters, index);
const auto values = pure.evaluate(300.0); // 温度单位 K。
// values.a: Pa m^6/mol^2；values.b: m^3/mol。
```

这是依赖调用方参数的接口片段，不提供虚构的真实物性。`index` 对应快照的请求顺序，越界抛 `out_of_range`。每个对象复制该组分 `PrPureRecord`、数据集 ID/版本和 `Applicability`，不保留源快照视图；来源、单位和人工数据标记不会被擦除。销毁原快照后仍可使用，增减/重排组分时从新快照重新准备对象。没有全局参数注册表或仅按组分数量命中的缓存。

`Pr76Pure<T>` 的 T 支持 float/double/long double；参数原来按 double 存储，转换到 T 时拒绝超出范围或把非零值变为零的转换，正常舍入不会补充原始数据精度。对象预计算 `sqrt(T_c)`、`a_c`、`b` 和 `kappa`；只读 getter 暴露系数与来源，不允许半更新或多字段赋值。可拷贝/移动构造，移动后的源对象不用于计算；来源引用不能从临时对象取得。

`evaluate` 接受同精度 T 或现有 `Dual<T,N>`，返回 `Pr76PureValues<Number>{a,b}`；混合基础精度和整数温度在接口约束阶段拒绝。其他满足结构约束的数值类型仍需提供保导数运算及 ADL sqrt，未验证其兼容性；没有高阶或参数拟合导数承诺。不同方向的种子可以任意设定，`b` 始终构造成同类型常量，全部温度导数为零。

## 3. 数学与数值契约

设 `q=1+kappa*(1-sqrt(T/Tc))`，实现 `a=a_c*q*q`。数学温度域为有限 `T>0`；Tc/Pc/omega 作为常量模型参数，不参与本轮求导。给定参数声明了温度区间时，在其包含端点的区间内求值，否则抛 `domain_error`；未知范围保持未知。接口不接收压力，不凭空构造压力或宣称已通过完整 T/p、组成或相态适用性检查。

原文叙述其 alpha 拟合使用从正常沸点到临界点的蒸气压数据。这不是对所有组分/温度的普适精度保证。本内核在未知经验范围或超临界温度上执行平方表达式，是明确的**代数延拓**而非物理验证。特别是 q<0 时，平方仍可计算，但不再声称满足非负主平方根 `sqrt(alpha)=q` 的原式形式；不自动裁剪 q、取绝对值、修改 alpha 或伪造有效性。需要禁止外推时，调用方必须提供经过核验的适用区间/领域规则。

普通浮点和 AD 共用表达式；只读取原值和导数作有效性检查，不提取普通数值重建 a，也不手工给 AD 注入解析斜率。通过 `using std::sqrt` 加未限定调用让 ADL 找到 AD sqrt。代码把 `sqrt(T/Tc)` 改写为 `sqrt(T)/sqrt(Tc)`，避免先形成极端温度比；a_c 不先计算 `(R Tc)^2`，a 也不先计算 `q^2`。这些是等价的求值顺序选择，会改变最后几位舍入，不是更换模型。

| 情形 | 行为 |
| --- | --- |
| T<=0、NaN/Inf 温度、非有限 AD 种子、已知温度区间外 | `domain_error`；不修改输入或准备对象 |
| 参数窄化越界/非零变零、预计算 a_c/b 非有限或非正、kappa 非有限 | `range_error`；不返回准备对象 |
| 温度表达式/最终 a 或其导数非有限，或非零 q 的正吸引参数下溢成零 | `range_error`；不返回半成品 |
| q 恰为零且可表示 | a 可为零，平方表达式的一阶导数为零；无额外奇点 |
| 分配失败或底层 AD 数学异常 | 原样传播，不裁剪或以默认数值冒充成功 |

这是经过检查的有限精度核，不是保证任意有限输入均能求得可表示最终结果的全动态范围算法。有些中间量不可表示时，即使重排或更高精度可能得到有限结果，本实现仍会拒绝；不能检测所有舍入误差或导数下溢。没有开启 fast-math、march=native 或强制 FMA，没有改变 AD 现有定义域/异常规则。

配置复制来源会分配；对支持的 T/Dual 的正常 evaluate 路径不分配、不使用可变缓存。常量计算一次复用，成本不含调用方遍历和 Jacobian 存储；没有性能基准，不宣称实时或优于成熟实现。

## 4. 导数判据与独立增量测试

以下是由已核验公式**独立求导**得到的验收判据，不是论文直接印出的导数表：

```text
s = sqrt(T/Tc), q = 1 + kappa*(1-s)
da/dT = -a_c*kappa*q / sqrt(T*Tc)
db/dT = 0
T = Tc: a = a_c, da/dT = -a_c*kappa/Tc, db/dT = 0
```

生产代码不调用解析导数函数。测试另用展开式 `a_c*((1+k)^2 - 2*k*(1+k)*sqrt(T/Tc) + k^2*T/Tc)` 及其解析导数；还用独立 Decimal 70 位计算核对固定参考：制造的 Tc=400 K、Pc=4e6 Pa、omega=0.2、T=100/400/900 K。数字来自纸面公式的独立十进制运算，不来自待测 C++ 输出，也不是实际流体或论文实验数据。

新增 `tests/thermodynamics/pr76`，11 个 CTest：printed_reference、analytic_temperature、chain_rule、fixed_jacobian、runtime_jacobian、snapshot_lifetime、temperature_domain、declared_bounds、algebraic_extension、numerical_range、headers。前十项覆盖三种基础浮点精度，header 项使用独立普通浮点翻译单元。

包含原文系数与临界温度恒等式、负/零/高 omega、任意/零种子与温度复合链式法则、固定 2x1 Jacobian，以及运行期组分快照和分块接口的必要集成。运行期测试在同一进程以 1->5->2->4 维输入消费新快照，5 维时包含尾块；不同组分分别拥有独立温度仅为导数矩阵测试，并非相平衡物理模型。来源寿命、温度域、已知/未知适用性、超范围转换与大/小尺度也有对应检查。

数值判据使用 `128*epsilon(T)` 相对容差、零绝对容差；预期为零时精确比较，以免微小导数被零替代。Release 不依赖 assert。没有真实临界常数表、实验精度、混合规则、EOS 根、闪蒸或性能验证；未做 OOM 注入、并发竞态、macOS/其他架构测试。

```sh
cmake -S tests/thermodynamics/pr76 -B build/thermo-pr76 -DCMAKE_BUILD_TYPE=Debug
cmake --build build/thermo-pr76 --target mpmc_thermodynamics_pr76_tests --config Debug
ctest --test-dir build/thermo-pr76 -C Debug -R '^thermo[.]pr76[.]' --output-on-failure --no-tests=error
```

## 5. CI 与下一步边界

沿用 Thermodynamics contracts 工作流身份，新增 contracts/pr76 选择器及 17 组路由回归。本次只选 pr76；原合同套件及 AD 各套件不构建、不运行。新数值核的 AD/契约集成测试属于本轮必要增量，不是旧套件全量重跑。后续共享 components/pr_parameters 或热力学构建变化选择两套；AD 公共头/构建变化选择新数值集成套件，旧 AD 工作流按原规则自行覆盖受影响 AD 测试。

官方 GCC Debug+ASan/UBSan、Clang Release、MSVC Release 是本增量的目标验证环境；实际结果以相应提交的 Actions 日志为准。原合同作业的编译器/选项不变；路由不替代语义审计，工作流变更影响旧编译条件时仍须补充必要测试。没有修改分支保护或必需检查规则。

下一步建议单独实现运行期组分集合的经典 PR 混合参数 a_mix(T,x)、b_mix(x)，消费本次纯组分核和显式 kij；先明确组成独立变量/归一化契约，新增解析组成/温度导数与组分顺序不变性测试，不进入 EOS 根或闪蒸。SW/CPA 的专属物理关联式仍须单独核验，不能把这个纯 PR76 核自动当作它们的完整实现。
