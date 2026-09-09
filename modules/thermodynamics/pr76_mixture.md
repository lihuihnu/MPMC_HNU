# PR76 运行期经典混合参数

本增量只计算 `a_mix(T,x)` 和 `b_mix(x)`，不计算压力、EOS 根、逸度、相稳定性或闪蒸。公共头为 `<mpmc/thermodynamics/pr76_mixture.hpp>`，继续使用 `mpmc::thermodynamics` 目标。普通浮点路径只需标准库和已有热力学头；使用 AD 时由调用方包含 `<mpmc/ad/math.hpp>`。未修改原 AD、参数契约、纯组分实现或其测试。

## 1. 原文依据与固定假设

依据用户提供的 *A New Two-Constant Equation of State*，PDF 第 2 页/期刊第 60 页式 (20)–(22)，已结合清晰页面图像核对求和下标及两个独立平方根：

```text
a_mix = sum_i sum_j x_i x_j a_ij
b_mix = sum_i x_i b_i
a_ij  = (1 - k_ij) sqrt(a_i(T)) sqrt(a_j(T))
```

论文使用 `delta_ij`（δᵢⱼ），本项目沿用 `kij`。这里的 x 是摩尔分数，不是质量分数。PDF SHA256 为 `0ea7cb1a8c6206f9e929d7d7730399bb6eecf310259cd7146d9eaadb4b4a6036`；来源 DOI 为 `10.1021/i160057a011`。仅实现上述混合关系，未复现论文相平衡表格，也未把 PDF/截图/全文加入仓库。

保留 `PR76/classical-vdw/constant-kij/no-translation` 参数配置及 [纯组分数值约定](pr76_pure.md)。所有组分使用同一个温度，kij 是已经显式录入的对称常数，k_ii=0 为结构零；不默认填充缺失 pair，不提供温度/组成相关 kij、体积平移、SW 或 CPA。源参数构造校验仍由 `PrParameterSet` 负责。

原文没有规定 C++ 输入接口、归一化容差或独立组成坐标。下述内容是本项目的明确接口约定，而非冒充论文原句。坐标转换同时参考 [NIST teqp 的 xN (in)dependent 说明](https://teqp.readthedocs.io/en/latest/derivs/compderivs.html#xn-in-dependent)，仅参考导数含义，不引入或运行该库。

## 2. 两个入口，两种导数含义

### 完整组成：evaluate_full

```cpp
mix.evaluate_full(temperature, full_mole_fractions, workspace);
```

输入长度必须恰好为 n，顺序与 `mix.parameters().components()` 完全相同。原值必须有限且逐项处于 [0,1]，AD 种子必须有限。使用补偿求和检查

```text
abs(sum(x_i) - 1) <= 64 * numeric_limits<T>::epsilon()
```

这是无量纲的绝对接受容差：`float` 约 7.63e-6，`double` 约 1.42e-14；`long double` 依平台精度。它只容纳有限舍入误差，不做科学参数校准。**不把 x 除以总和，不改原值或导数种子，不裁剪负分数，即使负数极小。** 接受接近 1 的和后，按原始向量计算，可能保留很小的非归一化残差；调用方需要精确闭合时使用约化入口或显式坐标转换。

各个 x 的种子独立传播，不要求种子和为零。它返回混合多项式的形式偏导，在所接受的组成点求值；不是自动施加 sum(x)=1 后的偏导，也不是输入合法性判断的导数。物理 simplex 的切向方向可由调用方显式播种，不能把一个任意完整梯度直接称为约束梯度。

### 约化组成：evaluate_reduced

```cpp
mix.evaluate_reduced(temperature, independent_mole_fractions, workspace);
```

输入长度必须为 n-1，对应当前快照的前 n-1 个组分。最后一个组分定义为

```text
x_last = 1 - sum(independent)
```

求和和相减都保留 `Number` 类型，由 AD 自动传播负的种子和。各输入及计算出的余量均须有限且处于 [0,1]，否则拒绝；不把略小于零的余量裁剪为零。此入口不做 x/sum(x) 归一化。对于单组分，组成输入是空 span，x_0 为常量 1。

通常以 `[T, x_0, ..., x_(n-2)]` 求 Jacobian，得到两个输出、n 个输入。更换最后的因变量组分会改变坐标图：约化梯度不能仅靠旧梯度列换序获得。

例如对 F 为 a_mix 或 b_mix，有

```text
(dF/dx_l)_reduced = (dF/dx_l)_full - (dF/dx_last)_full
```

这由代入余量并应用链式法则得到。零组成允许，但边界导数是多项式延拓或可行切向导数，不保证任意双向扰动都留在非负 simplex 内。实现不跳过零组成项，否则会丢失该组分的加入导数。

## 3. 准备对象、工作区和返回值

`Pr76Mixture<T>::from_parameters(parameters, limits={})` 从一个合法快照准备模型；T 可为 float、double、long double。复制完整有序快照以保留组分与 pair 来源，逐组分复用 `Pr76Pure<T>`，并预计算同顺序的 `1-kij`。不保留调用者的悬垂视图，不建立可变全局缓存。参数固定，不对 Tc/Pc/omega/kij 求导。

增减、替换或重排组分须先按 ID 创建新的参数快照，再准备新模型。只更改组成数组长度或只重排组成、未同步模型顺序是错误用法；数值向量本身不携带 ID，接口不能识别这种同长度错配。源记录与旧模型不被修改。对象可拷贝/移动构造但不可赋值；移动后源对象不再用于计算。参数引用仅从左值取得，寿命不超过模型。

`Pr76MixtureWorkspace<Number>` 仅保存 n 个根的临时向量，Number 必须与温度/组成类型相同。它不可拷贝/移动，不暴露可变存储；容量可连续复用。**每次求值都重新计算所有根与其导数**，不同 Jacobian 批次不能复用上一批的 AD 值。一个工作区仅供一个同步求值使用；并发任务必须使用不同工作区，并保持输入不变。

返回 `Pr76MixtureValues<Number>{a,b}` 按值拥有两个结果：a 的单位 Pa·m⁶/mol²，b 的单位 m³/mol；不是无量纲 A/B。任一错误抛异常，不返回部分结果。工作区在失败后可能含部分临时值，但下一次求值会在读取前覆盖它们，已经返回的结果不会被覆盖。

实现对角项直接用 a_i，交叉项遍历上三角并累加两次；使用 sqrt(a_i)*sqrt(a_j)，不先形成可能溢出的 a_i*a_j。准备对象储存 O(n²) 参数及元数据；热路径工作区 O(n) 个 Number，每次求值 O(n²) 标量运算（K 方向 AD 约 O(n²K)）。工作区增长/模型复制可能分配；容量足够后的正常内核求值不分配，但外层运行期 Jacobian 返回矩阵仍可能分配。未进行性能比较，不承诺实时延迟或全动态范围稳定。

## 4. 可微域、数值风险与异常

温度与温度种子沿用纯组分检查：T>0、有限，已声明温区受检查；未知域仍为未知。本核没有压力输入，不能判断压力、相态或物理准确性。

多组分 n>1 时，全部选中组分要求 **a_i(T)>0**，包括 x_i=0 或 kij=1 的特殊情形。sqrt(a_i(T)) 在 a_i=0 处可能出现尖点，当前核统一拒绝这种边界，不伪造有限导数，也不采用 epsilon 正则化。单组分不需要交叉平方根，保留原纯组分 a_i=0 的可用行为。对 q_i=1+kappa_i(1-sqrt(T/Tc_i))<0，使用现有平方公式的主平方根 sqrt(a_i)=sqrt(a_ci)*abs(q_i)，不能误用带符号的 q_i*q_j。外推仍只具有代数含义。

有限 kij 不保证 a_mix 正值；本核保留有限的负值/零值，不裁剪，也不把它解释为相稳定或可用 EOS 状态。未来物性/相稳定性层必须按其科学需求检查。

| 异常 | 条件 |
| --- | --- |
| `invalid_argument` | 完整/约化组成尺寸与模型不匹配，或误用移出后的空模型。 |
| `domain_error` | 组成、种子或余量非法；完整组成和不合格；温度违反纯组分约定；多组分有 a_i<=0。 |
| `range_error` | kij 窄化越界/非零变零，1-kij 不可表示；纯组分数值错误；平方根、系数或导数非有限，或 b_mix 不再正值。 |
| `length_error` | 模型准备超出配置/容器规模。默认逻辑上限不是服务内存预算。 |
| 其他异常 | 原参数、AD 数学和内存分配异常保留传播。 |

未检测所有中间下溢或灾难性抵消；不得将“返回有限值”解释为任意尺度的精度保证。Kahan 校验和不改变生产数学表达式；导数传播不经过普通浮点副本。

## 5. 独立导数判据与增量测试

下式由论文混合关系直接推导，未从被测程序输出拟合：令 r_i=sqrt(a_i)，g_i=da_i/dT，则在各 a_i>0、kij 常数时

```text
da_ij/dT = (1-kij) * (g_i*r_j/(2*r_i) + r_i*g_j/(2*r_j))
da_mix/dT = sum_i sum_j x_i*x_j*da_ij/dT
(db_mix/dT)_x = 0
(da_mix/dx_l)_full = 2*sum_j x_j*a_lj
(db_mix/dx_l)_full = b_l
```

约化组成导数按第 2 节取差。测试参考采用原始参数、原文常数字面量、独立纯组分解析斜率与完整双重求和；生产使用上三角。另有二元解析比例：同 Tc/omega，Pc_1=Pc_0/4、x=(1/4,3/4)、k01=1/8 时，a_mix=(95/32)*a_0，b_mix=(13/4)*b_0。

新增独立工程 `tests/thermodynamics/pr76_mixture`，15 个 CTest：二元解析比例、完整梯度、约化梯度、方向链式法则、零组成、4 组分全部 24 排列、运行期 1→5→2→4→1、归一化、温度与种子、pair 规则、主平方根、binary32 零吸引边界、所有权与失败恢复、大小尺度、公共头。其中 13 项覆盖 float/double/long double；零吸引边界采用明确的 float 舍入输入，公共头项为普通 double。24 排列是一个测试条目的子情形，不额外计数。

标准数值比较为 512*epsilon(T) 相对容差、零绝对容差；涉及相消的试例选择非零且可核验判据，精确结构零仍精确比较。测试只使用明确标记 synthetic_test 的人工参数；不是实际流体、实验数据或论文表格验证。源码内给出参考推导与错误诊断，Release 不依赖可禁用 assert。

```sh
cmake -S tests/thermodynamics/pr76_mixture -B build/thermo-mixture -DCMAKE_BUILD_TYPE=Debug
cmake --build build/thermo-mixture --target mpmc_thermodynamics_pr76_mixture_tests --config Debug
ctest --test-dir build/thermo-mixture -C Debug -R '^thermo[.]pr76_mixture[.]' --output-on-failure --no-tests=error
```

官方工作流只选择新增 `pr76_mixture` 套件，不配置/构建/运行旧 AD、参数契约或纯组分套件。以后修改纯组分核时纳入 pure+mixture，修改公共契约时纳入所有相关热力学套件；AD 变化纳入现有下游集成。新增路由回归属于本轮必要测试，路径规则不能替代语义审计。验证状态以实际提交的官方 runner 日志为准。

下一步可单独实现 PR 的显式压力 p(T,v,x)，先明确 v>b_mix 的定义域与单位，复用本混合核并配套解析/AD 偏导测试；不自动启动根求解或闪蒸。
