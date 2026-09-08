# 运行期维数的分块 Jacobian

公共头：`<mpmc/ad/runtime_differentiate.hpp>`。`value_and_jacobian_runtime<K>` 在输入数量 n、输出数量 m 于运行期确定时，同时取得函数值与完整一阶稠密 Jacobian。内部继续使用现有 `Dual<T,K>`；不修改 `dual.hpp`、`math.hpp`、固定维数 `differentiate.hpp` 或其测试，不创建 `DynamicDual`，不依赖其他项目模块或第三方库。

## 1. 设计与依赖边界

K 是一次传播的方向数，不是组分数或变量数量上限。使用运行期长度的输入/输出容器和固定宽度数值类型，避免为每个 n 重新编译，或立即重写全部动态标量运算。参考 Ceres 的 DynamicAutoDiffCostFunction 固定 Stride、多次求值策略 [R1]，独立实现；不复制、引入或执行 Ceres 代码。K 必须为正，具体最优批宽需要真实工作负载测量，本轮不作最优宽度或实时延迟保证。

依赖单向：新头仅包含 AD 自身 `dual.hpp` 和 C++20 标准库。旧头不反向包含新驱动，固定接口的一次回调语义不变。新的独立测试会调用 math.hpp 和固定接口作集成检查，但不重跑旧算术/数学/固定 Jacobian 套件。模型参数、组分身份、组成约束、单位、相态和服务任务版本不进入 AD。

## 2. 调用与结果

```cpp
#include <mpmc/ad/runtime_differentiate.hpp>
#include <span>
#include <vector>

using Number = mpmc::ad::Dual<double, 4>;
mpmc::ad::RuntimeJacobianWorkspace<double, 4> workspace;
const auto function = [](std::span<const Number> x, std::span<Number> y) {
    y[0] = x.front() + x.back();
};

// 同一个已编译接口与工作区可以接收不同长度；不需要增设模板实例。
std::vector<double> input(5, 2.0);
const auto first = mpmc::ad::value_and_jacobian_runtime<4>(
    function, input, 1, workspace, {100, 100, 10000});
// first.values = {4}; first.jacobian = {1, 0, 0, 0, 1}; 回调两次。
input.assign(2, 3.0);
const auto second = mpmc::ad::value_and_jacobian_runtime<4>(
    function, input, 1, workspace, {100, 100, 10000});
// second.values = {6}; second.jacobian = {1, 1}; 回调一次。
```

完整入口：

```text
value_and_jacobian_runtime<K>(function, inputs, output_count, workspace, limits = {})
```

T 由 `RuntimeJacobianWorkspace<T,K>` 推导；输入可从相同 T 的 vector、array 或 span 转换为只读 `std::span<const T>`。调用方必须提供有效的连续范围；不接受伪造指针/长度，也不静默进行不同浮点精度之间的容器转换。T 为无 cv 限定的内建浮点类型，并须支持 quiet NaN，用于检测未写入输出；支持范围仍以实际平台测试为准。

回调可概括为 `void(std::span<const Dual<T,K>>, std::span<Dual<T,K>>)`。n 和 m 必须为正；标量输出表示为 m=1。回调收到全部 n 个变量和恰好 m 个输出槽位，无法通过调整 span 的局部副本改变驱动的输出形状。常量输出必须显式构造成 `Number{constant}`。不同基础精度、方向数、普通浮点输出缓冲区、非 void 返回值、只接受可变输入以及仅能作为右值调用的回调会被约束拒绝。

返回 `RuntimeValueAndJacobian<T>`，独立拥有数据：

| 成员 | 语义 |
| --- | --- |
| `input_count` / `output_count` | 运行期 n / m。 |
| `values` | 长度 m 的 `std::vector<T>`，保存第一批原值。 |
| `jacobian` | 长度 m*n 的 `std::vector<T>`；第 i 行第 j 列在 `i*n+j`，连续行主序。 |

结果在所有批次成功后才返回；不与输入或工作区共享存储，工作区销毁后仍有效。与固定结果一样，这些是可修改的公开数据成员，调用方修改后自行保持尺寸一致；不承诺二进制 ABI 或网络序列化格式。

## 3. 分块算法与回调契约

n=10、K=4 时，三次回调分别播种第 0–3、4–7、8–9 列。每次原值和变量顺序完全相同，尾批两条多余方向为零，不增加虚假变量或零组分。成功调用次数为 `n/K + (n%K != 0)`，无额外探测/普通数值调用。

输入原值先检查并复制到私有工作区；不在后续批次重新读取外部输入视图。初始化一次全部零导数，此后只播种当前块、清除当前块，不每批清空全部 n*K 个方向。每批输出都先置为 NaN 哨兵，回调必须从头完整赋值每个输出；不能在尚未初始化的输出上使用 `+=`。这防止遗漏输出时误用上一批的值或零默认值。

回调不被复制、移动或保留；即使传入临时函数对象，每批也作为**同一个左值**调用。这不同于固定接口的单次完美转发。允许计数等不改变数学结果的状态，但必须满足重复求值契约：整个求导期间表达式、参数、输出数量与含义保持不变；不能因调用次数或播种方向改变原值、复用过期 AD 导数、推进时间步或改变相集合。

每批输出原值须与第一批精确相等，且符号零一致。不采用任意容差接受不同问题的原值。并行非确定性等导致原值变化时，该接口会拒绝；这不是函数数值精度容差。原值一致只是诊断防线，无法证明回调纯净，也无法识别剥离 `.value()` 后重建 Dual、在相同原值下篡改导数或使用过期缓存。输入/输出 span 仅在当前回调期间有效，不得保存后使用；C++20 span 本身并不提供越界检查 [R2]，回调须遵守边界和只读输入要求。

## 4. 校验、异常与失败恢复

此新驱动采用**有限输入、有限输出和有限活动方向导数**的严格契约。它不修改固定驱动允许非有限输出透传的既有行为；异常检查也不把不适用点变成可微点。

| 情形 | 行为 |
| --- | --- |
| n=0 或 m=0；空函数指针 | `std::invalid_argument`，回调不执行。 |
| n*m 不可表示、超过 vector 容量限制或调用方限额 | `std::length_error`，在分配与回调之前拒绝。 |
| NaN/±Inf 输入 | `std::domain_error`，在分配与回调之前拒绝。 |
| 输出未完整赋值、原值非有限或活动方向导数非有限 | `std::domain_error`；不返回半成品，不裁剪或改造成有限值。 |
| 不同批次的原值不一致，含 +0/-0 变化 | `std::runtime_error`；拒绝拼接不一致的矩阵。 |
| 同一工作区重入 | `std::logic_error`。此保护不是互斥锁，不允许并发共享工作区。 |
| 分配失败或回调抛异常 | 保留异常传播；工作区可用于后续新调用，已有返回结果不受影响。 |

正常可表示的零值/零导数，包括下溢后的零，不额外拒绝；非活动尾方向不作为结果列提取。回调的外部副作用不能回滚。失败后工作区内部内容未指定，但下一次调用重新初始化；不存在一个可被调用方误读为成功的“部分结果”缓存。没有注入实际 OOM 来验证 `std::bad_alloc`，此传播依赖标准容器/RAII 语义，不声称完成分配器故障测试。

`RuntimeJacobianLimits` 提供 `max_inputs`、`max_outputs`、`max_jacobian_entries`，均为包含端点的逻辑问题上限；默认只限制可表示范围和标准容器容量，**不是服务端安全配额**。对外部请求，服务层必须按资源预算显式设置上限。限额不包括回调自己的分配、分配器元数据、工作区已保留容量或同时存活的旧结果，不能把它当作进程总内存上限。

## 5. 所有权与成本

工作区拥有两块可复用的 AD 向量，尺寸缩小不会主动释放容量，重新增长到现有容量内不要求扩容；`input_capacity()` / `output_capacity()` 可观察容量。释放长期保留内存可结束工作区生命周期。工作区不可复制或移动，防止回调期间移动它使视图失效；并发任务使用各自工作区，不能把重入标志当作线程同步。

工作区规模 O(K*(n+m))，完整稠密结果 O(m*n)；目标函数的内部暂存另计。驱动的播种初始化/块更新为 O(n*K)，输出初始化/扫描与提取为 O(m*K*ceil(n/K))，另有每批原值计算的重复成本。工作区首次创建或容量增长时可能分配，每个返回结果仍分配 values 和 Jacobian；**不承诺整次调用零分配**，但没有逐标量动态梯度或驱动每批新增分配。没有性能基准，不能宣称分块一定比全方向传播更快。

这里提供运行期变尺寸的一阶稠密求导能力，不提供组分管理、隐式求解敏感度、稀疏组装、反向/二阶 AD、前端任务取消、缓存失效或实时延迟保障。可支持这些后续模块的数值入口，但不代表完整软件已具备它们。

## 6. 独立增量验证

```sh
cmake -S tests/ad/runtime -B build/ad-runtime -DCMAKE_BUILD_TYPE=Debug
cmake --build build/ad-runtime --target mpmc_ad_runtime_tests --config Debug
ctest --test-dir build/ad-runtime -C Debug -R "^ad[.]runtime[.]" --verbose --no-tests=error
```

新工程只消费 `mpmc::ad`，不配置原算术、初等函数、固定 Jacobian 或旧独立消费者。正式证据来自对应提交的 GitHub 官方 runner；本次仅选择 runtime 套件及受影响 CI 路由检查。原根 preset 不变。新公共模板使用必要的 GCC Debug + ASan/UBSan、Clang Release、MSVC Release 验证，不以普通容器运行结果替代官方 CI。

| 新增 CTest 后缀 | 判据 |
| --- | --- |
| `block_seeding` | K=1/3/4/8；n=1/2/3/4/5/8/9/10；完整原值、调用次数、尾批与单位矩阵。 |
| `rectangular_analytic` | 运行期 n/m 的非方阵函数族，独立手推逐元素导数。 |
| `fixed_crosscheck` | 5 输入、3 输出、K=3；解析 Jacobian 与固定接口双重核验。 |
| `constants_unused` | 常量零行与未使用输入零列跨三批保持正确。 |
| `resize_ownership` | 同一工作区 12→2→7→1→12，容量复用、结果独立与工作区销毁后的结果。 |
| `callback_contract` | 非复制回调、临时对象的重复左值调用、reference_wrapper、函数指针与空指针。 |
| `failure_recovery` | 第二批抛异常时保留原异常载荷，不覆盖已有结果；下一次调用正确。 |
| `input_validation` | 各输入位置的 NaN/±Inf 在分配与回调前拒绝，输入不变。 |
| `size_limits` | 空维数、三个限额、包含端点与巨大输出数；不伪造超长输入 span。 |
| `output_validation` | 遗漏第二批输出、非有限值/导数、不使用残留输出及失败恢复。 |
| `primal_consistency` | 随调用变化的原值、符号零变化不能拼成一个 Jacobian。 |
| `workspace_isolation` | 同工作区重入拒绝；不同工作区嵌套不污染种子。非并发竞态测试。 |
| `scaled_derivatives` | 跨批次的微小非零导数，绝对容差为零，不以零导数蒙混通过。 |
| `header_odr` | 新公共头首个/重复包含、双翻译单元与独立消费 AD。 |

前 13 个条目覆盖 float/double/long double，头文件链接项使用 double；接口约束另用编译期检查。解析函数是明示的人工构造数值测试，不是实验数据。短无量纲表达式使用 64 倍对应精度 epsilon 的绝对/相对容差；微小导数只有相对容差。Release 验证不依赖 assert，不复制旧初等函数逐项用例。

CI 的 runtime 头/测试路径只选择 runtime；固定 differentiate.hpp 变更需包括 jacobian 和 runtime；math.hpp 变更需包括 math、jacobian、runtime；dual.hpp 或共享构建变化仍包括所有受影响套件。未知依赖保守处理；CI 路由不代替语义审计，伴随代码的工作流修改影响其他套件的编译环境时必须补充必要验证。

## 7. 参考依据

- [R1] [Ceres DynamicAutoDiffCostFunction 官方源码](https://github.com/ceres-solver/ceres-solver/blob/master/include/ceres/dynamic_autodiff_cost_function.h)，本轮读取前 90 行，blob `05b13bc8da01239bfbb4d8784e30ecfc18aab3bd`；确认固定 Stride 分批思路与 BSD 条款，未全面审计或复制实现。
- [R2] [C++20 工作草案 N4861：span](https://timsong-cpp.github.io/cppwp/n4861/views.span)，用于视图所有权、连续范围、转换及索引前置条件。

项目自身许可证仍待负责人决定；参考成熟软件不构成采用它的源码、依赖或项目许可证。
