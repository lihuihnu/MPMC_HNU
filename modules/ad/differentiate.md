# 固定维数函数值与 Jacobian

公共头：`<mpmc/ad/differentiate.hpp>`。接口 `mpmc::ad::value_and_jacobian(function, inputs)` 独立于原算术与初等函数实现；只包含本 AD 模块的 `dual.hpp` 与 C++20 标准库，不依赖其他项目模块、第三方 AD 或矩阵库。调用方需要初等函数时自行包含 `math.hpp`。

## 1. 本增量的边界与设计

现有 `Dual<T, N>` 已能传播 N 个一阶方向，本接口只统一单位基播种、一次函数求值和导数提取。与让每个领域模块手写这些步骤相比，这样更容易固定索引和所有权契约；与新增动态梯度、求导图或矩阵库相比，不改变已验证的数值类型和依赖边界。

参考 autodiff 官方教程中同时返回向量函数值与 Jacobian 的接口分工，以及 Ceres 的多方向 Jet 思路；采用独立实现，没有复制、引入或运行这些库。标准调用规则参考 C++20 `std::invoke`，使用完美转发而非 `std::function` 类型擦除。原始入口见文末。

本次没有修改 `dual.hpp`、`math.hpp`、原算术/初等函数测试、原根与模块 CMake 或 Presets。新增独立测试工程不配置旧测试工程；仅为新套件扩展 CI 路由。未来修改基础 AD 应扩展到本接口，修改 `math.hpp` 应包含本套件中的数学表达式集成测试。

## 2. 类型、调用与矩阵布局

输入为 `const std::array<T, N>&`，其中 T 是无 cv 限定的内建浮点类型，N > 0。回调接收 `const std::array<Dual<T, N>, N>&` 或其按值副本，并按值返回 `std::array<Dual<T, N>, M>`，M > 0。两种维数均自动推导，不需要显式模板实参。

返回类型为 `ValueAndJacobian<T, N, M>`：

| 成员 | 类型与约定 |
| --- | --- |
| `values` | `std::array<T, M>`，第 i 项为 F_i(x)。 |
| `jacobian` | `std::array<std::array<T, N>, M>`；`jacobian[i][j] = ∂F_i/∂x_j`。 |
| `input_count` / `output_count` | 编译期常量 N / M；行对应输出，列对应输入原始顺序。 |

返回对象拥有自己的数值与矩阵；矩阵按嵌套行数组访问，不承诺可对首行指针进行跨行扁平遍历或提供稳定二进制 ABI。M=1 时仍返回长度为 1 的向量和一行 Jacobian，不另加标量重载。

回调以单位基播种的独立输入为参数，**一次调用同时传播所有 N 个方向**。输出长度通过编译期类型推导获得，不预先执行回调；不会为获取普通函数值再调用一次。支持函数指针、捕获 lambda、`std::ref`、不可复制及有左右值调用限定的函数对象，不复制或保留函数对象。空函数指针等不可调用的运行时状态由调用方保证不存在。

回调必须保留 AD 运算链。普通浮点输出、错误的 Dual 精度/方向数、零长度输入/输出、单个标量输出、输出引用或视图、仅能接收可变输入引用的回调会在接口约束处被拒绝。常数输出应写成 `Number{constant}`；故意提取 `.value()` 后重建 Dual、自行篡改种子或使用过期缓存，无法只靠返回类型自动识别，不得将其当作有效求导。

## 3. 可微域与异常

所有输入原值必须有限；NaN 与正负 Inf 抛出 `std::domain_error`，发生在任何回调调用之前。回调内的定义域、除零或其他异常原样传播，保留动态类型和有效载荷；不吞掉异常返回默认结果。仅在完整求值与提取完成后返回结果对象。

接口不通过参数修改调用方输入；播种数组只在本次调用内有效，回调不得保存其引用、指针或视图供后续使用。回调通过捕获修改外部对象的副作用无法回滚，线程共享的回调状态也由调用方同步；本接口没有共享可变工作区。

有限输入并不保证有限输出或导数。输出中的 NaN/Inf 按实际回调结果提取，不裁剪、不额外抛溢出异常、不自动修复；调用方必须检查适用性和结果有效性。这里不改变原算术与初等函数各自的异常策略。

求导针对选定表达式在当前点的普通一阶导数；可微域为该表达式各操作实际允许的可微域。接口不把分段分支变光滑，也不证明隐式迭代收敛解的导数正确。所有输入坐标在播种时视为独立变量；组成归一化、约束消元、变量单位及尺度转换属于调用方的数学建模，不在本接口偷偷加入。额外捕获的普通参数视为常数。

本层的播种存储/初始化为 O(N²)，提取和结果存储为 O(MN)，单次回调内部 AD 运算仍有 N 方向成本。无本层正常路径堆分配、虚调用或录制 tape；回调和异常路径不作零分配承诺。适合小型/局部稠密问题，不能将全局网格未知量数直接当作 N。未运行性能基准，不声称比手写或成熟 AD 库更快；不承诺 C++20 编译期求值。

## 4. 使用示例

```cpp
#include <mpmc/ad/differentiate.hpp>
#include <mpmc/ad/math.hpp>
#include <array>
#include <cmath>

const auto function = [](const auto& p) {
    using std::exp;
    using std::log;
    return std::array{p[0] * p[0] + p[1], p[0] * p[1], exp(p[0]) + log(p[1])};
};
const auto result = mpmc::ad::value_and_jacobian(function, std::array{2.0, 3.0});
// values = {7, 6, exp(2) + log(3)}
// jacobian rows = {{4, 1}, {3, 2}, {exp(2), 1/3}}
```

泛型数学调用采用 `using std::...` 后不限定命名空间的调用方式，使普通浮点与 AD 重载分别选用正确实现。这个 3 输出、2 输入示例已纳入新增解析 Jacobian 测试；是否通过以对应提交的官方 runner 日志为准。

## 5. 独立增量验证

```sh
cmake -S tests/ad/jacobian -B build/ad-jacobian -DCMAKE_BUILD_TYPE=Debug
cmake --build build/ad-jacobian --target mpmc_ad_jacobian_tests --config Debug
ctest --test-dir build/ad-jacobian -C Debug -R "^ad[.]jacobian[.]" --verbose --no-tests=error
```

这是独立工程，不会运行旧 `ad.dual`、`ad.math.*` 或原独立消费者。根 `ad-debug` preset 不变。正式执行证据来自 GitHub 官方 runner，本轮计划仅运行 jacobian 套件的 GCC Debug + ASan/UBSan、Clang Release、MSVC Release 作业；CI 路由自身的新增/受影响选择检查也是必要增量，不是旧数值测试重跑。

| 独立 CTest 条目后缀 | 独立判据 |
| --- | --- |
| `identity_seeds` | N=4 的单位基、输入顺序、调用次数和单位 Jacobian。 |
| `rectangular_analytic` | 3×2 非方阵，9 个输入点；手工解析行 `[2x,1]`、`[y,x]`、`[exp(x),1/y]`。 |
| `constants_unused_inputs` | 常量输出为零行；未使用输入为零列；N=3、M=2。 |
| `single_input_output` | N=M=1，x³ 在 x=-2 的值 -8、导数 12。 |
| `callable_forwarding` | 不可复制回调、左右值限定、`std::ref`、函数指针；恰当调用次数。 |
| `repeated_calls_ownership` | 无残留种子、无共享工作区、输入不变、返回数据独立持有。 |
| `exception_propagation` | 自定义异常的类型/载荷、数学域异常和输入不变。 |
| `nonfinite_inputs` | 每个输入位置分别使用 NaN/±Inf；回调调用数必须为 0。 |
| `output_passthrough` | 刻意构造非有限结果，确认提取不伪造有限输出；不是物理验证。 |
| `header_odr` | 公共头首个/重复包含，两个翻译单元链接，独立消费 AD 目标。 |

前九项覆盖 float/double/long double，公共头链接项使用 double。类型/维数/输出所有权约束另由 `static_assert` 编译检查。精确多项式和单位基使用精确比较；短无量纲数学表达式采用 `64*epsilon(T)*max(1,abs(expected))`，参考导数不是由待测 AD 结果回填。运行期测试不会被 NDEBUG 删除。不重做旧初等函数逐项测试，也不虚构未实现下游测试。

## 6. 原始参考

- [autodiff 官方教程：向量函数的 Jacobian 与同时取得函数值](https://autodiff.github.io/tutorials/#jacobian-matrix-of-a-vector-function)。仅参考接口分工，不使用其 Eigen 适配或源码。
- [Ceres 自动微分：多方向 Jet 与单位基](https://ceres-solver.org/automatic_derivatives.html)。现有 Dual 的参考背景，不新增 Ceres 依赖。
- [C++20 工作草案 N4861：std::invoke](https://timsong-cpp.github.io/cppwp/n4861/func.invoke)。用于核对回调转发与返回类型推导规则。
