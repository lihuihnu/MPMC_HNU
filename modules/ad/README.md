# 独立 AD 模块：首个算术增量

`mpmc::ad::Dual<T, N>` 是提供给其他数值模块的一阶前向自动微分数值类型。实现只有标准库头文件依赖，**不依赖 `core`、其他项目模块、Eigen、Ceres、autodiff 或 CppAD**。CMake 目标为 `mpmc::ad`；也可只复制 `include/` 并启用 C++20 使用。

本增量只建立固定维数的值/导数语义、种子和算术，不宣称已经完成全部 AD 能力。测试结果以对应提交的 GitHub Actions 日志为准；未做性能基准，也没有热力学或实验数据验证。

## 1. 设计取舍与参考

设计参考 Ceres 的 Jet 数值模型与算术接口、autodiff 的标量替换用法，以及 CppAD 的前向方向导数语义：[R1][R2][R3]。采用独立实现，不复制或打包这些项目的源码。已核对 Ceres 源文件中的 BSD-3-Clause 条款、autodiff 的 MIT 条款及 CppAD 的许可说明；项目自身许可证仍待负责人决定，引用不构成重新授权。

选择固定长度 `std::array` 而非动态梯度、表达式模板或计算图：本阶段容易审计所有权、自别名和链式法则，正常算术无需堆分配、虚调用、全局状态或录制 tape。一般 Dual 算术的工作量和存储随 `N` 线性增长；标量加减的原地操作只修改值。异常对象可能分配内存，不承诺错误路径零分配。没有基准证据前，不宣称比成熟库更快。

后续可以基于测量重构存储或求导策略，但不得改变既定导数语义而不说明。大规模全局问题应按需求评估方向播种、分块/稀疏传播或反向模式，不能把全局未知量数量直接等同于每个标量的 `N`。

## 2. 类型与接口契约

| 接口或约束 | 语义 |
| --- | --- |
| `Dual<T, N>` | `T` 为无 cv 限定的内建浮点类型；`N > 0`，在编译期确定。常规使用 `float`、`double`、`long double`。默认 `N = 1`。 |
| `Dual{}` / `Dual{value}` | 常数，导数全部为零；标量构造为 `explicit`。 |
| `Dual::variable(value, index)` | 第 `index` 个导数分量置 1，其余为 0。越界抛出 `std::out_of_range`，Release 中也检查。 |
| `Dual{value, gradient}` | 显式提供任意方向种子。各输入使用同一方向基；类型不会自动识别物理变量。 |
| `value()` / `derivative(i)` | 按值读取数值/单个导数；后者检查索引。没有到普通浮点数的隐式转换。 |
| `derivatives()` | 对左值返回只读数组引用；对临时对象返回拥有数据的数组，避免悬垂引用。左值引用的生命周期仍由持有者负责。 |
| 算术 | 支持一元 `+/-`、Dual 间及两侧标量的 `+ - * /`、对应复合赋值。自别名 `x *= x`、`x /= x` 有专门回归。 |
| 类型混合 | 不隐式混合不同 `T` 或 `N` 的 Dual。普通标量遵循到 `T` 的 C++ 转换规则，不自动提升精度；精度敏感的常数显式写成 `T{...}`。 |
| 替换与并发 | 用 `x = Dual{constant}` 清除旧导数；对象持有独立数据。不同对象可独立计算，同一可变对象的并发写入由调用者同步。 |

`N` 表示传播的方向数，不一定是物理输入数量。单位基种子给出梯度/Jacobian 各列；单方向种子可在一次计算中得到 `J*v`。所有读数、种子和合法算术可用于 `constexpr`；非法常量求值不能生成有效编译期结果。

### 算术与异常

对 `q = u/v`，采用 `q' = (u' - q*v')/v`，而不形成 `v*v` 或 `1/v`；乘法使用 `u'*v + u*v'`。复合运算在覆盖数据前保留需要的原值，并只读写同一导数分量，因此不需要整份梯度的额外工作区。

分母为 `+0` 或 `-0` 时抛出 `std::domain_error`，失败的 `/=` 保持原对象不变；不会用任意 epsilon 判断非零分母，也不会静默截断小数或夹紧导数。该保护是本模块的明确契约，不是要求底层浮点除法也抛异常。

对 NaN、Inf、溢出和下溢，保留所用浮点运算的行为，不生成虚假的有限结果；不承诺非有限输入有有效导数。乘积、中间量和导数仍可能溢出，固定表达式也可能存在相消；上述除法安排不是对任意动态范围的稳定性证明。禁止用 `fast-math` 或默认 flush-to-zero 改写已测试的数值语义。

本阶段**不提供** `exp/log/sqrt/pow`、比较/分支运算、动态维数、高阶/嵌套 AD、反向模式、稀疏 Jacobian 容器或 Eigen 适配。不能把 `std::log(x.value())` 等显式剥离数值的计算冒充 AD；所需初等函数将在后续小增量中加入并验证定义域。

## 3. 用法

```cpp
#include <mpmc/ad/dual.hpp>

using Number = mpmc::ad::Dual<double, 2>;
const auto x = Number::variable(2.0, 0);
const auto y = Number::variable(3.0, 1);
const auto f = x * x + x * y;
// f.value() == 10.0
// f.derivative(0) == 7.0; f.derivative(1) == 2.0
```

通用数值函数应以类型参数代替硬编码 `double`，并保留 AD 类型到表达式末端。以上例子也作为独立消费者测试编译并运行，不依赖平台根构建或测试辅助代码。

在使用者的 CMake 项目中：

```cmake
add_subdirectory(path/to/modules/ad ad-build)
target_link_libraries(your_target PRIVATE mpmc::ad)
```

这里只声明源树内使用方式，不承诺已经提供安装包、稳定二进制 ABI 或 `find_package` 导出配置。`modules/ad` 自身也可以单独配置；INTERFACE 库没有待编译源文件，实际编译由消费者触发。

## 4. 测试入口与证据边界

需要 C++20 编译器和 CMake 3.21 或更新版本，无测试框架下载。根目录入口为：

```sh
cmake --preset ad-debug
cmake --build --preset ad-debug
ctest --preset ad-debug
```

这些命令供复现使用；本项目的正式执行证据必须来自 GitHub 官方托管 runner。`tests/ad/dual_test.cpp` 包含 12 个命名用例组，使用正常运行期检查而不是可能被 `NDEBUG` 删除的 `assert`。CTest 将其注册为 `ad.dual`，标签 `ad;unit`。

| 验证范围 | 判据 |
| --- | --- |
| 构造、播种、存储与类型约束 | 常量导数为零；独立拷贝；临时对象生命周期；1/2/4 维；错误维数/类型在编译期拒绝。 |
| 算术与复合赋值 | 固定解析值；标量在运算符两侧；负值、零值；自别名与复合表达式。 |
| 导数正确性 | 2 输入、2 输出的多项式/有理函数；手工简化的解析 Jacobian；每种基础浮点类型的 16 个确定性输入点；独立 `J*v` 判据。 |
| 边界与数值尺度 | 索引越界；正负零除法和状态不变；`1e-200`、`1e200`、适用平台的次正规分母；不隐藏 NaN/Inf。 |
| 可组合性 | 公共头首个包含、重复包含、两个翻译单元链接，以及只引入 AD 模块的外部消费者。 |

无量纲普通算例使用 64 倍对应 `T` 的机器 epsilon 作为绝对与相对误差界，覆盖短表达式的舍入误差；极小的非零期望导数将绝对容差设为 0，防止错误的零导数被放过。参考导数来自独立手工推导，不从待测 AD 结果回填。常规浮点表达式只用于值传播交叉检查，不被称为独立物理验证。

首轮 CI 覆盖 Linux/GCC Debug + AddressSanitizer/UndefinedBehaviorSanitizer、Linux/Clang Release、Windows/MSVC Release；每个作业运行上述单元测试和独立消费者。它们用于验证新公共模板与初始跨编译器构建，不是通用全量平台矩阵；macOS、其他架构、其他编译器版本和 HPC 性能尚未验证。没有已实现的热力学/闪蒸下游，因此本次没有虚构下游测试。

CI 在相关路径的 PR 上运行，合入后不重复执行同一树；手动 `workflow_dispatch` 可复测指定 ref。代码变更应通过 PR 完成官方验证再集成，不能依靠直接推送 `main` 触发本工作流。本次没有修改分支保护；**将来设置必需检查前，必须把当前事件级路径筛选改成始终报告状态的门禁，避免纯文档 PR 永久等待。**

## 5. 原始参考入口

- [R1] [Ceres 自动微分说明](https://ceres-solver.org/automatic_derivatives.html)及 [`jet.h`](https://github.com/ceres-solver/ceres-solver/blob/master/include/ceres/jet.h)。审阅源文件 blob：`3672011fd6d45d1cc0c119a9ce02e56be3e5a7cd`。
- [R2] [autodiff 官方文档](https://autodiff.github.io/)及 [`forward/dual/dual.hpp`](https://github.com/autodiff/autodiff/blob/main/autodiff/forward/dual/dual.hpp)。审阅文件开头的许可与接口组织，blob：`31b65032fcecf642348b338eda83be956c2fc119`；未声称全面审计该库。
- [R3] [CppAD 一阶前向模式](https://cppad.readthedocs.io/latest/forward_one.html)及 [许可说明](https://github.com/coin-or/CppAD/blob/main/COPYING)。仅借鉴 `J*v` 语义，本轮没有运行 CppAD 对比测试。
