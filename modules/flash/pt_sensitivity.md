# PR76 汽液 PT 收敛解隐式灵敏度

## 1. 本增量的边界

基线为 `main@0981a8e4a743fd0875c7faf406a131c223cb05a6`。本增量只为**已经被现有 PR76 PT flash 接受的内部汽液两相状态**计算局部一阶收敛解灵敏度，不改变相稳定性、相分裂、Rachford–Rice、SSI、Gibbs 回溯、最终共同切平面复核、状态枚举或任何既有容差。

目标导数列采用约化总体组成坐标

```text
q = (p_Pa, T_K, z_0, ..., z_{N-2})
z_{N-1} = 1 - sum(z_0, ..., z_{N-2})
```

输出为 `betaV`、每个组分 `logK`、液相组成 `x` 和汽相组成 `y` 对上述坐标的 Jacobian。最后一个总体组分是当前有序参数快照中的最后组分；改变组分顺序会改变约化坐标图，但按组分 ID 映射后的物理切向导数应一致。

本实现**不提供**单相导数、跨泡点/露点导数、相消失导数、临界/根切换导数、零进料支持集变化导数、三相导数、SW/CPA 导数或全局光滑性证明。成功状态只表示当前相集合、当前 PR76 简单根分支和当前约化组成图上的局部导数可用。

## 2. 写入前审计：现有未知量和残差

现有 `iterate_pt_split` 的外层数值状态是 N 个 `logK_i`。给定总体组成 `z` 和 `logK`，`solve_rachford_rice` 求内部 `0<beta<1`，再构造满足物料平衡的 `x/y`；相性质回调计算两个候选分支的 `ln(phi)`。外层实际平衡残差为

```text
F_i = log(x_i) - log(y_i) + ln(phi_i^L) - ln(phi_i^V)
    = log(f_i^L/f_i^V), i=0..N-1.
```

因此收敛两相状态可以写成

```text
F(logK; q) = 0,
q = (p,T,z_reduced).
```

不需要把 SSI、Rachford–Rice 二分或 line search 本身作为可微程序。对于固定的已接受相集合，只需对这个收敛方程使用隐函数定理。

PR76 `Pr76Phase::evaluate_reduced` 已经支持 AD 数值类型，且 EOS Z 根导数本身就是对三次方程的局部隐式求导；它明确不对根二分过程求导。因此本增量复用已有热力学导数核，不复制 PR76 公式，也不改变 thermodynamics 接口。

## 3. RR 代数流形的局部隐式微分

定义

```text
K_i = y_i/x_i
F_RR(beta) = sum_i z_i (K_i-1)/(1-beta+beta K_i) = 0
c_i = x_i y_i / z_i
H   = sum_i (y_i-x_i)^2 / z_i
```

在严格正组成和内部两相点上，直接对 RR 方程求微分得到

```text
d beta = [sum_i c_i d(logK_i)
          + sum_i ((y_i-x_i)/z_i) dz_i] / H.
```

随后

```text
d log(x_i) = dz_i/z_i
             - ((y_i-x_i)/z_i) d beta
             - beta (y_i/z_i) d(logK_i)

d y_i = y_i [d log(x_i) + d logK_i].
```

生产实现只取已经接受的 `beta/x/y` 作为线性化基点；**不会调用或微分 RR 二分迭代**。最后一个 `x/y` 分量由 `1-sum(first N-1)` 在 AD 数值类型中构造，使相组成导数严格位于单纯形切空间。

## 4. 平衡方程的隐式求导

令

```text
A = partial F / partial logK
B = partial F / partial q.
```

通过现有分块 `value_and_jacobian_runtime<4>`，将 `logK`、`p/T` 和约化 `z` 的种子沿上述 RR 局部微分和 PR76 相性质导数传播，得到 A、B 以及 `beta/x/y` 对局部变量的偏导。随后只求解

```text
A * d(logK)/dq = -B.
```

再用链式法则得到 `d beta/dq`、`dx/dq`、`dy/dq`。这是**收敛解的局部隐式导数**，不是有限迭代映射的导数。

当前使用带部分主元的稠密 LU。计算 `||A||_inf` 和通过同一 LU 求得的 `||A^-1||_inf`，报告

```text
rcond_inf = 1 / (||A||_inf ||A^-1||_inf).
```

内部主元和 `rcond` 只承担双精度数值可逆性保护，不是严格条件数界。调用方可以通过 `minimum_reciprocal_condition` 提出更严格的局部条件要求。线性求解还独立检查归一化后向误差 `A X + B`；成功后重新检查每列相组成导数和为零以及微分后的逐组分物料守恒。

## 5. 公共接口与失效语义

包含：

```cpp
#include <mpmc/flash/pr76_sensitivity.hpp>
```

使用 `mpmc::flash_sensitivity` CMake 目标；基础 `mpmc::flash` 目标仍只承载原 PT flash 能力。灵敏度适配层额外依赖独立 `mpmc::ad`，依赖方向仍为 `flash sensitivity -> flash/thermodynamics + ad`，AD 不反向依赖任何项目模块。

```cpp
const auto split = mpmc::flash::solve_pr76_pt_vle(p_pa, t_k, z, evaluator);
const auto derivative =
    mpmc::flash::differentiate_pr76_pt_vle(split, evaluator);

if (derivative.sensitivity.status ==
    mpmc::flash::PtSensitivityStatus::success) {
    const double dbeta_dp = derivative.sensitivity.d_vapor_fraction(0);
    const double dx0_dT = derivative.sensitivity.d_liquid(0, 1);
}
```

灵敏度函数从 `split.solution.initial_stability` 读取 flash 实际使用的 `p/T` 和舍入范围内归一化后的 feed，避免调用方再次传入一个不一致的基点。结果也保存这个基点。公共状态：

| 状态 | 含义 |
| --- | --- |
| `success` | 已接受内部两相状态在当前根/相集合/坐标图上得到通过数值检查的局部一阶导数。 |
| `solution_not_accepted` | 输入不是现有完整 flash 接受的两相状态，或保存数据不能重现平衡残差。 |
| `phase_boundary` | 相分率进入灵敏度边界保护、相区分退化，或 RR 局部流形退化；不外推导数。 |
| `unsupported_feed_support` | v1 不对含严格零总体组分的支持集变化求导。 |
| `property_failure` | 当前 PR76 根/物性分支不能提供可靠局部导数。 |
| `ill_conditioned_equilibrium` | `dF/dlogK` 不能可靠分解或不满足要求的 `rcond`。 |
| `arithmetic_failure` | AD、线性求解或微分守恒/归一化检查出现不可表示结果。 |

默认 `minimum_derivative_phase_fraction = sqrt(epsilon(double)) = 2^-26`。它只是避免在相出现/消失附近冒充固定相集合导数的数值政策，**不改变**现有 flash 的 `minimum_phase_fraction=1e-10`，也不改变相数判定。

## 6. 独立导数参考

新增 Decimal(80) 参考不使用生产 `logK/RR/SSI/TPD` 坐标。参考未知量是

```text
s = (x_0,...,x_{N-2}, y_0,...,y_{N-2}, beta)
```

直接联立：

```text
N 个:     log(x_i)+ln(phi_i^L)-log(y_i)-ln(phi_i^V) = 0
N-1 个:   (1-beta)x_i + beta y_i - z_i = 0.
```

参考程序独立实现 Decimal PR76 相性质，用高精度中心差分只在**参考程序内部**形成 `G_s` 和 `G_q`，然后求

```text
G_s ds/dq = -G_q.
```

它再对 p、T 和约化 feed 分别施加明显更大的有限扰动，每次重新求解完整平衡，以第二条数值路径交叉检查隐式导数。参考侧有限差分不是生产实现。

首批两个内部状态：

1. 既有 nitrogen/ethane 二元回归：`p=7.6 MPa, T=270 K, z_N2=0.30`；参数和 `kij=0.08` 沿用现有独立二元回归。
2. PR #16 已验证 methane/ethane/propane 三元中间状态：`p=4 MPa, T=220 K, betaV≈0.5`，参数与三个显式零 `kij` 完全沿用已合并三元回归。

参考头共记录 61 个导数字段，并由 Decimal(80) 脚本重新计算到 `1e-30*max(1,abs(value))`。C++ 与独立参考的比较使用单独的导数回归预算，不修改任何 flash 停止条件或物性容差。

## 7. 仍未声称的能力

本增量不能推出相包络上的连续导数、临界点导数、跨相集合导数或“全局 flash Jacobian”。有限 TPD 仍不是全局稳定性认证。后续若要用于流动 Newton Jacobian，调用方必须读取灵敏度状态并在不可用区采用明确的相切换/降阶策略，而不能继续使用上一个固定相分支导数。
