# 根入口与库消费

项目根的 `MPMC_ENABLE_THERMODYNAMICS`、`MPMC_ENABLE_FLASH`、`MPMC_ENABLE_MESH`、`MPMC_ENABLE_DISCRETIZATION`、`MPMC_ENABLE_FLOW_DISCRETIZATION` 默认均为 `OFF`。AD 始终可用；启用 flash 会按既有 target 依赖带入 thermodynamics 和 AD，启用 discretization 会带入 mesh。依赖链可以独立或同时启用，不强制改写下层选项的缓存值。

| thermodynamics 选项 | flash 选项 | 根入口提供的库目标 |
| --- | --- | --- |
| OFF | OFF | `mpmc::ad` |
| ON | OFF | `mpmc::ad`、`mpmc::thermodynamics` |
| OFF 或 ON | ON | `mpmc::ad`、`mpmc::thermodynamics`、`mpmc::flash` |

| mesh 选项 | discretization 选项 | 额外提供的空间库目标 |
| --- | --- | --- |
| OFF | OFF | 无 |
| ON | OFF | `mpmc::mesh` |
| OFF 或 ON | ON | `mpmc::mesh`、`mpmc::discretization` |

启用 `MPMC_ENABLE_FLOW_DISCRETIZATION` 会提供 `mpmc::flow_discretization`、`mpmc::flow` 和 `mpmc::flow_thermodynamics`，并带入 AD、thermodynamics、mesh、discretization；不启用 flash 或 PETSc。即使这些下层选项仍为 `OFF`，必要库目标也会存在。已有 flow 模块同时定义两个消费契约：`mpmc::flow_thermodynamics` 传递 flow/thermodynamics/AD 使用要求；`mpmc::flow_discretization` 传递 flow/discretization/mesh 使用要求，本身不向消费者传递热力学头文件路径。

从项目根配置库目标：

```sh
cmake --preset flash-libraries
cmake --build --preset flash-libraries
# 仅热力学：
cmake -S . -B build/thermo-libraries -DBUILD_TESTING=OFF -DMPMC_ENABLE_THERMODYNAMICS=ON -DMPMC_ENABLE_FLASH=OFF
# 网格与离散库（不启用 PETSc）：
cmake -S . -B build/spatial-libraries -DBUILD_TESTING=OFF -DMPMC_ENABLE_DISCRETIZATION=ON
# flow 与离散桥接（仍不启用 PETSc）：
cmake -S . -B build/flow-libraries -DBUILD_TESTING=OFF -DMPMC_ENABLE_FLOW_DISCRETIZATION=ON
```

这些库当前是 header-only `INTERFACE` targets；配置/构建成功本身不代表编译了全部公共头或通过科学验证。`flash-libraries` preset 显式关闭空间库和 flow 桥接，不注册测试，也不提供 test preset；原有 `ad-debug` preset 保持不变。普通根构建中 `BUILD_TESTING=ON` 仍只注册 `ad.dual`，专项数值/物理回归继续使用[热力学](../../../modules/thermodynamics/README.md)、[flash](../../../modules/flash/README.md)、[mesh](../../../modules/mesh/README.md)、[flow](../../../modules/flow/README.md) 文档及 [discretization 注册入口](../../../tests/discretization/core/CMakeLists.txt) 中的独立工程。

嵌入其他 CMake 工程时，在 `add_subdirectory` 前设置所需选项；库消费者只链接所需最高层 target，不手工添加下层 include 路径：

```cmake
set(BUILD_TESTING OFF CACHE BOOL "Disable MPMC root tests")
set(MPMC_ENABLE_FLASH ON CACHE BOOL "Enable MPMC flash library")
add_subdirectory(path/to/MPMC_HNU mpmc)
target_link_libraries(my_solver PRIVATE mpmc::flash)
```

缓存变量由调用方管理；若现有构建目录已有值，使用 `-D` 显式覆盖。`BUILD_TESTING` 是共享 CMake 选项，宿主若需自己的测试可自行 `enable_testing()`。此入口不引入 MPI/PETSc、RPC、UI、网络下载或新的第三方依赖，也不改变各模块直接 `add_subdirectory` 的现有用法。

[消费方契约](CMakeLists.txt) 检查默认 AD 测试清单、五个开关的全部 32 种组合、重新关闭选项、库目标可见性和传递使用要求，并执行消费程序。原 16 种组合完整保留为 flow 桥接关闭时的子集。各消费者的头文件路径分别由自己链接的 target 提供：

| CTest | 启用条件 | 唯一链接的生产 target |
| --- | --- | --- |
| `build.root_libraries` | 始终 | 显式选项中的最高层 `mpmc::flash` / `mpmc::thermodynamics` / `mpmc::ad` |
| `build.root_spatial_libraries` | 显式启用 mesh 或 discretization | `mpmc::mesh` 或 `mpmc::discretization` |
| `build.root_flow_thermodynamics` | 启用 flow_discretization | `mpmc::flow_thermodynamics` |
| `build.root_flow_discretization` | 启用 flow_discretization | `mpmc::flow_discretization` |

这些用例只执行稳定 ID、自然变量布局、来源无关的类型使用与 synthetic 零传导系数探针，不重新拥有模块科学回归。唯一 CI owner 是既有 AD arithmetic 的 GCC Debug + ASan/UBSan、Clang Release、MSVC Release 矩阵；这是构建集成证据，不是新增物理模型验证。根构建、相关模块 CMake、消费程序实际依赖的公共头和本节命令变更均进入该 owner；现有科学测试仍按各自规则选测。完整复现：

```sh
python3 tests/build/root_libraries/verify.py --config Debug --compiler g++ --sanitizer ON
```
