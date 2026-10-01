# 根入口与库消费

项目根的 `MPMC_ENABLE_THERMODYNAMICS`、`MPMC_ENABLE_FLASH`、`MPMC_ENABLE_MESH`、`MPMC_ENABLE_DISCRETIZATION` 默认均为 `OFF`。AD 始终可用；启用 flash 会按既有 target 依赖带入 thermodynamics 和 AD，启用 discretization 会带入 mesh。两条依赖链可以独立或同时启用，不强制改写下层选项的缓存值。

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

从项目根配置库目标：

```sh
cmake --preset flash-libraries
cmake --build --preset flash-libraries
# 仅热力学：
cmake -S . -B build/thermo-libraries -DBUILD_TESTING=OFF -DMPMC_ENABLE_THERMODYNAMICS=ON -DMPMC_ENABLE_FLASH=OFF
# 网格与离散库（不启用 PETSc）：
cmake -S . -B build/spatial-libraries -DBUILD_TESTING=OFF -DMPMC_ENABLE_DISCRETIZATION=ON
```

这些库当前是 header-only `INTERFACE` targets；配置/构建成功本身不代表编译了全部公共头或通过科学验证。`flash-libraries` preset 显式关闭空间库，不注册测试，也不提供 test preset；原有 `ad-debug` preset 保持不变。普通根构建中 `BUILD_TESTING=ON` 仍只注册 `ad.dual`，专项数值/物理回归继续使用[热力学](../../../modules/thermodynamics/README.md)、[flash](../../../modules/flash/README.md)、[mesh](../../../modules/mesh/README.md) 文档及 [discretization 注册入口](../../../tests/discretization/core/CMakeLists.txt) 中的独立工程。

嵌入其他 CMake 工程时，在 `add_subdirectory` 前设置所需选项；库消费者只链接所需最高层 target，不手工添加下层 include 路径：

```cmake
set(BUILD_TESTING OFF CACHE BOOL "Disable MPMC root tests")
set(MPMC_ENABLE_FLASH ON CACHE BOOL "Enable MPMC flash library")
add_subdirectory(path/to/MPMC_HNU mpmc)
target_link_libraries(my_solver PRIVATE mpmc::flash)
```

缓存变量由调用方管理；若现有构建目录已有值，使用 `-D` 显式覆盖。`BUILD_TESTING` 是共享 CMake 选项，宿主若需自己的测试可自行 `enable_testing()`。此入口不引入 MPI/PETSc、RPC、UI、网络下载或新的第三方依赖，也不改变各模块直接 `add_subdirectory` 的现有用法。

[消费方契约](CMakeLists.txt) 检查默认 AD 测试清单、四个开关的全部 16 种组合、重新关闭选项、库目标可见性和传递使用要求，并执行消费程序。`build.root_libraries` 继续消费 AD/thermodynamics/flash；启用空间库时，独立的 `build.root_spatial_libraries` 只链接 mesh 或 discretization，不借用 AD/热力学的 include 路径。它执行稳定 ID 与 synthetic 零传导系数探针，不重新拥有 mesh/discretization 的科学回归。唯一 CI owner 是既有 AD arithmetic 的 GCC Debug + ASan/UBSan、Clang Release、MSVC Release 矩阵；这是构建集成证据，不是新增物理模型验证。根构建、相关模块 CMake、消费程序实际依赖的公共头和本节命令变更均进入该 owner；现有科学测试仍按各自规则选测。完整复现：

```sh
python3 tests/build/root_libraries/verify.py --config Debug --compiler g++ --sanitizer ON
```
