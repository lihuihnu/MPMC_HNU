# 本地与云端开发

开发规则见 [AGENTS.md](../AGENTS.md)。本页提供可执行入口，不定义另一套科学验收标准。

## 工作方式

1. 获取远程状态，检查当前分支、已有改动和适用规则。围绕一个问题写清预期行为、受影响接口和验收条件。
2. 云端开发继续使用既有 [GitHub CI](../.github/ci/README.md)。本地开发先检查环境，再按累计差异选测；不为了每轮本地试验反复推送。
3. 保留原有 CMake target、CTest、断言、参考值和容差。局部修复检查必要下游，性能改动同时检查正确性与测量条件。
4. 修改完整后检查差异、证据和未验证项，提交一个批次并推送对应 PR 分支。现有云端平台/专项验证继续执行；同步到分支与合入主线分别报告。

## 运行规模与资源分工

本地开发允许承担大网格、大算例和长时实际计算，并用于暴露规模相关问题、测量性能与峰值内存。先按本机内存/CPU 估算规模和并发，明确输入、版本、停止条件与验收判据；结果和日志放仓库外。`local.py` 的有限套件支持范围不限制使用模块原生命令开展这些运行。

云端可在资源约束下保留足以保证本次改动正确性与可行性的代表性部分，无需重复全部本地大算例。测试规模可以不同，物理假设、断言、容差、必要边界和失败语义必须一致；规模相关风险须由对应本地验证覆盖。提交报告分别列出本地实际规模/时间/峰值内存、云端执行范围及未验证项。现有 workflow 继续保留，变更其规模须单独审计，不静默削减必要覆盖。完整规则见 [AGENTS.md 的资源分工](../AGENTS.md#本地与云端的资源分工)。

## 环境准备

基础依赖为 Git、Python 3.10+、PyYAML 6.0.2、CMake 3.21+（含 CTest）及 C++20 编译器。Windows 使用 Visual Studio 2022 Build Tools 的 C++ 工具和 Windows SDK；CMake 可直接使用 Visual Studio generator，无需先把 `cl.exe` 加入普通终端 PATH。Linux/macOS 使用已有 C++ 工具链。初次采用新的编译器/配置应重新建立构建目录。

从仓库根部执行。环境、构建和日志路径由用户选择，必须放在检出目录外。以下是本项目 Windows 三目录布局的 PowerShell 示例：

```powershell
$agentDir = (Resolve-Path ../../03_agent_workspace).Path
python -m venv "$agentDir/venvs/development"
$python = "$agentDir/venvs/development/Scripts/python.exe"
& $python -m pip install -r .github/ci/requirements-local.txt
& $python -B .github/ci/local.py doctor
```

已有虚拟环境时直接设置 `$python` 并运行检查，无需重新创建。`doctor` 只检查基础工具可发现性和版本，完整可用性仍须实际 configure/build/CTest 证明。安装失败须保留错误，不自动替换科学依赖或降低测试条件。不要全局覆盖其他项目的 Python 包。

Linux/macOS 示例使用 `python3 -m venv /absolute/agent/venvs/development`，后续解释器为该环境的 `bin/python`。安装同一 requirements 后使用下面相同的 Python 入口；这里不宣称本地平台与 GitHub runner 的系统版本相同。

## 先检查范围

```powershell
git fetch origin
& $python -B .github/ci/local.py plan --base origin/main
# 只预览某条路径的影响，不执行测试：
& $python -B .github/ci/local.py plan --path modules/mesh/include/mpmc/mesh/topology.hpp
```

实际选测从 `merge-base(base, HEAD)` 比较到当前工作区，包括已提交、已暂存、未暂存、删除/重命名两侧和未忽略的未跟踪文件。工具不自动 fetch、stash、清理或修改分支；先由开发者确认 base。默认面向 main 的 PR，复用 `.github/ci/plan.py` 的累计 `ready_for_review` 选择及下游规则，不使用云端历史 checkpoint。

输出区分 `selected_cloud_gates`、`local_suites` 与 `cloud_only_gates`。未知源码/测试路径没有 owner 时直接失败，须先完善既有 CI ownership。治理检查总是包含；纯治理文档不会因此触发全部科学测试。中央 workflow 本身修改时，云端还会比较 job 语义；本地不会仅凭路径声明该比较完成，而是额外标记 `central_router_semantic_validation` 为待云端验证。

## 执行与结果

```powershell
# 运行当前累计改动中可在本地执行的测试；未支持的选中 Gate 明确报告为 partial。
& $python -B .github/ci/local.py run --base origin/main `
  --build-root "$agentDir/b" --config Release --jobs 2

# 聚焦开发时可明确选择原有套件；这不代表全部受影响测试已通过。
& $python -B .github/ci/local.py run --suite ad.math --suite thermo.pr76 `
  --build-root "$agentDir/b" --config Release --jobs 2
```

| 本地 suite | 复用内容 |
| --- | --- |
| `governance` | 原 `verify_workflows.py`，含选测、注册 CTest、工作流映射和本地入口回归 |
| `ad.arithmetic` | 原 AD 算术、standalone consumer、全部 64 种根可选库组合；此套件成本明显高于单个 AD 测试 |
| `ad.math`、`ad.jacobian`、`ad.runtime` | 各自原有 AD target 与 CTest 筛选 |
| `thermo.contracts`、`thermo.pr76`、`thermo.pr76_mixture`、`thermo.pr76_pt` | 原热力学套件；PT 同时重新计算既有独立 Decimal 参考 |
| `legacy_mesh_core`、`legacy_discretization_core` | 原 Mesh/Discretization 工作流中的 source、target 列表与 CTest 筛选；保留既有 Gate 名以便对照云端 |
| `flow_core`、`flow_discretization` | 原 Flow 工作流中的 source、target 列表与 CTest 筛选 |

目前本地入口执行 native Debug/Release、sanitizer OFF。不会把此配置冒充 GCC ASan/UBSan；Flow 的 GCC 专项独立 oracle 步骤仍由原云端 Gate 拥有。四个核心套件直接读取原 workflow 中已审计的构建/测试字段；格式变化无法明确解析时拒绝运行，要求审查本地适配，不执行任意工作流 shell。AD/热力学保留原套件映射及其额外消费/参考步骤。

构建目录按配置和套件隔离，使用稳定的短哈希子目录减少 Windows/MSBuild 长路径问题，默认并发 2。Windows 应优先使用示例中的短 `b` 目录；若外层路径仍很长，须选用更短的仓库外路径。每次运行在 `--build-root/runs/` 下生成唯一目录，记录命令、每套件日志、HEAD、工作区差异摘要、环境、结果和未支持 Gate。临时配置目录也放在 build root 下。测试期间不要编辑仓库；检测到源码或工作区内容变化时，结果标为失败，不能归到稳定版本。

退出码：`0` 为所声明本地范围通过，`1` 为配置/构建/测试/来源状态失败，`2` 为已执行的本地范围通过但仍有选中的 Gate 未支持。`plan` 成功只证明已生成计划。`--suite` 是明确缩小的开发验证范围，不能用它绕过正式验收中的必要下游。

失败后查看该运行日志，集中修复已确认问题，再执行受影响范围；不通过删测试、放宽容差或把失败转换成 skipped 获得成功。

## 与现有云端流程的边界

现有 63 个 workflow、触发、平台矩阵、测试命令、私有 runner 白名单、发布及签名语义保留。本地入口的回归由原中央 impact job 调用，不新增独立 CI 入口。

PETSc/MPI、外部格式样本、SW92/CPA/flash 专项、Clapeyron/ThermoPack、前端、gRPC、Android 和安装/发行等尚未接入本地适配。选中后列为 cloud-only，继续使用其既有入口；不能据此说项目没有这些测试。按实际开发需要逐项安装、固定和验证依赖后，才能扩展本地支持。基础工具安装不代表已具备所有专项环境。

所有本地结果还需按改动风险补足其他平台、sanitizer 和专项证据。本次流程不修复既有云端 checkpoint 复用偏差，也不更改远程分支保护。遇到该偏差应明确记录累计重验的具体提交和运行，不把人工补验当作调度已修复。
