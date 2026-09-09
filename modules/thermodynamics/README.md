# 热力学模块：有序组分与 PR 数据契约

本页描述 C++20 数据契约、校验与有序快照。后续已新增[PR76 纯组分 a(T)、b 数值核](pr76_pure.md)及解析/AD 温度导数测试，并增加[运行期经典混合参数](pr76_mixture.md)与完整/约化组成导数；已增加 [PT 候选相性质核](pr76_phase.md)：给定 p、T、相组成计算可用 Z 根、逐组分 ln(phi) 与简单根的局部隐式导数；尚无相稳定性与闪蒸。CMake 目标 `mpmc::thermodynamics` 仅依赖标准库；不依赖 AD、core、网络、数据库或第三方框架。所有数值输入仍须有独立科学依据，构造成功不等于物理模型已验证。

## 1. 模型范围与资料核验

本轮唯一接受的模型标识为：

```text
PR76/classical-vdw/constant-kij/no-translation
```

该标识固定为 1976 年 Peng–Robinson 原始 alpha 形式，经典二次吸引项/线性共体积混合规则，常数、对称的二元交互参数，且不做体积平移。[R1][R2] 不使用 PR78 的高偏心因子分段修正，也不自动选择其他 alpha、SW、CPA、PRSV 或温度相关 kij。模型 ID 和参数数据集的 `dataset_id/revision` 是不同的身份，不能相互代替。

审阅了 Peng 作者上传原文的网页转录中模型参数及混合规则所在部分（含式 18、20–22）、NIST teqp 官方参数接口及 whitson 官方说明。此前只有网页转录；用户随后提供了清晰的 PR76 PDF。本轮已对照渲染页核验期刊第 60 页的纯组分式 (9)、(10)、(12)、(13)、(17)、(18)，采用印刷系数与单独声明的现代 SI 气体常数。文件摘要、公式及核验范围见[数值核说明](pr76_pure.md)；不将这次局部核验说成全文全部表格/勘误或真实物性已验证。随后 PT 增量对照原文第 60 页式 (4)–(8)、(19) 核验根与逸度表达式，明确采用精确 1±sqrt(2) 而不是原文印刷的 2.414/0.414；数值约定见 [PT 说明](pr76_phase.md)。本次没有录入任何真实组分参数表。

SW 是 Søreide–Whitson 含水体系的 PR 修改，不是所有“改进 PR”的统称。whitson 官方水物性说明区分水相/非水相参数，并涉及温度、盐度及水 alpha 修改；该现代页面同时引用 Yan 等后续关联式，不能整体当作 1992 原版 SW 参数表。[R3] 后续应明确选择原始 SW 还是有单独版本号的扩展。

CPA 不能简单共享 PR 参数数值。已审阅 DTU 作者稿第 2 节经典 CPA 组成与参数说明：该实现以 SRK 项加缔合项为基础，需要 a0、b、c1，以及缔合能、缔合体积和位点方案。[R4] 这篇文章还研究 crossover 扩展；本项目没有据此选用该扩展。未来 CPA 的立方项、位点/交叉规则及能量单位均需单独核验；不能把 PR 由 Tc/Pc 得到的 a/b 自动当作 CPA 拟合参数。

## 2. 公共数据层：components.hpp

`Component`、`Provenance`、`SourcedScalar`、`Applicability` 是**未验证的输入记录**，允许表达缺失；经过 `OrderedComponents::select` 或 `PrParameterSet::create` 校验后才形成可消费快照。没有“填写字段即科学正确”的隐含保证。

| 字段/类型 | 规则 |
| --- | --- |
| `Component.id` | 非空、无空白的可见 ASCII，大小写敏感。是稳定身份键，不是数组下标、显示名称或自动解析的 CAS。 |
| `display_name` | 非空白标签，可使用 UTF-8，可以重名；不参与参数匹配。 |
| `kind` | 必须显式为 `pure` 或 `pseudo`；不从名称猜测纯物质、假组分、水或相态。 |
| `definition` | 必需来源记录；假组分应指向具体表征方案与版本。同名 C7+ 不保证为同一组分。 |
| `molar_mass` | 可选，提供时必须有限且 >0，单位 `kilogram_per_mole`。PR76 数据契约不强制它存在。 |
| `OrderedComponents` | 从 catalog 按请求 ID 顺序创建拥有数据的快照，提供 `size/at/index_of/items`，不暴露可变元素。 |

目录和请求顺序均拒绝重复 ID；不存在的 ID、非法 ID 不会被修剪、转小写、去重或用其他组分替代。会校验整个传入目录，包括未选择的组分；目录可以包含暂时不用于本次请求的有效组分。ID 的化学真实性/别名等价性由数据提供者负责，本轮不创建全球组分注册库。

`Provenance` 必须记录 `kind/reference/revision/locator/acquisition/usage_terms`，分别为来源类别、文献或文件身份、版本、定位、取得方式和许可/使用声明。`assumption` 与 `synthetic_test` 还必须有非空白说明 `note`。文献 DOI、数据库行号、本地数据内容摘要等可以作为信息，不会在构造期间联网解析或验证真实性、权限或内容摘要。

`SourcedScalar` 还必须有存储单位 `unit`、`original_unit` 和 `conversion`。数值**已在调用前转换到规定单位**；例如 MPa 输入由导入层转换到 Pa 后，记录原单位与转换步骤。这里校验目标单位标签和字段域，不解析单位字符串、不执行隐式换算，也不能检测“错误数值被标记为正确单位”。同精度 `double` 数组只是当前参数存储，不是强制未来 EOS 状态使用 double；未来参数拟合求导需单独设计。

来源类别包括文献、数据库、用户提供、明确假设和人工测试。默认 `DataPolicy::ordinary` 拒绝 `synthetic_test`；测试须明确传入 `allow_synthetic_tests`，且原类别不会被擦除。这个开关是防误用措施，不是安全授权系统；用户提供或文献来源也不会自动获得“独立验证通过”状态。许可字段完整不等于已获重新分发许可。

## 3. PR 专属数据层：pr_parameters.hpp

| 输入字段 | 必需性、单位及数值规则 |
| --- | --- |
| `model_id` | 必须精确匹配上面的 PR76 配置；裸 `PR`、空值、SW/CPA 均拒绝。 |
| `dataset_id/revision` | 两项均须非空白；参数发生变化时提供者必须更新版本。 |
| `applicability.declaration` | 必需；记录整个数据集的适用性来源及局限。 |
| `PrPureRecord.component_id` | 按 ID 关联目录；每个被选择的组分必须恰好有一个纯组分参数记录。 |
| `critical_temperature` | 必需 `SourcedScalar`；K，有限且 >0。 |
| `critical_pressure` | 必需 `SourcedScalar`；Pa，有限且 >0。 |
| `acentric_factor` | 必需 `SourcedScalar`；无量纲，有限；负值与零允许，不擅自设置通用经验上下限。 |
| `PrBinaryRecord.first_id/second_id` | 不同的已知组分 ID；当前 PR 配置使用无序组分对，两方向重复也拒绝。 |
| `kij` | 每个被选择的异组分对必需，有限、无量纲、可为零/负数；不擅自要求 [0,1]。 |

使用 `optional` 区分缺失与合法零值；显式数据对象的默认数值为 NaN 而非零。缺少交互参数不会使用零、上一配置参数或其他 EOS 的参数。采用零 kij 假设时仍须提交有来源、版本和说明的显式记录。

**对角线 kij=0 是此 PR 混合配置的结构约定，不是测得或拟合的参数**，故不接受自配对记录。成功结果在所选异组分对上没有缺项。记录可以覆盖比当前请求更大的组分集合；但已提交的任何记录若重复、ID 未知、单位或值非法，都会报错，而不是因暂未选择就忽略。缺少未选择组分的整条参数记录允许；选择它时再检查完整性。

生成的 `PrParameterSet` 保存有序 `components`、有序 `pure_records`、上三角遍历顺序的 `binary_records`（含原始来源），以及连续只读数组 `critical_temperatures_k/critical_pressures_pa/acentric_factors/kij_matrix`。矩阵索引 `i*n+j`，两轴均遵循请求顺序。`kij(i,j)` 检查两轴；只读 span 的 `operator[]` 不额外检查下标，使用者必须遵守长度。

重新排序不要求调用方手动排列 Tc/Pc/omega 或 kij 两轴。新增/删除/替换组分都通过新构造完成，不提供容易造成半更新的 `resize_components` 或独立矩阵修改器。原快照不受影响；构造失败不发布部分结果。快照可拷贝/移动构造，但禁止赋值，避免多字段复制失败后留下失配状态；移动后的源对象仅应销毁，不继续用于计算。视图只能从左值取得，寿命不超过所属快照；不承诺跨线程修改源记录安全。

初始化使用基于 ID 的有序查找表，随后把数值集中为连续数组；字符串/来源处理不放入未来 EOS 内循环。储存全矩阵需 O(n²)，每对绑定包含对数级查找。构造及拥有数据的快照会分配内存；没有性能基准，也未实现缓存、二进制 ABI 或序列化。

## 4. 适用范围不是布尔“科学有效”

`Applicability` 的 T/p 区间可选，给出时要求有限、正值、下界不大于上界。缺省表示未知，不表示无限范围，也不以临界点自动推断范围。

`assess(T, p)` 只返回 `inside_declared_bounds`、`outside_declared_bounds` 或 `unknown`：任一已知轴越界即 outside；两轴都有界且都在区间内才为 inside；其余为 unknown。端点包含在内，非法状态数值报错。这个函数不计算 EOS，也不证明相态、组成、盐度或临界区准确性；声明来源的 note/locator 应解释其他约束。

当前在数据集层记录整个集合的保守适用范围，不自动交叉合并每条物性来源的温压域，也不因选择了子集就扩大范围。未来存在不同相别/关联式/拟合域时应按模型细化，不伪造缺失区间。数值校验、声明范围内与独立实验验证是三个不同层次。

`ContractLimits` 可约束目录/所选组分数、输入 pair 记录数和生成矩阵元素数；构造前先除后乘检查 n*n 的表示/容器上限。默认只有表示限制，**不是服务端资源预算**，导入层还须限制文件大小、字符串字节和总记录数。错误返回 `ContractError`（派生自 `invalid_argument`），含稳定 `code()`、字段路径和原因；索引越界为 `out_of_range`，内存分配异常原样传播。不捕获错误后伪装成空参数集。

## 5. SW 与 CPA 的复用边界

| 共用部分 | 模型专属部分，不能误复用 |
| --- | --- |
| 组分 ID、纯/假组分定义、运行期有序选择、可选摩尔质量 | SW 水组分身份与相别能力必须显式确认，不能由排列位置猜测。 |
| 来源/版本/许可/换算记录、缺失与零的区分、错误诊断 | SW 的相别/温度/盐度相关 kij 和水 alpha 应有自己的关联式记录及版本。 |
| 不可半更新快照、稳定的 ID 绑定与只读数据视图 | CPA 的立方参数、能量基准、位点数量/类别、交叉缔合规则及参数另行定义。 |
| 声明适用范围的表达方式 | 每个模型的数据集、拟合范围与验证证据独立；PR 的 symmetric constant-kij 规则不进入公共组分层。 |

SW 可以在来源和定义兼容的前提下复用 PR 纯组分记录的字段表示；不会通过把 `model_id` 改名就将 `PrParameterSet` 当作 SW。CPA 不被要求伪造它本不需要的 PR 临界属性或从 omega 推导拟合参数。暂不创建 SW/CPA 空类或参数字典；真实需求到来时允许把已证实共用的绑定实现提取出来，不预设继承体系。

前端增减组分时，应提交有序 ID 与完整新配置，并建立新快照。即使数量相同，替换 ID、修改假组分定义或参数版本仍是新配置；未来缓存键应覆盖有序身份与全部相关版本，不能只用 n。摩尔组成状态、独立变量坐标、任务版本和缓存失效不属于本增量。

## 6. 使用与必要增量测试

包含 `<mpmc/thermodynamics/components.hpp>` 或 `<mpmc/thermodynamics/pr_parameters.hpp>`；外部项目用 `add_subdirectory(path/to/modules/thermodynamics thermo)` 并链接 `mpmc::thermodynamics`。输入记录应由调用方从经审计资料构建，最小调用为：

```cpp
// catalog: vector<Component>; order: vector<string>; input: PrParameterInput.
// All records must already contain explicit SI units and provenance.
const auto parameters = mpmc::thermodynamics::PrParameterSet::create(
    catalog, order, input);
const auto temperatures = parameters.critical_temperatures_k();
const double interaction = parameters.kij(0, 1); // Requires at least two selected components.
```

这是依赖调用方已有输入的片段，不是有物理数据的运行示例。完整的人工数据构造及默认拒绝测试见 `tests/thermodynamics/contracts/contracts_test.cpp`。

```sh
cmake -S tests/thermodynamics/contracts -B build/thermo-contracts -DCMAKE_BUILD_TYPE=Debug
cmake --build build/thermo-contracts --target mpmc_thermodynamics_contract_tests --config Debug
ctest --test-dir build/thermo-contracts -C Debug -R '^thermo[.]contracts[.]' --output-on-failure --no-tests=error
```

新增 13 个独立 CTest 条目：有序身份、增减替换、4 组分全部 24 种排列、身份错误、缺失参数、重复/非法 pair、值/单位、来源策略、模型/数据版本、范围声明、所有权/失败恢复、限额、公共头。参考为手工指定的人工 ID 映射与矩阵，不从被测输出生成期望值；无需浮点近似容差，因为这里只做校验和数据搬运。两个公共头分别以首个/重复包含编译，再与主测试跨翻译单元链接；Release 使用正常检查而非可禁用 assert。

正式执行仅用 GitHub 官方托管 runner：工作流 `Thermodynamics contracts` 已增加按依赖选择 contracts/pr76/pr76_mixture/pr76_pt 的路由；独立配置所选工程，不构建或运行无关旧套件。新增公共接口与所有权需要必要 GCC/Clang/MSVC 编译覆盖及 Linux ASan/UBSan。仅文档改动不启动该工作流；启用必需检查前仍须设计始终回报状态的外层门禁，本轮不修改分支保护。运行结果以实际提交的 Actions 日志为准。

原契约增量没有进行真实物性验证、AD 导数集成、跨 EOS 参数转移验证、序列化/前端测试、并发竞态、OOM 注入、macOS/其他架构或性能基准。旧 AD 代码、测试、构建、工作流及 AGENTS 均保持不变。

## 7. 参考与待取得资料

本次独立实现，不复制成熟软件源码、参数表或论文全文，不新增第三方依赖或项目许可证。

- [R1] Peng, D.-Y.; Robinson, D. B. (1976). *A New Two-Constant Equation of State*, 15(1), 59–64. [DOI:10.1021/i160057a011](https://doi.org/10.1021/i160057a011)。[作者上传的全文转录入口](https://www.researchgate.net/publication/231293953_New_Two-Constant_Equation_of_State)。原契约增量审阅转录；随后取得用户提供 PDF 并核验纯组分公式，定位与 SHA256 见 [pr76_pure.md](pr76_pure.md)。未复核表格数值或全部勘误。
- [R2] [NIST teqp — General cubics](https://pages.nist.gov/teqp-docs/en/main/models/cubics.html)。核对 Tc/Pc/omega 输入、SI 示例和模型/参数分离。其接口也允许非对称矩阵；本契约选择对称常数 kij 是明确限定，不宣称所有立方模型都要求如此。
- [R3] [whitson 官方 Water bot 说明](https://manual.whitson.com/methods/water-bot/)。用于核验相别/温度/盐度相关扩展边界，未移植表中的数值。原始 SW：Søreide & Whitson (1992), [DOI:10.1016/0378-3812(92)85105-H](https://doi.org/10.1016/0378-3812(92)85105-H)，原文尚需补齐。
- [R4] Vinhal, Yan & Kontogeorgis (2020), *Modeling the Critical and Phase Equilibrium Properties of Pure Fluids and Mixtures with the Crossover Cubic-Plus-Association Equation of State*, [DOI:10.1021/acs.jced.9b00492](https://doi.org/10.1021/acs.jced.9b00492)。[DTU 作者稿](https://backend.orbit.dtu.dk/ws/files/199155034/Vinhal_et_al_Article_New_clean.pdf)，已查看 PDF 第 3 页（文内第 2 页）经典 CPA 参数说明；不将 crossover 扩展选为项目 CPA 版本。原始 CPA 1996 年论文 [DOI:10.1021/ie9600203](https://doi.org/10.1021/ie9600203) 的全部公式与参数表未核验。

纯组分实现所需 PR76 原文已由用户提供并核验；实际组分应用前仍需要 Tc/Pc/omega、所选 pair 的 kij 与可追溯验证资料。缺乏数据只阻塞依赖它的数据录入/物理验证，不能由合成测试数据替代。SW/CPA 的原文、位点/关联式方案与参数版本在各自实现前单独补齐。
