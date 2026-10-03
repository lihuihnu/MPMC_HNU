# 外部网格兼容性与二维/三维独立读回

唯一 C++ owner 是本目录的 `mpmc_mesh_external_compatibility` target。原有五路径参数入口继续验证固定上游 deal.II/OPM 样例；新增 `--emit-2d <目录>` 只生成待检文件。`verify_2d_readers.py` 调用这个 producer，再用官方 Gmsh 4.15.2 API 和 VTK 9.7.1 `vtkXMLUnstructuredGridReader` 读取输出，完全不调用 MPMC 解析器。它不是新增 CTest，也不由 mesh core target 重复编译或执行。

## 判据与边界

输入为人工构造的 SI 平面网格：独立三角形、顺时针非矩形凸四边形、三角形/四边形混合共享边。包括稀疏且乱序的顶点/单元 ID、显式/生成 face ID、边界组、重叠区域组、点与单元的标量和多分量字段。每例输出原生 MSH/VTU 及 canonical Gmsh→VTU、VTU→Gmsh，共 12 个文件、6 份转换报告。

- Python 中独立写出输入和解析期望，不读取 producer 生成的期望值。三角形面积/质心为 `1 m², (2/3,1/3) m`；梯形为 `5 m², (19/15,14/15) m`；混合例额外三角形为 `1 m², (3,2/3) m`。从第三方读到的坐标以三角形扇分解重新积分，并与解析值比较。
- Gmsh 验证节点坐标、稳定节点/单元/面 ID、线/三角形/四边形类型、循环连接、共享边、PhysicalNames、边界与多区域组成员。
- VTK 验证坐标、三角形/四边形类型、循环连接、UInt64 单元 ID、PointData/CellData 的字段名、分量数、数值和实体对应。ID 使用整数访问器，避免经过 double。VTU 用保留数组 `mpmc_global_vertex_id` 编码稳定节点 ID；官方读取器按 UInt64 整数逐点核验，并通过 FieldData 面身份表核验边/面绑定。
- 数值采用绝对/相对容差 `1e-12`（坐标/质心 m、面积 m²、字段按声明单位）；样例量级为 1–300，容差用于十进制文本和浮点积分舍入。类型、整数 ID、连接和分组精确匹配。
- 跨格式输出同时检查实际可读内容与 loss report：Gmsh 明确丢失字段，VTU 明确丢失组和面标签。节点和 face ID 均已保留；报告不再含 `vtu.face_ids_remapped` 或 `vtu.vertex_ids_remapped`，官方 VTK 同时确认实际实体绑定和文件点顺序。无标签/字段、默认编号不误报及同集合错绑的实际往返回归由既有 `mesh.core.exchange_io` 唯一拥有，覆盖二维混合单元与三维金字塔的 16 个场景。VTU 的自定义 `mpmc_*` 单位/来源 XML 属性不是 VTK 标准字段语义，本测试不宣称 VTK 理解这些元数据，也不证明所有转换损失都已穷尽。
- 五个负对照分别移除单元 ID 字段、改错字段值、反转连接、改错组名及破坏 XML；均须被拒绝。破坏 XML 的控制会产生预期的 VTK parser 错误诊断，只有所有拒绝判据满足才打印总 PASS。

范围为现有二维线性 XY 平面、MSH 4.1 ASCII 和单 Piece VTU ASCII。三维范围见下节；未覆盖 binary/appended/compressed、高阶/曲面、大规模性能、Gmsh GUI/CAD 操作或科学求解验证。旧公共真实样例验证仍保留，不能与本项混称。

## 三维四类线性单元的独立读回

同一 C++ target 的 `emit_3d_exports.cpp` 由 `--emit-3d <目录>` 调用。`verify_3d_readers.py` 独立定义以下人工解析输入；每个文件放两个同类、互不连接的单元，第二个由第一个缩放 2 倍并平移 `(10,-5,2) m` 得到。单元文件顺序与点块顺序相反，ID 为 `31/901`，顶点采用乱序稀疏 ID，两个基面保留显式 `701/702`，其余面按契约生成。这样的输入能检查字段/标签是否对应到了正确实体。

| 单元 | 基础几何（m） | 基础/变换单元解析体积（m³） |
| --- | --- | --- |
| 四面体 | 从原点沿三个坐标轴伸出 2、3、4 | 4 / 32 |
| 六面体 | 2 × 3 × 4 长方体 | 24 / 192 |
| 三棱柱 | 两直角边 2、3，柱高 4 | 12 / 96 |
| 金字塔 | 2 × 3 矩形底，底中心正上方高度 4 | 8 / 64 |

四类各输出原生 MSH/VTU 与两个 canonical 适配方向，共 **16 文件、8 份损失报告**。这些是从人工构造的格式输入对象经 writer/adapter 生成的文件，不能混称为先读取这些文件再跨格式转换的额外输入链。

- 官方 Gmsh 检查节点/体单元/面 ID、单元类型和节点顺序、三角/四边面连接、边界和多体区域 Physical Groups 的名称及成员。利用官方 `getIntegrationPoints(type, "Gauss2")` 与 `getJacobian(element, points)` 对 Jacobian 行列式积分，每个积分点须为正，所得体积与上述解析值比较；不读取或复用 MPMC 的几何计算结果。
- 官方 VTK 检查 point 顺序/坐标、单元类型/节点顺序、UInt64 单元 ID、PointData 温度标量和三分量位置（m），以及 CellData 标量/双分量字段。`vtkCellSizeFilter` 从读回单元独立计算体积并逐单元对照解析值。
- 数值容差沿用绝对/相对 `1e-12`，分别用于坐标（m）、体积（m³）、字段（声明单位），适用于当前解析样例的文本/浮点舍入；类型、整数 ID、连接及分组精确比较。规则依据：[Gmsh 积分与 Jacobian API](https://gmsh.info/doc/texinfo/gmsh.html)、[VTK CellSizeFilter](https://vtk.org/doc/nightly/html/classvtkCellSizeFilter.html)。这不是扭曲单元任意积分精度的证明。
- canonical Gmsh→VTU 必须报告组、面标签及 face 身份损失；VTU→Gmsh 必须报告字段未写入。VTU 输入对象不携带其格式尚未编码的边界标签。仍不声称 VTK 理解自定义单位/来源属性。
- **6 个负对照**：改错单元类型、反转局部连接、交换两个单元的字段值、移除 ID、保留分组名称但交换面成员、破坏 XML。每个都必须被读取检查拒绝，拒绝原因记入结果；损坏 XML 的解析错误是预期控制。

同一入口还验证以下两种三维混合单元共享面；两例的原生 MSH/VTU 与两个 canonical 方向合计增加 **8 文件、4 报告**，不调用 MPMC 读取器或几何 helper。

| 混合样例 | 共享面 | 唯一面 / 内部面 / 边界面 | 单元体积 / 总体积（m³） |
| --- | --- | --- | --- |
| 四面体 + 三棱柱 | z=4 m 的三角面，面积 3 m² | 8 / 1 / 7 | 4、12 / 16 |
| 六面体 + 金字塔 | z=4 m 的四边面，面积 6 m² | 10 / 1 / 9 | 24、8 / 32 |

- 下部柱体仍采用上表几何，上部四面体顶点为 `(0,0,8) m`，金字塔顶点为 `(1,1.5,8) m`。上、下单元分别为 31/901；共享面显式 ID 703 无边界标签，边界面 701/702 保留不同组，上部单元属于两个区域组。点/单元字段与身份损失报告沿用已有逐实体检查。
- Gmsh 从官方读取的体单元连接与独立列出的线性参考面循环构造面关联；VTK 直接枚举读回的 `vtkCell.GetFace`。两者都按实际节点身份建立键，**禁止按相同坐标合并节点**。精确检查共享面节点、两个相邻单元、每个面的关联数（1 或 2）和外边界数量；Gmsh 的 Physical Group 成员必须落在重建的真实外边界，共享面只序列化一次。
- 由读回坐标计算平面多边形面积向量，以凸单元内部的顶点均值确定朝外方向；共享面的两侧单位法向分别为 `(0,0,-1)` 和 `(0,0,1)`。面积、法向、逐单元及总体积沿用 `1e-12` 绝对/相对容差，Gmsh 正 Jacobian 与 VTK 体积检查独立保留；这不是验证文件中存在法向字段。
- 增加 **4 个负对照**（每个混合样例 2 个）：VTU 复制接口节点并只重接上部单元，坐标/字段不变且官方 VTK 仍测得原解析体积，必须因缺少共享面而失败；MSH 仅把内部面加入现有边界组，必须因内部面误标边界而失败。拒绝原因明确断言，不能让无关解析错误冒充拓扑检查。

三维范围包括上述同类双单元与两种共形混合单元、平面面的解析样例；非共形接口、非平面面、任意畸变/退化、高阶、binary 和科学求解不在本项范围。生产网格库无需新依赖。

使用下一节同一隔离环境与 producer，把脚本替换为 `verify_3d_readers.py` 并指定新的输出目录即可复现。总 PASS 应为 `files=24 reports=12 negative_controls=10`；`result.json` 还保存各文件逐单元的实际体积、8 个混合输出的面关联/面积/法向、断开接口负对照仍保持的体积、共享辅助脚本 hash、负对照拒绝原因，以及包括负对照文件在内的 SHA256。云端在原 external job 中执行一次，单独保存 `mesh-3d-independent-readback` artifact；二维命令及 artifact 均保留。

## 官方生成、MPMC 导入和导出的三维完整链路

`verify_3d_chains.py` 复用上述两种混合解析几何，由官方 Gmsh API 和 `vtkXMLUnstructuredGridWriter` 分别写出 MSH 4.1 ASCII、单 Piece VTU ASCII 输入。它们是官方写入器序列化的人工指定网格，不是 CAD 自动划分质量的验证。输入不经过 MPMC writer，也不做 XML 清洗。4 个正常输入分别经 MPMC 导入、`make_mesh_exchange_document`、两个 canonical writer 输出，形成 **8 条链路和 8 份转换报告**；官方读取器检查输入和最终输出。

- Gmsh 输入使用稀疏乱序节点标签、混合类型的体单元块和重叠体区域组，只显式写两个有标签的外边界面。其余面（包括内部共享面）必须由导入器构建。VTK 输入反转点存储顺序、按 901/31 排列单元，带精确 UInt64 单元 ID、点标量/三分量字段和单元双分量字段。
- 同一 C++ target 新增测试专用 `--convert-3d <gmsh|vtu> <输入路径> <输出前缀>`。`convert_3d_file.cpp` 在 canonical 导出前将真实导入结果写为 `.import` 快照；Python 独立核对顶点、单元、全部 face↔cell 关联、边界分类/标签、显式面 ID、体积及面面积/owner 外法向。快照是待检证据，不提供期望值。最终文件再次检查几何、类型、节点身份、分组与字段实体绑定；合法实体重排不视为错误。
- 正常共享面的计数、3/6 m² 面积、16/32 m³ 总体积和 `1e-12` 容差沿用上节。检查只按节点身份建立邻接，不按相同坐标合并。VTU 节点 ID 由 UInt64 数组读回，稳定 face ID 由完整 FieldData 表读回并与实际面的节点关联核对，转换报告必须与身份绑定的实际变化一致。
- 同格式链路要求声明的网格/字段/分组语义保持；Gmsh→VTU 明确报告组、面标签损失，保留 face ID，VTU→Gmsh 报告字段损失。源 VTU 不提供单位属性或边界组，不能把导入生成的默认元数据冒充外部文件原有信息。`lossless` 不表示 XML 字节或可重算范围缓存原样保留。
- 另由官方写入器生成 **4 个接口断开输入**（同坐标不同节点身份），通过相同导入和两种导出，形成 **8 条断开对照链、8 份额外报告**。它们是合法网格，转换应成功；导入快照与最终文件必须均保留 0 内部面，体积不变，并明确拒绝“应有一个共享面”的断言。静默焊接或错误拒绝合法输入都会失败。

这条真实链路暴露并修复了共享 VTU ASCII 解析器对官方写入器所附 `vtkDataArray/L2_NORM_RANGE` 的拒绝。现在仅接受尾部一个结构完整、索引 0/1、两个有限值的范围缓存，再只解析原数值载荷；未知键（包括 `UNITS_LABEL`）、错误属性/索引、重复块、额外数字/嵌套及不完整 XML 仍拒绝。缓存不参与几何或字段计算，也不作为科学单位。二维 `mesh.core.vtu_ascii_roundtrip` 和 `mesh.core.vtu_ascii_invalid` 增加公共入口回归，原严格 ASCII/binary 拒绝不变；三维由上述官方文件直接覆盖。依据：[VTK XML InformationKey 结构](https://docs.vtk.org/en/v9.3.1/design_documents/IOXMLInformationFormat.html)、[固定版本官方 writer 源码](https://github.com/Kitware/VTK/blob/v9.7.1/IO/XML/vtkXMLWriter.cxx)。

本测试唯一 owner 仍为 external compatibility 工程，无新增 CTest 或 workflow。`verify_3d_chains.py` 在同目录显式复用 `verify_3d_readers.py` 的解析几何/独立读取及面判据、`verify_2d_readers.py` 的数值比较；三者及 producer、公共 importer/writer 依赖均由中央路由覆盖。原二维/三维导出验证和公共上游样例保留。Windows MSVC 与云端 Ubuntu GCC 使用同一依赖、输入域及容差；生产头变更还由既有核心/离散/PETSc 下游选测。

本地使用下一节同一 producer/隔离环境，将脚本替换为 `verify_3d_chains.py`，输出到新目录。PASS 为 `inputs=4 chains=8 reports=8 disconnected_chains=8`。`result.json` 分列正常和断开链，保存导入/最终拓扑、体积、共享面面积/法向、身份损失判定、拒绝原因，以及输入、快照、输出和报告的 SHA256、脚本/工具版本。云端由同一 external job 调用一次，保存 `mesh-3d-chain-readback` artifact。本项仍不覆盖非共形接口、高阶、binary、任意畸变或科学求解。

## 官方生成、MPMC 导入和导出的二维完整链路

`verify_2d_chains.py` 复用二维混合样例（逆时针三角形 + 顺时针非矩形凸四边形），以官方 Gmsh API 和 VTK writer 各写一个输入，实际导入 MPMC 后经 `make_mesh_exchange_document` 分别导出 MSH/VTU，由官方软件读取输入和输出。共 **2 个正常输入、4 条链、4 份报告**。官方文件直接用于导入，包括 VTK writer 附加的 `L2_NORM_RANGE`，不清洗 XML。Gmsh/VTK 在这里序列化人工指定的网格，不代表 CAD 自动划分验证。

- 正常网格应有 **6 条唯一边、1 条内部边、5 条边界边**；三角形/四边形面积分别为 **1 / 5 m²**，总面积 6 m²，面积质心分别为 `(3,2/3)` / `(19/15,14/15) m`。共享边长度为 `sqrt(5) m`，四边形朝向三角形的单位法向为 `(2,1)/sqrt(5)`，另一侧相反。Gmsh 官方 Jacobian 积分、VTK `vtkCellSizeFilter` 面积和独立三角扇面积/质心分别与解析值比较；节点身份、类型、循环方向与边关联精确检查，数值沿用 `1e-12` 绝对/相对容差。
- Gmsh 输入保留稀疏乱序节点 ID、两条显式外边界边 701/702、边界组和重叠二维区域组；其余边由 MPMC 构建。VTK 输入反转点顺序，单元按 901/31 排列，包含 UInt64 单元 ID、点标量/三分量字段及单元双分量字段。读回检查字段与实际节点/单元身份的绑定，不按坐标去重。
- 同一 C++ target 的 `convert_2d_file.cpp` 通过 `--convert-2d <gmsh|vtu> <输入路径> <输出前缀>` 调用实际 importer，在统一导出之前写出原始 `.import` 快照。检查全部 cell↔edge 双向关系、显式/生成边身份、面积和面积质心、边长/边中点、owner 外法向及 Gmsh 边界标签。二维 `VtuImportResult` 没有边界元数据，快照以 `-1` 标示缺失；边界计数由实际关联数推导，不能冒充输入已提供标签。
- 同格式报告须为 `lossless`；Gmsh→VTU 必须报告组、边界标签损失，保留 face ID，VTU→Gmsh 必须报告字段损失。VTU 节点身份直接读回，边身份通过实际读取的 FieldData 与稳定端点 ID 核验，并逐实体对照报告；不把相同 ID 集合当作绑定未变。
- 官方 writer 另生成 **2 个断开输入**，仅给三角形复制共享边两端节点，坐标不变、身份不同。形成 **4 条对照链及 4 份额外报告**。转换必须成功，导入和最终文件均保持 **7 条边、0 内部边、7 边界边**且面积不变；最终读回必须因共享边关联缺失而拒绝连接网格契约，防止静默焊接。

唯一执行 owner 仍为 external compatibility job，无新增 CTest/workflow。脚本显式复用同目录 `verify_2d_readers.py` 的二维解析判据，`verify_3d_readers.py` 的官方 VTU reader，以及 `verify_3d_chains.py` 中维度无关的官方 VTU 写入、关联/字段/报告检查；VTK 类型映射由二维调用方显式提供，三维默认行为保留。共享依赖和新增文件的直接命中、删除/重命名、生产头间接命中由中央路由回归覆盖。旧公共样例、二维/三维导出和三维完整链均保留。

本地将下节脚本替换为 `verify_2d_chains.py`，使用新的输出目录。总 PASS 为 `inputs=2 chains=4 reports=4 disconnected_chains=4`；`result.json` 保存 4 份导入快照对应的度量、8 份最终输出结果、控制拒绝原因、全部 24 个输入/快照/输出/报告文件的 SHA256，以及工具版本和 4 个 oracle 脚本 hash。中央/手动入口在原 job 中各按其触发模式执行一次，保存 `mesh-2d-chain-readback` artifact。范围仍为 XY 平面共形线性网格、MSH 4.1 ASCII 和单 Piece VTU ASCII，不覆盖曲面、高阶、binary、任意畸变、大规模性能或科学求解。

## 稳定节点身份的统一 VTU 契约

二维和三维都通过 `PointData/mpmc_global_vertex_id` 保存 Points 顺序对应的稳定节点 ID。写出固定为 UInt64、单分量；读入接受 UInt64 全范围（含 0）或非负 Int64。只有数组缺省时才采用旧 1..N 编号；数组存在时严格检查格式、类型、分量数、整数解析范围、值数、可选 NumberOfTuples 与唯一性。身份保留名不能成为普通科学字段，也不能放错 association。生产依赖不变，仍为标准库；具体名称是 MPMC 对标准 VTK 数组的约定，不声称第三方软件自动赋予它全局身份语义。

- 既有二维/三维所有导出文件都通过官方 VTK 的 `vtkUnsignedLongLongArray.GetValue` 逐点核验身份和顺序，不用 double 访问器。原负对照保留；复制接口节点时为其分配新的稳定 ID，确保仍然因共享面断开而失败。
- 既有完整链中的无节点 ID 官方 VTU 输入继续验证旧文件兼容。另在 `verify_2d_chains.py` 增加 2 个、`verify_3d_chains.py` 增加 4 个官方 VTU 输入：分别复用连通/断开混合几何，反转点存储，使用 `0`、`2^64-1`、`2^53+1` 及稀疏小整数。链路为官方 VTK 写入 → MPMC 原始导入快照 → 统一 VTU 导出 → 官方 VTK 读回，精确检查身份/坐标/连接/字段绑定、几何不变量与 lossless 报告。相同坐标不同身份不得焊接。转换工具附带生成的 Gmsh 文件不计入这 6 条 VTU 身份链的独立验收。
- `result.json` 在原 `cases` 与计数之外单列 `vertex_identity_cases`、`vertex_identity_inputs`；原完整链 PASS 计数仍指原 Gmsh/旧 VTU 链，新增身份链逐例输出独立 PASS。仅计上述旧链与节点链（含附带产物）为二维 36、三维 72 个文件；后续面身份链另加证据，全部保留文件和脚本 hash。
- `mesh.core.exchange_io` 同一 owner 增加两个公共 importer 的精确整型/旧格式/往返回归，每维 16 个非法节点数组场景和 10 个保留字段冲突；原 16 个身份损失场景保留并按新契约更新。节点和面身份始终保留；二维稀疏/置换节点 ID 及存储重排后，逐实体核对读回身份绑定，报告与实际一致。

本切片主要在本地执行核心、必要离散下游和所有独立读回脚本；保留现有 CI owner、依赖选测和专项入口，Draft 不启动云端测试。完整功能收尾再执行原多平台/sanitizer/PETSc 验收，不将本地 native 通过冒充这些证据。

## 稳定面身份的统一 VTU 契约

在 `UnstructuredGrid/FieldData` 写出三个单分量 UInt64 数组：`mpmc_global_face_id`、`mpmc_face_vertex_offsets` 和 `mpmc_face_vertex_ids`。offsets 为每条记录累计结束位置，不含起始 0；vertices 保存稳定节点 ID，二维宽度为 2、三维为 3 或 4。全表覆盖内部/边界面，记录顺序与节点顺序不构成身份；面 owner/外法向仍由拓扑和几何建立。该应用约定使用 [VTK XML FieldData](https://docs.vtk.org/en/latest/design_documents/IOXMLTimeInFieldData.html) 机制，不增加主单元、不声称 VTK 内建面身份或任意 filter 自动维护表。

- 生产二维/三维共用解析器/写入器，再通过各自既有构建器的显式面注释恢复 ID。旧文件没有整个表时继续生成面 ID；部分表、重复记录/ID、非法整型/偏移、未知节点、不存在的面、缺失面或错误 association 拒绝。UInt64 全范围（含 0）及非负 Int64 直接整型解析。
- 唯一核心 owner 仍是 `mpmc_mesh_exchange_io_tests` / `mesh.core.exchange_io`（mesh core 工程）。新增每维 52 个非法表场景、非负 Int64 读入、旧文件回退；原 16 个往返场景覆盖稀疏 ID、重排、0、UInt64 最大值及超过 2^53 的整数。身份字段冲突覆盖 point/cell 各 5 个保留名。
- 所有旧导出读回和完整链新增官方 VTK `GetFieldData()` 的精确整型读取，并以原生 `GetEdge()`/`GetFace()` 独立重建拓扑覆盖。Gmsh→VTU 实际保留原显式及生成面 ID；组/标签损失继续报告。原官方无表输入仍覆盖兼容。
- 原节点身份链继续保留；面身份另加二维 2 条、三维 4 条官方写入→MPMC 原始导入快照→统一 VTU 导出→官方读回链，涵盖连通/同坐标不同身份的断开混合网格。节点及面 ID 都含 0、2^64-1 和超过 2^53 的整数；点、单元、表记录顺序被打乱，面节点列表循环移位。检查 ID→节点集合→邻接单元、几何量、外法向及科学字段绑定，报告必须 lossless。
- 每条新链再加 4 个负对照：重复面 ID、未知节点、缺失面须被 MPMC 拒绝；交换两个面 ID 的合法文件可被官方 VTK 打开，但独立绑定判据必须拒绝。每个 `result.json` 单列 `face_identity_cases`/`face_identity_inputs`、精确绑定和负对照原因；原 PASS 计数仍指旧链，新链逐例输出 PASS。附带 Gmsh 文件不计入这些高位/零身份 VTU 链的独立验收。

共享 oracle 继续由现有 external compatibility target/job 拥有；没有新增 CTest、workflow 或生产依赖。几何判据、容差和原科学下游不变；这是格式与离散几何验证，不是物理模型验证或性能测量。遵循本地优先与完整功能云端验收规则。

## 隔离依赖与复现

选择官方读取器而非继续自读自写，能直接检查第三方软件是否理解输出；不引入额外通用格式转换层。依赖仅存在于测试 venv，生产 `mpmc::mesh` 仍只有标准库。固定完整 Python wheel 依赖闭包，禁止测试期间静默升级。Gmsh Python 包为官方 SDK（GPL-2.0-or-later），VTK 为 BSD-3-Clause；不复制其实现或把库链接进生产 target。VTK 的 matplotlib 依赖一并固定，但验证不使用绘图后端。

官方依据：[Gmsh API/格式文档](https://gmsh.info/doc/texinfo/gmsh.html)、[Gmsh SDK 包](https://pypi.org/project/gmsh/)、[VTK Python API](https://docs.vtk.org/en/latest/api/python.html)、[VTK 包](https://pypi.org/project/vtk/)。

在仓库根目录，Linux 示例（Python >= 3.12；`build_root` 和 `output_dir` 必须在仓库外，输出目录须尚不存在）：

```bash
build_root=/absolute/agent-workspace/mesh-external
output_dir=/absolute/local-artifacts/mesh-readback-run
sudo apt-get update
sudo apt-get install -y python3-venv libglu1-mesa
python3 -m venv "$build_root/venv"
"$build_root/venv/bin/python" -m pip install --only-binary=:all: -r tests/mesh/external_compatibility/requirements-readers.txt
cmake -S tests/mesh/external_compatibility -B "$build_root/build" -DCMAKE_BUILD_TYPE=Release
cmake --build "$build_root/build" --target mpmc_mesh_external_compatibility --parallel 2
"$build_root/venv/bin/python" -B tests/mesh/external_compatibility/verify_2d_readers.py --producer "$build_root/build/mpmc_mesh_external_compatibility" --output-dir "$output_dir"
```

Windows 使用同一 requirements 和脚本；venv Python 路径改为 `venv\Scripts\python.exe`，MSVC 构建加 `--config Release`，producer 路径为 `build\Release\mpmc_mesh_external_compatibility.exe`，不需要上述 apt 步骤。中文输出路径由 Python 处理，C++ producer 在该目录中以相对路径写文件。

成功要求返回码 0 且总 PASS 为 `files=12 reports=6 negative_controls=5`。`result.json` 记录真实读取器版本、平台、Python、producer/oracle SHA256、逐文件结果、输出文件 SHA256 及失败原因。CI 另记录提交 SHA 并保存结果与文件为 artifact。本地运行需同时记录 checkout 的提交/dirty 状态及构建命令；文件 hash 不替代源码版本。

PR 中央入口与原专项手动入口执行完全相同的步骤，仍由既有 external compatibility job 拥有，官方 Ubuntu runner、10 分钟上限和旧公共样例断言不变。失败返回码直接阻断 job，不能把成功生成文件算作读回成功。
