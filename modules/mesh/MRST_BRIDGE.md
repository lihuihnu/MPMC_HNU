# MRST / FaceMesh 交换契约 v1

普通 MRST 二维有限体积网格和三维显式面网格使用带版本的单文件 HDF5 交换，VTU 用于通用软件交换与显示。CSV 文件夹仅保留为历史迁移入口。HDF5 是容器，以下 schema 是 MPMC 约定，不是任意 `.mat` 文件或 MRST 官方文件标准。

## 接口和边界

| 入口 | 能力 | 限制 |
| --- | --- | --- |
| `FaceMesh` / `validate_face_mesh` | polygon/polyhedron、显式共享面、方向、64 位稳定 ID、typed arrays、NNC | 传输校验不认证几何或离散方法 |
| `read_face_mesh_hdf5` / `write_face_mesh_hdf5` | 全部 v1 数组保留；可选 HDF5 1.14.6 | 未知版本、额外条目、软/外部链接、隐式类型转换拒绝 |
| `read_face_mesh_vtu` / `write_face_mesh_vtu` | 单 Piece ASCII；polygon(7)、polyhedron(42)，并读 tri/quad/tet/hex/wedge/pyramid | 不支持 binary/appended/compressed、多 Piece、高阶；不四面体化 |
| `make_face_mesh` / `import_face_mesh_vtu` | 旧 Gmsh/GRDECL/VTU 文档迁入面契约，保留字段来源、标签、组、逻辑角点数组 | VTU 按显式元数据版本分派，不在解析失败后宽松回退 |
| `export_face_mesh_legacy` | 标准线性单元回到已有 Gmsh/GRDECL writer | 一般多面体/分裂面 unsupported；成功仍报告完整面 schema 等元数据遗漏为 lossy |
| MATLAB `export_mpmc_mesh` / `import_mpmc_mesh` | `G`、数值 rock、reference geometry、indexMap、稳定 ID、NNC 双向交换 | 未定义 struct 字段拒绝；粗网格/嵌入面/专用裂缝/LGR 未定义 |

旧三格式 API 保持原契约，一般 MRST 网格走 FaceMesh。新 exporter 拒绝覆盖文件。I/O 中断可能留下不完整文件，调用者应使用新路径并核验返回状态。

## HDF5 schema

所有 dataset 为 **rank 1**；矩阵按行展开。points 为 `x0,y0,z0,x1,y1,z1,...`，MATLAB 使用 `reshape(A',[],1)`，避开多维布局差异。整数 ID 不经过 double。

| 路径 | 类型 | 含义 |
| --- | --- | --- |
| schema_version / dimension | UInt64，各 1 项 | 版本 1；维数 2/3 |
| length_unit / z_convention / provenance | UInt8 UTF-8 字节 | 单位必须 m；Z 为 depth/elevation，不自动翻转 |
| points | Float64，3×节点数 | 2D 为 XY 平面，Z=0 |
| node_ids / face_ids / cell_ids | UInt64 | 各类内唯一，允许 0、超过 2^53、UINT64_MAX |
| face_offsets / face_nodes | UInt64 | 面→节点 CSR；零基局部引用；offsets 首 0、末等于连接数 |
| cell_offsets / cell_faces | UInt64 | 单元→面 CSR；同上 |
| signs | Int64 | 与 cell_faces 对齐，仅 ±1 |
| arrays/NAME/values | Float64/Int64/UInt64 | 行优先数值，有限实数 |
| arrays/NAME/components | UInt64，1 项 | 每行分量数大于零；允许零行 |
| arrays/NAME/location / unit | UInt8 UTF-8 | node/face/cell/incidence/nnc/metadata；单位不可空 |

数组名只允许 ASCII 字母、数字和下划线。arrays 组即使为空也必须存在，每个属性组恰有四个 dataset。C++ 检查类型、rank、CSR、ID 唯一性、范围和关联数量。默认 dataset payload 读取预算 8 GiB，可由 library 参数改动；这**不是整个进程 RSS 上限**，校验/解析还需临时空间。

3D 面环右手法向、2D 有向边右侧法向为正；sign=+1 表示该方向朝本单元外，-1 表示反向。共享面两侧符号相反，每侧最多一个单元；孤立面拒绝。MRST neighbors 第一列对应正侧、第二列负侧，0 为外部。`G.cells.faces(:,2)` 是标签，不是方向。

约定属性：

- `mrst_index_map`：cell/UInt64/1 分量，保留 MRST **一基**活动→逻辑映射，与 stable cell ID 分离。`cart_dims` 为 metadata/UInt64，检查正值、维数、乘积溢出和 indexMap 范围。
- `cell_face_tags`、`face_tags`、`mrst_*_global`：保留原值。global 若对应完整逻辑网格而非活动实体，关联 metadata，不截断或错绑。
- `reference_face_areas/normals/centroids`、`reference_cell_volumes/centroids`：保留原始参考值。五个名称均为保留名，读取/写入前按名称强制 `location/components/unit/type`：face areas=`face/1/(m|m2)/Float64`、face normals=`face/3/(m|m2)/Float64`、face centroids=`face/3/m/Float64`、cell volumes=`cell/1/(m2|m3)/Float64`、cell centroids=`cell/3/m/Float64`（括号内按 2D/3D）。任何 association 或元数据不一致均在解释参考值前拒绝。normals 为面积向量，重算不覆盖参考。
- `rock_*`：cell 属性；perm 为 m²、poro 为 1。其他未定义单位记 unspecified，物理计算前须明确。
- `nnc_cells`：nnc/UInt64/2 分量，零基单元端点。其他 nnc_* 与其行数一致。NNC 单独进入相邻图，不伪造几何面/法向。

MATLAB stable IDs 保持 UInt64；恢复 MRST 所需 double 索引时拒绝精度损失。重写时传回 import 返回的 metadata，保存来源和扩展属性；改变/重排实体后须同步字段。原 G.type 历史保存在 provenance，重建 G.type 标明 import_mpmc_mesh。不推导源文件缺失的 inactive 几何/岩性。

## VTU 面语义

标准 polygon/polyhedron 承载单元，mpmc_global_vertex_id / mpmc_global_cell_id 承载稳定 ID。dataset-level FieldData 使用：

- `mpmc_fm_version`=1、dimension、z、provenance（后面三者同 mpmc_fm_ 前缀，文本为 UInt64 字节值）。
- `mpmc_fm_face_ids/face_offsets/face_nodes/cell_offsets/cell_faces/signs`，对应 HDF5 数组。
- 每个属性由 `mpmc_fm_values_NAME`（原类型和分量数）、`mpmc_fm_location_NAME`、`mpmc_fm_unit_NAME` 定义。

node/cell 属性同时写 PointData/CellData；读入时可见值与 FieldData 必须精确相同，标准单元面环也必须与 sidecar 一致，不能借元数据掩盖损坏 connectivity。支持旧 faces/faceoffsets 以及 VTK 9.7 的 face_connectivity/face_offsets/polyhedron_to_faces/polyhedron_offsets。

外部新文件没有该契约时面 ID 新建；无损承诺从完整 schema 建立后开始。任意外部过滤器可能删除 FieldData，不能因为文件仍可显示就声称元数据无损。

## 几何和计算准备

`prepare_face_mesh_geometry` 独立于传输：检查有向边界闭合/连通；MRST 面间边分段不同但共线时，仅在验证中建立共同细分，不改原面环和 ID，非共线缺口拒绝。3D 以面节点算术均值构造共享三角扇，标量面积为三角片面积之和，面积向量为有向向量之和；非平面面两者不能混同。反向折叠扇形、零面积拒绝。

体积和质心由局部原点下有向四面体及一阶矩计算，保留负子四面体贡献，不 abs/裁剪。2D 使用有向边界面积积分。检查正体积及面积向量闭合。参考数组存在才比较，报告实际比较的标量数量；数量零不代表与 MRST 参考一致。

容差：u 为 double epsilon，X 为最大坐标绝对值（至少 1 m），L 为局部单元尺度。共线距离界为 `u*(64*X+4096*L)`。面积向量闭合界为 `4096*u+64*u*X*edge_length/surface_area`，二维相应使用边数。正体积阈值为 `4096*u*L^dimension`。参考比较默认 `1e-8*max(1,abs(actual),abs(reference))`，按相应 SI 量度量；面积/体积另加坐标量化前向误差项 `64*u*X*perimeter` / `64*u*X*surface_measure`。界由坐标精度和运算尺度确定，不按源残差调整。报告 reference_error 为该 max(1,...) 归一化误差，可能大于 1e-8 而落在量化误差项内，并非承诺纯相对误差始终小于 1e-8。

本版本适用域为上述三角扇表面定义；不证明非相邻单元无相交，不为任意自交网格/折叠面/所有凹面认证。非平面面的单一 centroid/area vector 不能充当任意线性通量的精确积分，离散层须选择相应积分规则。

`face_mesh_topology` 按需构建既有 Topology 四个关系；LocalIndex/CSR 仍为 32 位，超容量拒绝。`face_mesh_cell_graph` 构建对角、面邻接和 NNC 的对称标量图，供编号/稀疏结构准备。它不是通用 PDE 的最终 stencil，不含传导系数、TPFA 认证、未知量或数值装配；这些继续由 discretization/solver 拥有。

## 构建与使用

HDF5 为可选 target `mpmc::mesh_hdf5`；默认 `mpmc::mesh` 仍只依赖标准库，不自动联网、不要求 MATLAB。准备脚本固定 HDF5 1.14.6 提交 `7bf340440909d468dbb3cf41f0ea0d87f5050cea`，安装到仓库外。第三方源码不入库，HDF5 许可证/notices 随安装保留。

从仓库根目录运行 PowerShell：

```powershell
$agentDir = (Resolve-Path '../../03_agent_workspace').Path
python tests/mesh/core/prepare_hdf5.py --root "$agentDir/deps/mesh-hdf5"
cmake -S tests/mesh/core -B "$agentDir/b/face-mesh" `
  -DMPMC_MESH_WITH_HDF5=ON -DHDF5_USE_STATIC_LIBRARIES=ON `
  -DHDF5_ROOT="$agentDir/deps/mesh-hdf5/install"
cmake --build "$agentDir/b/face-mesh" --config Release --target mpmc_mesh_core_all --parallel 2
ctest --test-dir "$agentDir/b/face-mesh" -C Release --output-on-failure
$convert = "$agentDir/b/face-mesh/mesh-library/Release/mpmc_mesh_convert.exe"
& $convert input.h5 output.vtu --geometry --graph
& $convert output.vtu returned.h5 --geometry --graph
```

二维 VTU 导入加 `--2d`。GRDECL 导入要求 `--grdecl-si` 明确数值为 m/m²，其他单位使用原 library 缩放接口。HDF5 自带 dimension。Linux 同样配置，加 `-DCMAKE_BUILD_TYPE=Release`，单配置可执行文件路径没有 Release 子目录。

CLI 退出码：0 请求步骤通过；1 输入/输出/转换错误；2 传输已写出、几何失败。返回 2 时保留文件/诊断，仍可构建 graph，不可冒充可求解网格。HDF5 文件名使用 UTF-8，Windows 中文路径有回归覆盖。

MRST startup 后，从仓库根目录：

```matlab
addpath('modules/mesh/matlab');
G = computeGeometry(cartGrid([3 2 2], [3 2 2]));
rock = struct('poro', .2*ones(G.cells.num,1), 'perm', 1e-12*ones(G.cells.num,1));
export_mpmc_mesh('grid.h5', G, rock);
[Gback, rockback, metadata] = import_mpmc_mesh('returned.h5');
export_mpmc_mesh('mrst-return.h5', Gback, rockback, metadata);
```

实际文件使用仓库外绝对路径；示例文件名仅表示接口。

## 验证和依据

核心测试覆盖解析五边形/立方体、非平面面、平移稳定性、共线分段/NNC、错误方向/身份/CSR/参考值、旧格式适配、HDF5 超预算/空数组/不覆盖。独立 VTK 9.7.1 生成共享边 polygon/quad 和共享面五棱柱/hex，经 C++→HDF5→VTU→官方 VTK 读写→C++，检查精确数据、解析面积/体积、连接及 14 个负对照。MATLAB 本地 `verify_mrst_bridge.m` 覆盖 cartGrid 2D/3D、PEBI、大 ID、rock、空/非空 NNC 和错误字段。

`mrst_folder_audit.py` 仅迁移历史 numeric CSV，记录源 hash 和未迁移条目，不把 history/type 等任意原文件冒充无损迁移。大型地质模型留本地；云端执行代表性小网格的同样不变量、独立读回和原有多平台/sanitizer 检查，不依赖 MATLAB 许可证。命令见[外部测试](../../tests/mesh/external_compatibility/README.md)。

2026-10-04 本地 MRST 2024b / MATLAB R2022b：13 份历史网格全部 schema 数据集精确往返，11 份通过几何，2 份 Norne 在零基单元 118 边界不闭合被拒绝几何认证；源数据未改动。SPE11C 的 2,718,810 单元、8,215,725 面、2,778,528 节点完成参考值比较、往返和 18,911,652 项相邻图。DQ、SAIGUP、Statfjord 另经 MATLAB 读回和 computeGeometry 复算。这是特定文件的几何/传输证据，不是储层物理验证或任意 MRST 网格认证。

依据：[SINTEF MRST core/grid 文档](https://www.sintef.no/contentassets/2551f5f85547478590ceca14bc13ad51/core.html)（网页较旧，公式另核对本地 MRST 2024b computeGeometry.m）、[VTK XML 格式](https://docs.vtk.org/en/latest/vtk_file_formats/vtkxml_file_format.html)、[HDF5 1.14.6](https://github.com/HDFGroup/hdf5/releases/tag/hdf5_1.14.6)。调用独立软件，不复制其实现；地质数据不上传远程。
