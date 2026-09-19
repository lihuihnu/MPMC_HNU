# Mesh 模块：基础契约

`mpmc::mesh` 面向后续多相多组分流动离散，负责网格拓扑、几何、字段、求解自由度布局、文件 I/O 与并行分区元数据。网格层不得依赖 thermodynamics、flash、physics、runtime、前端或具体流动方程；PETSc/MPI 只允许出现在可选适配层，公共核心头文件不得泄漏 PETSc 类型。

> 当前状态：core topology/index、2D Cartesian topology/geometry、`FaceBoundarySnapshot`、`DenseFieldSnapshot`、`DofLayout`、`PartitionSnapshot`、`DofNumberingSnapshot` 与 `SharedEntityPlan` 已建立；Gmsh MSH 4.1 ASCII 与 VTU ASCII 已有 2D import/export round-trip，GRDECL 已有 raw corner-point parser 与无-fault active-cell shared-face processor。processor 现在进一步按 `source_logical_cell_ids` 投影 `PORO/PERMX/PERMY/PERMZ` 为 processed-topology `DenseFieldSnapshot`，并为每张 quad face 计算 centroid、area、owner 与 owner-relative unit normal。processed cell `GlobalEntityId` 和 source mapping 仍保留原 logical cell ID。`2×1×1` processed GRDECL shared topology 已实际送入 serial DMPlex：现有 topology adapter 现在严格接受 2D triangle/quad 或 3D hexa/quad 两类输入，3D 路径保持 `cell->face->vertex` 的部分插值表示并验证 hexa/quad/point strata、cone/support 与 stable cell/face/vertex identity；`vertex_coordinates_m` 也已写入 serial 3D DMPlex coordinate section/local Vec，只有 vertex points 拥有 3 个 coordinate DoFs。测试从 PETSc 持有的 coordinates 与 cone/support/closure 独立重算 11 张 quad face 的 centroid/area/owner-relative unit normal 及 2 个 hexa cell volume，并逐项与 `FaceGeometry3D`/`cell_volumes_m3` 对照。当前 3D topology depth=2、dimension=3，仍不伪造 edge；rooted `2×1×1` processed GRDECL 现在也已进入真实 2-rank `DMPlexDistribute(overlap=0)` 与 `DMPlexDistributeOverlap(depth=1)`：root rank 独占完整 25-point source DAG，non-root 初始为空，rooted DM dimension 从 serial source 广播而不再硬编码 2D；stable cell/face/vertex identity 通过 migration SF 重建并逐 stable ID 验证唯一 owner，3D coordinate section/local Vec 由 PETSc 随 distribute/overlap 自动迁移并按 stable vertex ID 对照原 `vertex_coordinates_m`。本 Gate 尚未迁移 `PORO/PERM*` 或 `FaceGeometry3D`。当前 GRDECL processing 仍只验证无 fault 的 `2×1×1` 与 `1×1×2`，fault split/pinch/NNC 继续显式拒绝或后置。core 仍不依赖 PETSc/MPI；不含 VTU binary/appended/compressed、Gmsh binary/high-order/3D、GRDECL fault/NNC processing、Mat、残差、求解器或流动物理。

## 1. 目标

本 PR 的最终目标是形成一个基本可用于科研流动计算的网格基础设施：

- 支持 1D/2D/3D 结构化与非结构化网格，并能表示常用 line/triangle/quad/tetrahedron/hexahedron/wedge/pyramid 单元；
- 使用稳定全局实体 ID、紧凑本地索引、连续数组/CSR 邻接表达 cell/face/edge/vertex 拓扑，避免每实体独立堆分配；
- 支持顶点坐标以及 cell/face 的质心、体积/面积、外法向、边界标记和必要邻接几何量的输入或一致计算，并显式拒绝退化/非法几何；
- 提供 location-aware 字段系统，可在 cell/face/edge/vertex 上存取标量、向量或张量数据；首批标准字段包含孔隙度、渗透率和传导率，同时允许后续物理模块注册自定义字段；
- 提供求解变量的拓扑 DoF 布局、局部/全局索引、owned/ghost 映射与连续存取接口，为离散残差、Jacobian 和线性/非线性求解器提供稳定索引；
- 支持常见科研网格文件的导入/导出，并保留可追溯的边界/物理组和字段数据；
- 提供并行 partition/overlap 元数据与可选 PETSc 适配，使内部拓扑、DoF 和 ghost 语义可直接映射到 DMPlex/PetscSection/PetscSF；
- 保持 C++20、语义清晰、模块边界明确、可独立测试和跨平台构建。

## 2. 核心数据契约

拓扑、几何、字段和求解布局必须分离。拓扑结构是网格身份；几何与属性可以更新，但不得通过隐式重排破坏稳定 global ID。公共索引至少区分：

- `GlobalEntityId`：跨分区稳定，用于文件、重分区和全局所有权；
- `LocalIndex`：进程内紧凑连续，用于热路径存取；
- entity kind：cell / face / edge / vertex；
- ownership：owned / ghost；
- partition owner rank 与 global-to-local/local-to-global 映射。

邻接使用连续 offset + index 数组（CSR 风格）或等价的紧凑结构，不允许把 `std::vector<std::vector<...>>`、链表或每实体多态对象作为生产热路径存储。常见邻接至少包括 cell->face、face->cell、face->vertex、cell->vertex；反向关系只能在确有消费方时物化，避免无条件重复存储。

`Topology` 是这一层的不可变 owning snapshot。`vertex/edge/face/cell` 的 global ID 分别连续保存；同一 `EntityKind` 内 ID 必须唯一，稳定实体身份定义为 `(EntityKind, GlobalEntityId)`，不同 kind 可保留各自来源编号。每个已物化 relation 只能占用一个 source-kind/target-kind 槽位，并且其 CSR 行数、目标计数必须与 snapshot 中对应实体数组严格一致。未物化的反向或派生 relation 不会自动生成。

2D Cartesian builder 使用零基、确定性编号：vertex/cell 按 `(j,i)` row-major；所有竖向 face 先编号，再按 `(j,i)` 编号横向 face。`cell->vertex` 固定为 logical lower-left/lower-right/upper-right/upper-left，`cell->face` 固定为 left/right/bottom/top；这些只是拓扑顺序，不代表已经计算几何法向。

## 3. 几何与字段

几何层至少提供：

- vertex coordinates；
- cell centroid 与 measure（1D 长度、2D 面积、3D 体积）；
- face centroid、measure 和相对 owner cell 的一致方向 normal；
- boundary/interface 标识；
- 非有限坐标、重复/越界 connectivity、零或负 measure、非法 orientation 的确定性错误。

当前 2D Cartesian geometry 已实现前三项的最小基线：x/y 轴坐标必须是有限、严格递增的 SI 米值；cell area 以 m²、face length 以 m 保存。face owner 取 canonical `face->cell` relation 的首个 cell，unit normal 从 owner cell centroid 指向 face centroid，因此边界 face 为 owner 的外法向，内部 face 则从 owner 指向另一侧。builder 会逐项核对四类 Cartesian relation，而不是仅凭实体数量假定 topology 兼容。

`CornerPointGeometry3D` 是当前 3D corner-point 的最小 owning snapshot：每个逻辑 cell 独占 8 个 corner vertices，顺序固定为 `{LLL,HLL,LHL,HHL,LLH,HLH,LHH,HHH}`（I 最快），坐标单位为 m，cell volume 为 m³。该类型刻意不把 faulted/split corner 强制合并为 shared vertices/faces；GRDECL parser 只物化 `cell->vertex`，以免在尚未做 fault/NNC processing 时伪造共享拓扑。volume 使用固定 5-tetra decomposition 检查 orientation：任何显著负 signed tetra 都视为 flipped/inverted 并拒绝；active cell 只要出现退化 tetra 或零总体积即拒绝，inactive cell 可保留退化 geometry 并记录 `volume=0`，但负体积翻转仍拒绝。这样 `ACTNUM=0` 不会导致 logical identity、corner coordinates 或属性被丢弃。

`process_active_corner_point_grid()` 是 raw corner geometry 到共享 active topology 的第一层处理，只针对无 fault 精确重合界面。processed cells 按 source logical I-fastest 顺序过滤 active cells，cell `GlobalEntityId` 继续使用 source logical ID，并在 `ActiveCornerPointGrid::source_logical_cell_ids` 中再显式保存 processed-local→source-global 映射。active corners 以 exact `(x,y,z)` 为 identity，所有 unique coordinates 先按字典序排序后生成 vertex IDs `1..N`；每个 hexa cell 固定生成 6 个 quad face slots，face identity 是 4 个 merged vertex local indices 的排序 key，所有 face keys 再按字典序生成 face IDs `1..N`。shared face 的 `face->cell` 支持数必须为 2，boundary 为 1；超过 2 个 cells 或非逻辑邻居共享整张 face 均拒绝。processor 还会先逐个 I/J/K 邻接检查 active-active logical neighbors 的四角坐标集合是否 exact 相同；不相同即明确报 unsupported fault/split geometry，不会悄悄变成两张 boundary faces。

processed field projection 不依赖 raw/processed local index 偶然相同：先用 raw cell `GlobalEntityId -> raw LocalIndex` 建映射，再对每个 `source_logical_cell_ids` 回查 `PORO/PERMX/PERMY/PERMZ`，按原 component layout 复制数值并原样保留 field metadata，最后用 processed `Topology` 重建 `DenseFieldSnapshot`。当前 raw GRDECL 四字段均为 cell scalar，因此 active filtering 后 entity count 与 processed cell count 严格一致。

`FaceGeometry3D` 保存与 processed face ordering 对齐的 centroid[m]、area[m²]、owner cell 与 owner-relative unit normal。每张 quad 使用 canonical cyclic 4-vertex order，并以对角线 `(v0,v2)` 分成三角形 `(0,1,2)` 与 `(0,2,3)`；face area 是两个 triangle area 之和，centroid 是按 triangle area 加权的 triangle-centroid 平均，normal 由两个 triangle area-vector 之和归一化。零面积 triangle、near-zero resultant area-vector 或两个 triangle normals 显著反向的 folded quad 都拒绝。owner 取 canonical `face->cell` 首个 cell；normal 最后用 `face centroid - owner cell 8-corner mean` 定向，保证点向 owner 外侧。当前仍未生成 3D edge entities，也未处理 fault split/pinch/NNC。

`FaceBoundarySnapshot` 与 geometry 独立，只消费 `Topology::face->cell`：一个相邻 cell 定义为 boundary，两个定义为 interior，0 个或多于 2 个都拒绝。每个 face 对齐保存 1-byte `FaceClassification` 与 32-bit `PhysicalTag`；tag `0` 保留为 untagged，非零 tag 只允许出现在 boundary face，同一 tag 可重复用于一个 physical group。这里不解释 tag 的任何压力/流量/壁面/井/材料语义。

字段系统必须记录 entity location、component count、数值类型语义与单位/来源元数据，不允许仅靠字符串猜测布局。

当前通用 `DenseFieldSnapshot` 只冻结存储契约，不内置任何具体物理字段：location 目前支持 vertex/face/cell，数值类型固定为 `double`，component count 必须大于零；底层采用 entity-major / component-interleaved 连续布局 `values[entity * component_count + component]`。构建时使用 `Topology` 对齐实体数量，拒绝长度不匹配、size overflow、NaN/Inf 与未声明单位。metadata 保存稳定 field ID、原样 unit 字符串，以及 source kind/reference/revision/locator；本层不猜测单位、不做转换，也不把 source metadata 当成真实性证明。缺失值尚无契约，因此不能用 NaN 代替。

后续 porosity、permeability、conductivity 等科研属性应建立在该通用容器之上，再分别定义 location、component layout、单位和物理有效域；本增量尚未添加这些规则。I/O 遇到未知单位、缺少必要分量或数组长度不匹配时仍不得静默补值或重排。

## 4. 网格来源与文件 I/O

本 PR 合并前至少覆盖以下实用入口：

1. 内建 Cartesian 生成：当前已有 2D `nx × ny` topology builder，以及由任意严格递增 x/y 轴坐标生成非均匀 Cartesian metric geometry；后续仍需 1D/3D；
2. Gmsh MSH：当前已建立 4.1 ASCII 2D import/export 基线，覆盖 `$MeshFormat/$PhysicalNames/$Entities/$Nodes/$Elements`，支持 2-node line、3-node triangle、4-node quad。import 允许 sparse/out-of-order node/element tags，并按 stable tag 确定性压紧 local ordering；调用方显式提供 `coordinate_scale_to_m`，仅接受可映射到 `Geometry2D` 的 XY-plane mesh。vertex `GlobalEntityId` 使用 Gmsh node tag，cell 使用 2D element tag；显式 line element 对应 face 使用其 element tag，缺失的内部 cell edge 由 importer 重建并从 `max(all element tags)+1` 起确定性生成 face ID。export 以 SI metre 写出 coordinates，为每个 vertex/face/cell 构造 point/curve/surface entity，并把 face 与 cell `GlobalEntityId` 原样写回 line/2D element tag，因此 generated internal face ID 也能 round-trip 保留。boundary `PhysicalTag`、surface multi-Physical-Group membership 与 `PhysicalNames` 均写回；curve entity 仍最多允许一个 Physical Group，因为 `FaceBoundarySnapshot` 每 face 只有一个 `PhysicalTag`。parametric nodes、binary、high-order、concave quad、非 XY-plane、3D、post-processing field sections 与通用 Gmsh data export 仍显式不支持；
3. VTK XML UnstructuredGrid：当前已建立 `.vtu` ASCII 单-Piece import/export round-trip。import 要求 `VTKFile type="UnstructuredGrid"`，`Points` 为 Float32/Float64 三分量且 z≈0，`Cells` 为标准 `connectivity + offsets + types`；connectivity/offsets 接受 Int32/Int64，types 接受 UInt8，cell type 只支持 VTK_TRIANGLE=5 与 VTK_QUAD=9。point/cell DataArray 当前只映射 Float32/Float64 为 `DenseFieldSnapshot`，component count 原样保留。稳定 cell identity 使用保留的 `CellData/DataArray Name="mpmc_global_cell_id"`（UInt64/Int64）；外部文件若缺失则按 cell 文件顺序生成确定性 `1..N`，export 总是写回 UInt64 identity array。vertex ID 当前按 point order 生成 `1..N`，faces 从 cell edges 确定性重建。字段数值走标准 PointData/CellData；为满足 core 强制的 unit/source provenance，export 额外写 `mpmc_unit/mpmc_source_kind/mpmc_source_reference/mpmc_source_revision/mpmc_source_locator` XML attributes，import 在存在时恢复，缺失时生成明确的 file-import 默认 metadata。export 使用 Float64 points/fields、Int64 connectivity/offsets 与 UInt8 types。binary、appended、compressor、多 Piece、3D/high-order、非三角/四边形以及整数型物理字段均显式不支持；
4. reservoir corner-point：当前已建立最小 GRDECL parser，要求 `SPECGRID/COORD/ZCORN/ACTNUM/PORO/PERMX/PERMY/PERMZ` 各出现一次并以 `/` 终止，支持 `N*value` repeat、Fortran `D` exponent 与 `--` 行注释；未知关键词（包括 `INCLUDE`）显式拒绝。`SPECGRID` 只接受 `NUMRES=1`、Cartesian `F`；COORD 必须是 `6*(NX+1)*(NY+1)`，ZCORN 是 `8*NX*NY*NZ`，ACTNUM/属性均为 `NX*NY*NZ`，logical cell index 固定 `i + NX*(j + NY*k)`。COORD pillars 按 `(i,j)`、I-fastest 顺序解释，ZCORN 使用 `{LLL,HLL,LHL,HHL,LLH,HLH,LHH,HHH}` corner ordering。调用方必须显式提供 `coordinate_scale_to_m` 与 `permeability_scale_to_m2`；parser 不猜 FIELD/METRIC。所有 logical cells 都进入 `Topology`，cell `GlobalEntityId=logical_index+1`，每 cell 的 8 个 corner vertex IDs 也确定性保留；`ACTNUM` 仅作为 activity mask，不压缩 cell ordering。`PORO` 强制 `[0,1]`、unit=`1`，PERM 值必须非负并按显式 scale 转成 `m2`。当前不构造 shared faces、fault intersections、pinch/NNC、MAPAXES、GRIDUNIT、LGR 或 export。 在其上当前已有最小 active-cell shared-face processing：无 fault active cells 可生成共享 vertex/face `Topology`；inactive logical cells 不进入 processed computational topology，但其原 logical IDs 仍可通过 processed mapping 追溯。

若某格式只实现声明的子集，必须在解析器入口与文档中显式拒绝未支持特性；不得“读取成功”后丢失高阶节点、物理组、inactive cell 或字段。

## 5. 求解变量与拓扑索引

最小 `DofLayout` 已实现 local scalar indexing 基线。每个 `DofVariable` 只包含稳定 ID、location 与 component count；location 当前支持 cell/face/vertex，component count 必须大于零，变量 ID 必须唯一。布局固定为 `[cell block][face block][vertex block]`；每个 location 内按 `entity-major -> variable declaration order -> component` 排列。同一实体上的变量 DoF 因此连续，且 location block 也连续。

布局提供每个 location 的 entity count、DoFs-per-entity、block offset、block scalar count，以及 `scalar_offset(variable_index, entity, component)` 的 O(1) 热路径映射；字符串 ID 到 variable index 的查询只用于控制路径。该顺序刻意接近后续 `PetscSection` 的 point/field 组织，但当前没有 PETSc 类型或依赖。

在其上，`DofNumberingSnapshot` 已建立 partition-aware scalar numbering。全局顺序与 `DofLayout` 保持同构：`[cell block][face block][vertex block]`，每类实体按 contiguous global entity ordinal，再按该实体上的 variable/component 顺序编号。strict serial 下 global entity ordinal 就是 local entity index，因此每个 global scalar index 与现有 local scalar offset 完全一致且连续。generic local 下，`PartitionSnapshot` 只知道本 rank 可见实体及 owner rank，无法单独决定跨 rank 连续 ordinal；因此 builder 必须额外接收 `(GlobalEntityId, GlobalEntityOrdinal)` 记录和每类 global entity count，并逐项与 partition 的 stable ID/local ordering 校验。不同 rank 对同一 GlobalEntityId 使用同一 ordinal 的全局一致性，必须由未来 distributor/MPI 层保证，当前 local snapshot 不冒充能够验证远端数据。

`DofNumberingSnapshot` 不为每个 DoF 复制 global/owner metadata，而按 entity 保存 ordinal 与 owner rank，再利用固定 DoFs-per-entity 算术映射。local→global 为 O(1)；global→local 使用按 ordinal 排序的 `LocalIndex` permutation 二分查找。它提供 local/global DoF 总数、owned/ghost local DoF 数和逐 local scalar ownership。尚未实现 constraint DoF、跨 rank communication plan 或 PETSc section/SF；该层仍不知道 pressure、saturation、composition 等具体物理意义。

## 6. 并行与 PETSc 边界

核心 `mpmc::mesh` 保持 MPI/PETSc 可选。当前 `PartitionSnapshot` 已冻结 serial/local ownership 基线：snapshot 记录 `local_rank` 与 `rank_count`，每类 entity 复制 `Topology` 中的稳定 `GlobalEntityId` 并按同一 local ordering 保存 `owner_rank`。ownership 不重复存储，而由 `owner_rank == local_rank` 推导为 owned，否则为 ghost，从而不存在 owner/ownership 两份状态漂移。local→global 为 O(1) 连续数组访问；global→local 保存一份按 global ID 排序的 `LocalIndex` permutation 并二分查找，避免每实体树节点或哈希桶。global ID 只要求在同一 entity kind 内唯一，所以查询始终携带 `EntityKind`。

严格 serial builder 固定 `rank_count=1`、`local_rank=0`，所有本地 vertex/edge/face/cell 都由 rank 0 拥有且 ghost count 必须为零。generic local snapshot 允许 owned/ghost 在 local ordering 中交错，不要求 owned-first；它只描述“本 rank 当前可见的 local topology”，不声称掌握全局所有实体。`DofNumberingSnapshot` 也遵循这个边界：它可以消费外部提供的 global entity ordinal 来形成连续 global scalar IDs，但不会从 owner rank 猜测缺失的远端实体顺序。

`SharedEntityPlan` 在这一 ownership 基础上冻结非 MPI halo exchange 元数据。输入 canonical `SharedEntityLink` 同时描述 owner 与 ghost 两端：`EntityKind`、稳定 `GlobalEntityId`、`owner_rank/owner_local`、`ghost_rank/ghost_local`。对本 rank 而言，若本地是 owner，则生成 send entry；若本地是 ghost，则生成 receive entry，并校验本地 stable ID、owned/ghost 状态以及 `PartitionSnapshot::owner_rank`。每个 local ghost 必须恰好出现一次；同一邻居上同一 `(EntityKind, GlobalEntityId)` 不得重复。owned entity 可以出现在多个不同 neighbor 的 send list 中，以支持一对多 halo mirror。

存储采用 `NeighborExchangeRange[] + flat send[] + flat receive[]`：neighbor range 按 rank 排序，send/receive 条目按 entity kind、global ID、local index 确定性排序。对 receive range，range.rank 就是 owner rank，entry.remote_local 就是 owner-local index；对 send range，range.rank 是 ghost rank，entry.remote_local 是远端 ghost-local index。这样 owner rank 不需要在每个条目里重复存储。当前仍未实现 halo depth、真实 send/recv buffer、MPI communicator、collective symmetry check 或跨 rank 自动构造；synthetic 2-rank fixture 只验证两份本地 plan 在相同 canonical links 上的 owner/ghost 对称性。

可选 `mpmc::mesh_petsc` 已建立第一条真实适配基线，但仍与 core 分离。当前 `create_section_mapping()` 使用连续 chart `[cell points][face points][vertex points]`，将 `DofVariable` 声明顺序映射为 `PetscSection` fields，并显式设置 point-major，使 PETSc point offset 与现有 `DofLayout::entity_offset()` 一致；`DofNumberingSnapshot` 同时转换为 PETSc `PetscInt` 宽度的 local-to-global scalar map，若全局编号超出当前 PETSc index 宽度则拒绝，而不是截断。

`create_entity_sf()` 仍可按单一 `EntityKind` 建立 `PetscSF`，其 remote root index 直接等于 core `owner_local`。为服务完整 local `PetscSection`，新增 `create_point_sf()`：adapter 在真实 MPI communicator 上 `Allgather` 各 rank 的 cell/face/vertex local counts，据此把 `SharedEntityPlan` 的 kind-local `owner_local` 转换为远端 `[cell][face][vertex]` flattened point index，而不向 core contract 塞入 PETSc 专用 point base。adapter 会核对 communicator 的 rank/size、layout/partition entity counts，以及 leaf 的 ghost/owner/GlobalEntityId 后再创建 graph。

2-rank gate 现在进一步调用 `PetscSectionCreateGlobalSection(localSection, pointSF, ...)` 与 `PetscSFSetGraphSection(sectionSF, localSection, globalSection)`。PETSc global section 的 owned offset 使用 PETSc 自身的并行 ownership layout；ghost point则保存负编码 `-(owner_offset+1)`，所以它的数值顺序并不强制等于 core 的 `[cell][face][vertex] + GlobalEntityOrdinal` 编号。回归先由 owned PETSc offsets 建立 `PETSc-global-offset -> DofNumberingSnapshot::GlobalDofIndex` 映射，再验证 owned 正 offset 与 ghost 负 offset 都解析到本地 `local_to_global` 指向的同一 core GlobalDofIndex。section-SF 对完整 10-DoF local array 执行 global-layout→local-layout broadcast，包含 cell point 上的 2-component `cell.primary` 与额外 cell DoF。

在此基础上，adapter 新增 `create_section_vecs()`、`global_to_local()` 与 `local_to_global_add()`。global Vec 使用 communicator 上的 `VECMPI`，其每 rank local size 严格等于 global section owned storage；local Vec 使用 `PETSC_COMM_SELF` 的 `VECSEQ`，size 等于 local section 全 storage。global→local 直接以 section-SF + `MPIU_SCALAR/MPI_REPLACE` 将 owner storage 广播到 owner/ghost local slots。local→global ADD 不把已有 global 值覆盖掉：先把 local `PetscScalar` contributions 通过 section-SF + `MPIU_SUM` reduce 到临时 global Vec，再用 `VecAXPY` 加回目标 global Vec。synthetic 2-rank fixture 显式确认每个 core global DoF 恰有两份 local copy，并用不同 rank 的 contribution 编码验证 owner 与 ghost 各参与一次、没有重复计数或 ownership 错位。

新增 `create_serial_dmplex_topology()` 作为 topology-only DMPlex 基线。它不让 PETSc 从 cell list 自动插值并重新生成 face，而是直接用 core `cell->face` 作为 cell cone、`face->vertex` 作为 face cone，再调用 `DMPlexSymmetrize()` 生成 support、`DMPlexStratify()` 建立 strata、`DMPlexComputeCellTypes()` 推导 quad/segment/point 类型。当前 point numbering 明确固定为 `[cells][faces][vertices]`，并返回独立 `DMPlexPointIdentity[]`：stable `GlobalEntityId` 保留为 core 的 64-bit 类型，不压入 `DMLabel/PetscInt`。构建前会拒绝 edge entity、缺失四类必要 relation、非 4-face/4-vertex cell、非 2-vertex face，以及 `cell->face` / `face->cell` 或 cell vertex closure 不一致。

synthetic serial 回归使用现有 `2×1` Cartesian quad connectivity，但把 vertex/face/cell stable IDs 分别提升到 5/6/7×10^9 量级，严格核对 DMPlex chart/height/depth strata、quadrilateral/segment/point cell type、cell cone、face cone、face support、vertex support，以及每个 PETSc point 的 64-bit core identity。该 serial builder 本身仍只在 `PETSC_COMM_SELF` 上工作，不注入伪造坐标或 metric geometry。

其上新增 rooted distribution gate。`create_root_dmplex_topology()` 在 rank 0 复用已验证的 serial adapter，再把同一 DAG 放到 `PETSC_COMM_WORLD`；其他 rank 初始 chart 为空。测试把 partitioner 固定为 `PETSCPARTITIONERSIMPLE`，调用 `DMPlexDistribute(..., overlap=0)` 后要求两个 rank 各拥有一个 cell，并从 distributed DM 的 point SF 重建每个 local cell/face/vertex 的 owner rank。stable identity 不存入 `DMLabel/PetscInt`，而由 `migrate_dmplex_identities()` 给 source points 建立每点两列 `uint64` 数据（kind + GlobalEntityId），通过 `DMPlexDistributeData()` 沿 migration SF 迁移到目标 DM，再按目标 cell/face/vertex strata 重新生成 kind-local `LocalIndex`。

对 overlap=0，回归对 2 个 cell、7 个 face、6 个 vertex 的每一个 stable ID 跨 rank 统计 owner 数，要求严格为 1；共享 closure points 的 point-SF leaves 必须能逐项解释为 `PartitionSnapshot` ghost。随后调用 `DMPlexDistributeOverlap(..., 1)`，再次沿 overlap migration SF 迁移 stable identity。对本 2×1 fixture，每个 rank 必须看见两个 cell，其中一个 owned、一个 ghost。测试从 overlap DM point SF 的 remote `(rank, point)`、各 rank strata ranges 与迁移后的 GlobalEntityId 构造 canonical `SharedEntityLink`，经 collective 汇总后交给现有 `SharedEntityPlan::create()`；最终要求 plan 的 receive 数等于 core ghost 数、两 rank send/receive 对称，且每个 point-SF leaf 的 owner rank 与 `PartitionSnapshot` 完全一致。PETSc point SF 的 `nroots` 按 DMPlex point-index space 的上界解释，而不是简单 chart size，这一差异已由真实 3.19.6 runner 定点修正。

distributed geometry gate 在同一个 rooted/distribute/overlap 链上加入真实 PETSc coordinates。`attach_root_geometry2d_coordinates()` 用 `DMSetCoordinateDim(2)`、`DMSetCoordinateSection()` 与 `DMSetCoordinatesLocal()` 把 `Geometry2D::vertex_coordinates_m` 写到 source DMPlex；只有 vertex points 拥有 2 个 coordinate DoFs，cell/face coordinate DoF 为 0。PETSc 的 DMPlex distribution 会随 migration SF 迁移 coordinate section/local Vec，因此 distributed 与 overlap DM 都直接通过 `DMGetCoordinateSection()` / `DMGetCoordinatesLocal()` 读取迁移后的坐标，而不是测试侧重新复制。

serial processed-GRDECL 3D geometry gate 使用独立 `attach_serial_vertex_coordinates_3d()`：要求 serial、dimension=3、vertex identity 与 coordinate count 一一对应且坐标有限，随后以 `DMSetCoordinateDim(3)`、3-DoF vertex coordinate section 和 local Vec 写入 SI metre 坐标；cell/face coordinate DoF 保持 0。由于当前 3D DMPlex 刻意保持 `cell->face->vertex` partial interpolation（depth=2）而没有 edge，gate 不调用要求 fully interpolated mesh 的 `DMPlexComputeCellGeometryFVM()`。测试只从 PETSc coordinate section/local Vec 与 DMPlex cone/support/closure 重建 quad 三角剖分 metric 和闭合表面体积，再把 face centroid/area/owner-relative normal 与 cell volume 逐项对照 core snapshot；因此没有通过再次调用 core metric helper 来形成循环验证。

rooted processed-GRDECL 3D distribute/overlap gate 继续复用同一 partial-interpolated topology，但 `create_root_dmplex_topology()` 现在从 root serial DMPlex 读取并广播真实 DM dimension，因此 2D/3D rooted paths 共用同一入口。`attach_root_vertex_coordinates_3d()` 只允许 root rank 提供完整 processed vertex coordinates/identities，non-root 必须为空；随后在 source DM 上建立 3-DoF vertex coordinate section/local Vec。官方 PETSc 3.19.6 2-rank gate 使用 `PETSCPARTITIONERSIMPLE` 将两个相邻 hexa 分到不同 rank，overlap=0 时每 rank 恰好一个 owned cell，depth-1 overlap 后每 rank 可见两个 cells（1 owned + 1 ghost）。identity 仍通过现有 migration SF 迁移；cell/face/vertex 每个 stable ID 全局恰有一个 owner。distributed 与 overlap DM 的坐标直接从 PETSc 自身 coordinate section/local Vec 读取，并按 stable vertex ID 与 root processed `vertex_coordinates_m` 对照；本 Gate 刻意不调用 boundary/property/face-geometry migration。

geometry 回归使用非均匀轴 `x={0,1.25,3.75}`、`y={-2,2}`，避免单位网格掩盖错误。每个 distributed/overlap vertex 先通过 stable `GlobalEntityId` 对照原 `Geometry2D` 坐标；每个 face 再从 DMPlex face cone 的两端 vertex coordinates 独立计算 Euclidean length；每个 cell 从 transitive closure 提取四个唯一 vertex，以坐标均值重算 centroid，并按围绕 centroid 的角排序后用 shoelace 公式重算 area。结果分别对照 `Geometry2D::face_length_m()`、`cell_centroid_m()` 和 `cell_area_m2()`，所以验证的是 PETSc 分发后的 metric geometry，而不只是 topology/identity。

distributed boundary/property gate 继续复用两段真实 migration SF。`migrate_face_boundary_snapshot()` 为每个 source face 建立 2×`uint32` payload：`FaceClassification + PhysicalTag`，经 `DMPlexDistributeData(..., MPI_UINT32_T, ...)` 迁移后按目标 face kind-local ordering 重建新的 `FaceBoundarySnapshot`。classification 也随 source 一起迁移，而不是在 overlap=0 的 local support cardinality 上重新猜；因此跨 rank 的内部 interface 即使某 rank 本地只看到一个 adjacent cell，也不会被错误改写成 physical boundary。

`migrate_dense_field_snapshot()` 以 source `DenseFieldSnapshot` 的 location/component_count 为 point-section DoF，使用 `MPI_DOUBLE` 迁移 entity-major components，并按迁移后的 `DMPlexPointIdentity` 重建实际 core `DenseFieldSnapshot`。gate 同时覆盖 cell scalar（1 component）、face field（2 components）和 vertex field（3 components）。字段 metadata 保持字段级 provenance，不作为 per-point payload 重复发送：`id`、`unit`、`FieldSourceKind::synthetic_test`、`reference`、`revision`、`locator` 在重建 snapshot 中逐项保持不变。测试执行严格两段链路：root snapshot → distribute snapshot → overlap snapshot；第二段只消费第一段重建结果。每个 target value/tag 最终都按 stable `GlobalEntityId` 回查 root reference，因此 owned 与 ghost copy 必须完全一致。

Gmsh importer 下游 gate 使用一个同一 surface 内相邻的 triangle+quad fixture：5 个 sparse/out-of-order node tags、5 个 physical boundary line elements、2 个 surface cells，内部共享 edge 不在 `$Elements` 中，因此 importer 必须生成 stable face ID 203。core 回归核对 topology relation widths、generated internal face、PhysicalTag、surface Physical Group、显式 SI scale 后的 centroid/area/face length，并覆盖 binary/high-order/missing-node/multi-boundary-group/nonplanar 拒绝路径。round-trip gate 进一步把每个 cell 的 surface membership 扩为 `{21,22}`，执行 `import(scale-to-m=2) -> export(SI metre) -> import(scale-to-m=1)`，逐项比较 vertex/face/cell stable IDs、四类 CSR relations、全部 `Geometry2D` metric、`FaceBoundarySnapshot`、`PhysicalNames` 与逐-cell surface Physical Groups。export 会为所有 faces（包括内部 face）写出 line element，因此 generated face 203 不会再次分配新 ID。随后同一导入结果在 2-rank PETSc gate 中由 rank 0 进入 rooted DMPlex，`PETSCPARTITIONERSIMPLE` 分发后仍保持一个 triangle 与一个 quad；coordinates、FaceBoundarySnapshot 与一个带 provenance metadata 的 cell `DenseFieldSnapshot` 分别经过 distribute 与 depth-1 overlap migration，并按 stable `GlobalEntityId` 在 owned/ghost 两侧重新核验。

VTU gate 使用同样的相邻 triangle+quad 2D 形状，`Points` 5 个点、`Cells` offsets=`{3,7}`、types=`{5,9}`，并包含一个 2-component point field 与一个 2-component cell field。cell identity 专门使用 `9223372036854775813`（超过 signed int64 上限）与 `7000000003`，验证 UInt64 不被缩窄；point field provenance 中含 XML `&amp;` 转义，验证 metadata escape/unescape。round-trip 执行 `VTU ASCII -> core -> VTU ASCII -> core`，逐项比较 vertex/face/cell IDs、四类 CSR relations、`Geometry2D`、point/cell field layout/value/metadata。invalid gate 覆盖 appended、compressor、unsupported cell type、non-planar coordinate、malformed offsets、duplicate cell IDs 与 reserved field-name collision。

GRDECL gate 使用 `2×1×1` axis-aligned corner-point fixture，第二个 logical cell `ACTNUM=0`。fixture 用 `8*0 8*1` 验证 repeat expansion、`PERMZ` 用 `D` exponent 验证数值解析，并显式把 coordinate scale 设为 2、permeability scale 设为 `1e-15 m²/source-unit`；两 cells 最终 volume 均为 8 m³，inactive cell 的 stable ID=2 与 PORO/PERM 值仍完整保留。另有独立 `1×1×1` degeneracy gate：active zero-thickness cell 必须拒绝，inactive zero-thickness cell 保留 identity/corners/property 且 volume=0，inactive flipped top/bottom 仍拒绝。invalid gate 还覆盖 radial `T`、非法 ACTNUM、缺 keyword、PORO>1、negative permeability、unknown INCLUDE 与非法 unit scale。

active face processor gate 使用两种正交共享方向：`2×1×1` 两 active cells 验证 I-interface，`1×1×2` 两 active cells 验证 K-interface；两者都必须从 16 个 cell-local corners 压成 12 unique vertices、从 12 张 cell-local faces 压成 11 faces，并且只有 1 张 `face->cell` support=2 的 shared face。相同 raw input 重复处理必须得到完全相同的 Topology 与 vertex ordering。所有单位立方体 face 都要求 area=1，normal 必须为单位向量且与 `face centroid-owner corner mean` 点积为正；I-interface 明确核验 centroid=`(1,0.5,0.5)`、normal=`(+1,0,0)`，K-interface 核验 centroid=`(0.5,0.5,1)`、normal=`(0,0,+1)`。双-active fixture 同时核对 PORO/PERMX 投影值与 metadata；activity mapping gate 将第一个 logical cell 设 inactive、只保留 logical cell 2，要求 processed local cell 0 的 topology cell ID 与 source mapping 都仍为 2，并核对 PORO/PERMX/PERMY/PERMZ 均来自 raw logical cell 2。另有 split-face fixture：两个 active logical neighbors 的共享 pillar face z 值故意错开，raw GRDECL geometry 仍可解析，但 processor 必须拒绝 unsupported fault split；全 inactive input 也拒绝，因为没有 computational cells 可输出。

当前 PETSc gate 固定在官方 `ubuntu-24.04` runner 的 PETSc 3.19.6。尚未实现通用 partition policy、超过 depth-1 的 overlap、Gmsh binary/high-order/3D、Gmsh NodeData/ElementData/ElementNodeData export、constraint DoF、真实 halo buffer abstraction、Mat integration 或残差/Jacobian；这些不得由当前 adapter/importer/exporter 冒充完成。

后续适配层仍可负责：

- 在已通过的 rooted processed GRDECL 3D distribute/overlap identity + coordinate gate 上，下一步应单独迁移 `PORO/PERMX/PERMY/PERMZ` processed cell fields，并按 stable cell `GlobalEntityId` 验证 owned/ghost 值与 metadata；`FaceGeometry3D` 继续作为下一条独立数据链，不与 property migration 合并；
- 在已有 point/global/section SF 与 Vec 基线上加入 constraints 与稳定 Mat integration；
- 使用 PETSc 的分发/overlap 机制验证 partition 与 ghost；
- 保持 PETSc 对象生命周期和错误码不穿透到核心网格接口。

PETSc 可使用其已支持的 partitioner；核心模块不复制 ParMETIS/PT-Scotch 算法。没有 PETSc 时，仍必须能够构建、导入、检查、索引和以显式 partition plan 构造本地/ghost 视图。

## 7. 性能与内存原则

- 实体 ID、connectivity、几何和数值字段使用连续存储；
- indexed read 不分配内存、不做字符串查找、不进行隐式拓扑构造；
- 解析阶段允许临时工作区，finalize 后释放不再需要的重复结构；
- 64-bit global ID 与紧凑 local index 分离，避免所有本地邻接无条件使用 64-bit；
- 性能结论必须来自固定 workload 的 benchmark，不以单次 CI wall time 宣称加速；
- 首个稳定实现建立内存/访问基线，后续无依据不得显著回退。

## 8. 本 PR 合并条件

本 PR 只有在以下条件同时满足时才能从 Draft 转为可合并：

- 核心拓扑/几何/字段/DoF/partition API 已实现且文档与代码一致；
- Cartesian、Gmsh、VTK/VTU 和声明的 GRDECL 子集达到本文件定义的 I/O 范围，并有 round-trip 或独立解析 fixture；
- 几何一致性、索引越界、非法 connectivity、退化单元、字段维度/单位和 unsupported-format 行为有明确测试；
- 至少完成二维与三维、结构与非结构、边界标签、inactive cell、属性字段和 mixed-cell 回归；
- PETSc 可选适配完成，并在官方 Linux runner 上通过至少 2-rank 的 DMPlex/Section/ghost/ownership 集成测试；无 PETSc 的核心构建保持独立；
- 核心测试覆盖 GCC Debug + ASan/UBSan、Clang Release、MSVC Release；PETSc/MPI gate 只运行实际受影响的平台；
- 热路径无每实体堆对象、无隐藏字符串查找、无不必要 connectivity 复制；固定大网格 benchmark 记录访问吞吐和内存基线；
- 最终 diff 审计确认没有引入 flow residual、Darcy flux、thermodynamic closure、solver、runtime 或 frontend 逻辑。

## 9. 本 PR 不做

本 PR 不实现有限体积/有限元残差、Darcy/Forchheimer 通量、井模型、热力学、相态、时间推进、Newton/Krylov 求解、AMR、动态重分区、GPU 数据布局、可视化 UI 或前端编辑器。上述能力只能在网格契约稳定后由独立增量消费。
