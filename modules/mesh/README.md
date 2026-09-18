# Mesh 模块：基础契约

`mpmc::mesh` 面向后续多相多组分流动离散，负责网格拓扑、几何、字段、求解自由度布局、文件 I/O 与并行分区元数据。网格层不得依赖 thermodynamics、flash、physics、runtime、前端或具体流动方程；PETSc/MPI 只允许出现在可选适配层，公共核心头文件不得泄漏 PETSc 类型。

> 当前状态：core topology/index 已形成不可变 `Topology` snapshot：统一拥有按 entity kind 分组的 global-ID 数组与实际 `CsrAdjacency` relations；实体 count 直接由 ID 数组长度确定，构造时校验同类 global-ID 唯一性、relation source/target count 和 kind-pair 唯一性。`EntityKind`、64-bit `GlobalEntityId`、32-bit `LocalIndex` 与紧凑 CSR 仍是底层表示；几何、字段、I/O、DoF、partition 与 PETSc 仍未实现。

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

## 3. 几何与字段

几何层至少提供：

- vertex coordinates；
- cell centroid 与 measure（1D 长度、2D 面积、3D 体积）；
- face centroid、measure 和相对 owner cell 的一致方向 normal；
- boundary/interface 标识；
- 非有限坐标、重复/越界 connectivity、零或负 measure、非法 orientation 的确定性错误。

字段系统必须记录 entity location、component count、数值类型语义与单位/来源元数据，不允许仅靠字符串猜测布局。首批标准科研字段：

- porosity：cell scalar，标准无量纲；
- permeability：cell scalar/diagonal/symmetric tensor，标准 SI 为 m^2；
- conductivity：cell/face 可配置 scalar 或 tensor，单位必须由调用方显式声明；
- 用户自定义 cell/face/vertex field。

I/O 遇到未知单位、缺少必要分量或数组长度不匹配时不得静默补值或重排。

## 4. 网格来源与文件 I/O

本 PR 合并前至少覆盖以下实用入口：

1. 内建 Cartesian 生成：1D/2D/3D，支持 uniform 与按轴非均匀坐标；
2. Gmsh MSH：至少支持当前常用 4.1 网格的导入，并能导出仓库支持的实体、physical tags 与字段子集；
3. VTK UnstructuredGrid：至少一种标准 VTK/VTU 路径可完成几何、拓扑和 cell/point fields 的 round-trip；
4. reservoir corner-point：至少支持 Eclipse 风格 GRDECL 的核心 `SPECGRID/COORD/ZCORN/ACTNUM` 导入，并能读取常用 `PORO/PERMX/PERMY/PERMZ` 属性。

若某格式只实现声明的子集，必须在解析器入口与文档中显式拒绝未支持特性；不得“读取成功”后丢失高阶节点、物理组、inactive cell 或字段。

## 5. 求解变量与拓扑索引

`DofLayout`（最终名称可在实现审计时调整）必须能够把一个或多个变量绑定到 cell/face/edge/vertex，并给出：

- 每实体 DoF 数；
- local contiguous offset；
- owned/ghost 视图；
- global DoF numbering；
- field/component 到标量槽位的 O(1) 或摊还 O(1) 映射；
- 只读 span/view 访问，避免热循环临时分配。

该层只管理布局，不知道 pressure、saturation、composition 等具体物理意义。

## 6. 并行与 PETSc 边界

核心 `mpmc::mesh` 保持 MPI/PETSc 可选。并行语义必须先在仓库自己的数据结构中成立：partition vector、owner rank、local owned/ghost numbering、overlap/halo、跨 rank 共享实体和一致 global ID。

可选 `mesh_petsc` 适配层负责：

- 从核心拓扑创建或填充 DMPlex；
- 把 DoF layout 映射到 PetscSection；
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
