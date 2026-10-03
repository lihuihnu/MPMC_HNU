# 外部网格兼容性与二维独立读回

唯一 C++ owner 是本目录的 `mpmc_mesh_external_compatibility` target。原有五路径参数入口继续验证固定上游 deal.II/OPM 样例；新增 `--emit-2d <目录>` 只生成待检文件。`verify_2d_readers.py` 调用这个 producer，再用官方 Gmsh 4.15.2 API 和 VTK 9.7.1 `vtkXMLUnstructuredGridReader` 读取输出，完全不调用 MPMC 解析器。它不是新增 CTest，也不由 mesh core target 重复编译或执行。

## 判据与边界

输入为人工构造的 SI 平面网格：独立三角形、顺时针非矩形凸四边形、三角形/四边形混合共享边。包括稀疏且乱序的顶点/单元 ID、显式/生成 face ID、边界组、重叠区域组、点与单元的标量和多分量字段。每例输出原生 MSH/VTU 及 canonical Gmsh→VTU、VTU→Gmsh，共 12 个文件、6 份转换报告。

- Python 中独立写出输入和解析期望，不读取 producer 生成的期望值。三角形面积/质心为 `1 m², (2/3,1/3) m`；梯形为 `5 m², (19/15,14/15) m`；混合例额外三角形为 `1 m², (3,2/3) m`。从第三方读到的坐标以三角形扇分解重新积分，并与解析值比较。
- Gmsh 验证节点坐标、稳定节点/单元/面 ID、线/三角形/四边形类型、循环连接、共享边、PhysicalNames、边界与多区域组成员。
- VTK 验证坐标、三角形/四边形类型、循环连接、UInt64 单元 ID、PointData/CellData 的字段名、分量数、数值和实体对应。ID 使用整数访问器，避免经过 double。VTU 未编码稳定 vertex/face ID；此处按样例的唯一坐标建立映射，不声称这些 ID 被保留。
- 数值采用绝对/相对容差 `1e-12`（坐标/质心 m、面积 m²、字段按声明单位）；样例量级为 1–300，容差用于十进制文本和浮点积分舍入。类型、整数 ID、连接和分组精确匹配。
- 跨格式输出同时检查实际可读内容与 loss report：Gmsh 明确丢失字段，VTU 明确丢失组和面标签。样例的稀疏 vertex/face ID 无法在当前 VTU 往返中保留，报告还必须含 `vtu.vertex_ids_remapped`、`vtu.face_ids_remapped`；官方 VTK 同时确认文件点顺序保持。无标签/字段、默认编号不误报及同集合错绑的实际往返回归由既有 `mesh.core.exchange_io` 唯一拥有，覆盖二维混合单元与三维金字塔的 16 个场景。VTU 的自定义 `mpmc_*` 单位/来源 XML 属性不是 VTK 标准字段语义，本测试不宣称 VTK 理解这些元数据，也不证明所有转换损失都已穷尽。
- 五个负对照分别移除单元 ID 字段、改错字段值、反转连接、改错组名及破坏 XML；均须被拒绝。破坏 XML 的控制会产生预期的 VTK parser 错误诊断，只有所有拒绝判据满足才打印总 PASS。

范围为现有二维线性 XY 平面、MSH 4.1 ASCII 和单 Piece VTU ASCII。未覆盖三维独立读回、binary/appended/compressed、高阶/曲面、大规模性能、Gmsh GUI/CAD 操作或科学求解验证。旧公共真实样例验证仍保留，不能与本项混称。

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
