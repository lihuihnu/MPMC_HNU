"""Read MPMC exports through official Gmsh/VTK APIs, not MPMC parsers.

Synthetic SI inputs and analytic expectations are authored here independently
of the C++ producer. No expected metrics/connectivity are read from its output.
Run in the pinned requirements-readers.txt environment (Python >= 3.12).
"""
import argparse
from datetime import datetime, timezone
import hashlib
import importlib.metadata
import json
import math
from pathlib import Path
import platform
import subprocess
import sys

import gmsh
from vtkmodules.vtkCommonCore import vtkVersion
from vtkmodules.vtkIOXML import vtkXMLUnstructuredGridReader


class VerificationError(RuntimeError):
    pass


def require(condition, message):
    if not condition:
        raise VerificationError(message)


def close(actual, expected, message):
    require(math.isfinite(actual) and math.isclose(actual, expected, rel_tol=1e-12, abs_tol=1e-12),
            f'{message}: {actual} != {expected}')


def case_data(name):
    # (type, cyclic local vertices, area, area centroid). Quad is clockwise.
    triangle = (2, (0, 1, 2), 1., (2/3, 1/3))
    quad = (3, (3, 2, 1, 0), 5., (19/15, 14/15))
    if name == 'triangle':
        return [50, 7, 90], [(0., 0., 0.), (2., 0., 0.), (0., 1., 0.)], {901: triangle}
    coords = [(0., 0., 0.), (3., 0., 0.), (2., 2., 0.), (0., 2., 0.)]
    if name == 'quad':
        return [50, 7, 90, 12], coords, {901: quad}
    return [50, 7, 90, 12, 110], coords + [(4., 0., 0.)], {
        31: (2, (1, 4, 2), 1., (3., 2/3)), 901: quad}


def analytic_geometry(points, expected):
    # Independent triangle-fan integration, not the producer's shoelace code.
    origin = points[0]
    area = mx = my = 0.
    for i in range(1, len(points)-1):
        b, c = points[i:i+2]
        signed = ((b[0]-origin[0])*(c[1]-origin[1]) -
                  (b[1]-origin[1])*(c[0]-origin[0])) / 2
        area += signed
        mx += signed * (origin[0]+b[0]+c[0]) / 3
        my += signed * (origin[1]+b[1]+c[1]) / 3
    close(abs(area), expected[2], 'analytic cell area')
    require(area != 0., 'zero signed area')
    close(mx/area, expected[3][0], 'analytic centroid x')
    close(my/area, expected[3][1], 'analytic centroid y')
    return (mx/area, my/area)


def expected_faces(name, ids, cells):
    keys = {tuple(sorted((ids[row[j]], ids[row[(j+1) % len(row)]])))
            for _, row, _, _ in cells.values() for j in range(len(row))}
    annotations = {tuple(sorted((ids[0], ids[1]))): 701,
                   tuple(sorted((ids[1], ids[4 if name == 'mixed' else 2]))): 702}
    result = {}
    generated = 10000
    for key in sorted(keys):
        if key in annotations:
            result[annotations[key]] = key
        else:
            generated += 1
            result[generated] = key
    return result


def gmsh_elements(dimension, entity=-1):
    kinds, tags, nodes = gmsh.model.mesh.getElements(dimension, entity)
    result = {}
    for kind, block_tags, block_nodes in zip(kinds, tags, nodes, strict=True):
        kind = int(kind)
        require(kind in {1, 2, 3}, f'unexpected Gmsh element type {kind}')
        width = {1: 2, 2: 3, 3: 4}[kind]
        require(len(block_nodes) == width*len(block_tags), 'Gmsh connectivity width')
        for i, tag in enumerate(block_tags):
            tag = int(tag)
            require(tag not in result, 'duplicate Gmsh element tag')
            result[tag] = kind, tuple(int(v) for v in block_nodes[i*width:(i+1)*width])
    return result


def verify_gmsh(path, name, with_groups):
    ids, coords, expected_cells = case_data(name)
    gmsh.clear()
    gmsh.logger.start()
    try:
        gmsh.open(str(path))
        node_tags, xyz, _ = gmsh.model.mesh.getNodes()
        require(len(xyz) == 3*len(node_tags), 'Gmsh node coordinate count')
        actual_nodes = {int(tag): tuple(float(v) for v in xyz[3*i:3*i+3])
                        for i, tag in enumerate(node_tags)}
        require(len(actual_nodes) == len(node_tags) and set(actual_nodes) == set(ids), 'Gmsh stable node IDs')
        for tag, point in zip(ids, coords, strict=True):
            for actual, expected in zip(actual_nodes[tag], point, strict=True):
                close(actual, expected, f'Gmsh node {tag}')
        actual_cells = gmsh_elements(2)
        require(set(actual_cells) == set(expected_cells), 'Gmsh cell IDs')
        for tag, expected in expected_cells.items():
            require(actual_cells[tag] == (expected[0], tuple(ids[i] for i in expected[1])),
                    f'Gmsh type/cyclic connectivity for {tag}')
            analytic_geometry([actual_nodes[v] for v in actual_cells[tag][1]], expected)
        actual_faces = gmsh_elements(1)
        wanted_faces = expected_faces(name, ids, expected_cells)
        require(set(actual_faces) == set(wanted_faces), 'Gmsh explicit/generated face IDs')
        for tag, key in wanted_faces.items():
            require(actual_faces[tag][0] == 1 and tuple(sorted(actual_faces[tag][1])) == key,
                    f'Gmsh face connectivity {tag}')
        expected_groups = {
            (1, 11): ('inlet & lower', {701}), (1, 12): ('outlet upper', {702}),
            (2, 21): ('rock & sand', set(expected_cells)),
            (2, 22): ('selected region', {next(iter(expected_cells))}),
        } if with_groups else {}
        groups = {tuple(map(int, item)) for item in gmsh.model.getPhysicalGroups()}
        require(groups == set(expected_groups), 'Gmsh physical-group keys')
        for (dim, tag), (label, members) in expected_groups.items():
            require(gmsh.model.getPhysicalName(dim, tag) == label, 'Gmsh physical-group name')
            actual = set()
            for entity in gmsh.model.getEntitiesForPhysicalGroup(dim, tag):
                actual.update(gmsh_elements(dim, int(entity)))
            require(actual == members, f'Gmsh physical-group membership {dim}:{tag}')
        require(len(gmsh.view.getTags()) == 0, 'Gmsh fields unexpectedly serialized')
        require(not any(line.startswith(('Error:', 'Warning:')) for line in gmsh.logger.get()),
                'Gmsh reader warning/error')
    except VerificationError:
        raise
    except Exception as error:
        raise VerificationError(f'Gmsh read failed: {path.name}: {error}') from error
    finally:
        gmsh.logger.stop()


def check_array(array, width, expected):
    require(array is not None, 'missing VTK field')
    require(array.GetNumberOfComponents() == width and array.GetNumberOfTuples() == len(expected),
            'VTK field shape')
    for row, values in enumerate(expected):
        for component, value in enumerate(values):
            close(array.GetComponent(row, component), value, 'VTK field value/association')


def vtk_vertex_ids(grid, required=True):
    array = grid.GetPointData().GetArray('mpmc_global_vertex_id')
    if array is None and not required:
        return list(range(1,grid.GetNumberOfPoints()+1))
    require(array is not None and array.IsA('vtkUnsignedLongLongArray'), 'VTK UInt64 vertex IDs required')
    require(array.GetNumberOfComponents() == 1 and array.GetNumberOfTuples() == grid.GetNumberOfPoints(),
            'VTK vertex ID shape')
    ids = [int(array.GetValue(i)) for i in range(grid.GetNumberOfPoints())]
    require(len(set(ids)) == len(ids), 'VTK duplicate vertex identity')
    return ids


FACE_ARRAYS = ('mpmc_global_face_id','mpmc_face_vertex_offsets','mpmc_face_vertex_ids')


def vtk_face_ids(grid, required=True):
    """Typed official VTK arrays, checked against native VTK cell faces/edges."""
    arrays = [grid.GetFieldData().GetArray(name) for name in FACE_ARRAYS]
    if all(array is None for array in arrays) and not required:
        return {}
    values = []
    for array in arrays:
        require(array is not None and array.IsA('vtkUnsignedLongLongArray') and
                array.GetNumberOfComponents() == 1, 'VTK UInt64 face table required')
        values.append([int(array.GetValue(i)) for i in range(array.GetNumberOfTuples())])
    ids, offsets, vertices = values
    require(len(ids) == len(offsets) and len(set(ids)) == len(ids), 'VTK face ID/offset count or uniqueness')
    point_ids = vtk_vertex_ids(grid,required=False)
    native = set()
    for i in range(grid.GetNumberOfCells()):
        cell = grid.GetCell(i)
        dim = cell.GetCellDimension()
        require(dim in (2,3), 'VTK face identity cell dimension')
        for j in range(cell.GetNumberOfEdges() if dim == 2 else cell.GetNumberOfFaces()):
            face = cell.GetEdge(j) if dim == 2 else cell.GetFace(j)
            native.add(tuple(sorted(point_ids[face.GetPointId(k)] for k in range(face.GetNumberOfPoints()))))
    faces, begin = {}, 0
    for tag, end in zip(ids,offsets,strict=True):
        require(begin < end <= len(vertices), 'VTK face offset range')
        row = tuple(sorted(vertices[begin:end]))
        require(len(set(row)) == len(row) and row in native, 'VTK face table references actual topology')
        faces[tag] = row
        begin = end
    require(begin == len(vertices) and len(set(faces.values())) == len(faces) and
            set(faces.values()) == native, 'VTK complete unique face coverage')
    return faces


FACE_TAG_ARRAY = 'mpmc_face_physical_tag'


def vtk_face_tags(grid, required=True):
    faces = vtk_face_ids(grid,required=required)
    array = grid.GetFieldData().GetArray(FACE_TAG_ARRAY)
    if array is None and not required:
        return {tag:0 for tag in faces}
    require(array is not None and array.IsA('vtkUnsignedIntArray') and
            array.GetNumberOfComponents() == 1 and array.GetNumberOfTuples() == len(faces),
            'VTK scalar UInt32 face physical tags required')
    tags = {tag:int(array.GetValue(i)) for i,tag in enumerate(faces)}
    point_ids = vtk_vertex_ids(grid,required=False)
    supports = {}
    for i in range(grid.GetNumberOfCells()):
        cell = grid.GetCell(i)
        dim = cell.GetCellDimension()
        for j in range(cell.GetNumberOfEdges() if dim == 2 else cell.GetNumberOfFaces()):
            face = cell.GetEdge(j) if dim == 2 else cell.GetFace(j)
            key = tuple(sorted(point_ids[face.GetPointId(k)] for k in range(face.GetNumberOfPoints())))
            supports[key] = supports.get(key,0)+1
    require(all(supports[row] == 1 or tags[tag] == 0 for tag,row in faces.items()),
            'VTK interior face must be untagged')
    return tags


def verify_vtu(path, name, with_fields):
    vertex_ids, coords, expected_cells = case_data(name)
    reader = vtkXMLUnstructuredGridReader()
    errors = []
    reader.AddObserver('ErrorEvent', lambda *_: errors.append('ErrorEvent'))
    reader.SetFileName(str(path))
    require(reader.CanReadFile(str(path)) == 1, 'VTK cannot recognize file')
    reader.Update()
    require(not errors and reader.GetErrorCode() == 0, 'VTK parse/pipeline error')
    grid = reader.GetOutput()
    require(grid.GetNumberOfPoints() == len(coords) and grid.GetNumberOfCells() == len(expected_cells),
            'VTK mesh counts')
    require(vtk_face_ids(grid) == expected_faces(name,vertex_ids,expected_cells), "VTK stable edge bindings")
    require(vtk_face_tags(grid) == {tag:({701:11,702:12}.get(tag,0) if not with_fields else 0)
            for tag in expected_faces(name,vertex_ids,expected_cells)}, 'VTK physical edge tag binding')
    require(vtk_vertex_ids(grid) == vertex_ids, 'VTK stable vertex IDs/order')
    actual_points = [grid.GetPoint(i) for i in range(grid.GetNumberOfPoints())]
    source_indices = []
    for point in actual_points:
        candidates = [i for i, expected in enumerate(coords)
                      if all(math.isclose(a, b, rel_tol=1e-12, abs_tol=1e-12)
                             for a, b in zip(point, expected, strict=True))]
        require(len(candidates) == 1, 'VTK point coordinate/identity mapping')
        source_indices.append(candidates[0])
    require(len(set(source_indices)) == len(coords), 'VTK duplicate/missing point')
    # Conversion identity analysis relies on the writer retaining file point
    # order; official VTK must observe that same order, even for sparse tags.
    require(source_indices == list(range(len(coords))), 'VTK file point order changed')
    pd, cd = grid.GetPointData(), grid.GetCellData()
    id_array = cd.GetArray('mpmc_global_cell_id')
    require(id_array is not None and id_array.IsA('vtkUnsignedLongLongArray'), 'VTK UInt64 cell IDs required')
    require(id_array.GetNumberOfComponents() == 1 and id_array.GetNumberOfTuples() == len(expected_cells),
            'VTK cell ID shape')
    # Typed integer accessor: GetTuple/GetComponent would pass IDs through double.
    cell_ids = [int(id_array.GetValue(i)) for i in range(grid.GetNumberOfCells())]
    require(len(set(cell_ids)) == len(cell_ids) and set(cell_ids) == set(expected_cells), 'VTK cell ID values')
    for i, tag in enumerate(cell_ids):
        expected = expected_cells[tag]
        require(grid.GetCellType(i) == {2: 5, 3: 9}[expected[0]], 'VTK cell type')
        cell = grid.GetCell(i)
        local_ids = [cell.GetPointId(j) for j in range(cell.GetNumberOfPoints())]
        require(tuple(source_indices[j] for j in local_ids) == expected[1], 'VTK cyclic connectivity')
        analytic_geometry([actual_points[j] for j in local_ids], expected)
    require({pd.GetArrayName(i) for i in range(pd.GetNumberOfArrays())} ==
            ({'temperature', 'point_vector', 'mpmc_global_vertex_id'} if with_fields else {'mpmc_global_vertex_id'}), 'VTK point field names')
    require({cd.GetArrayName(i) for i in range(cd.GetNumberOfArrays())} ==
            ({'mpmc_global_cell_id', 'marker', 'cell_pair'} if with_fields else {'mpmc_global_cell_id'}),
            'VTK cell field names')
    if with_fields:
        check_array(pd.GetArray('temperature'), 1, [(300.+j,) for j in source_indices])
        check_array(pd.GetArray('point_vector'), 3,
                    [(coords[j][0], coords[j][1], coords[j][0]-2*coords[j][1]) for j in source_indices])
        markers = [3.1 if tag == 31 else 90.1 for tag in cell_ids]
        check_array(cd.GetArray('marker'), 1, [(value,) for value in markers])
        check_array(cd.GetArray('cell_pair'), 2, [(value, -value) for value in markers])


def check_report(path, expected_codes):
    lines = path.read_text(encoding='utf-8').splitlines()
    require(lines and lines[0] == 'lossy' and set(lines[1:]) == expected_codes,
            f'conversion loss report mismatch: {path.name}: {lines}')


def negative_controls(directory):
    # Well-formed but wrong outputs must fail too; opening a file is not a pass.
    controls = [
        ('missing_cell_id', 'triangle.vtu', 'mpmc_global_cell_id', 'missing_cell_id', verify_vtu),
        ('wrong_field', 'triangle.vtu', '300 301 302', '999 301 302', verify_vtu),
        ('wrong_connectivity', 'triangle.vtu', '0 1 2', '0 2 1', verify_vtu),
        ('wrong_group_name', 'triangle.msh', 'rock & sand', 'wrong region', verify_gmsh),
        ('broken_xml', 'triangle.vtu', '</UnstructuredGrid>', '</WrongGrid>', verify_vtu),
    ]
    output = directory / 'negative-controls'
    output.mkdir()
    for label, source, old, new, verify in controls:
        text = (directory / source).read_text(encoding='utf-8')
        require(old in text, f'negative-control mutation did not apply: {label}')
        path = output / (label + Path(source).suffix)
        path.write_text(text.replace(old, new, 1), encoding='utf-8')
        try:
            verify(path, 'triangle', True)
        except VerificationError:
            print(f'[PASS] independent.negative.{label}')
        else:
            raise VerificationError(f'corrupt export escaped the oracle: {label}')
    return len(controls)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--producer', type=Path, required=True)
    parser.add_argument('--output-dir', type=Path, required=True,
                        help='new directory outside the checkout; retains files and result.json')
    args = parser.parse_args()
    directory = args.output_dir.resolve()
    directory.mkdir(parents=True, exist_ok=False)
    producer = args.producer.resolve(strict=True)
    report = {'started_utc': datetime.now(timezone.utc).isoformat(),
              'platform': platform.platform(), 'python': sys.version,
              'gmsh': gmsh.__version__, 'vtk': vtkVersion.GetVTKVersion(),
              'packages': {name: importlib.metadata.version(name) for name in ('gmsh', 'vtk', 'numpy')},
              'producer_sha256': hashlib.sha256(producer.read_bytes()).hexdigest(),
              'oracle_sha256': hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
              'passed_files': [], 'status': 'failed'}
    initialized = False
    try:
        # Relative output keeps the C++ argv filesystem boundary portable on Windows.
        run = subprocess.run([str(producer), '--emit-2d', '.'], cwd=directory,
                             capture_output=True, text=True, timeout=60, check=True)
        report['producer_stdout'] = run.stdout
        gmsh.initialize(['mpmc-independent-reader'], readConfigFiles=False)
        initialized = True
        gmsh.option.setNumber('General.Terminal', 0)
        for name in ('triangle', 'quad', 'mixed'):
            for suffix, check_file, original_payload in (
                ('.msh', verify_gmsh, True), ('.vtu', verify_vtu, True),
                ('_from_vtu.msh', verify_gmsh, False), ('_from_gmsh.vtu', verify_vtu, False),
            ):
                filename = name + suffix
                check_file(directory / filename, name, original_payload)
                report['passed_files'].append(filename)
                print(f'[PASS] independent.read.{filename}')
            check_report(directory / (name + '_from_vtu.msh.report'), {'gmsh.fields_not_serialized'})
            check_report(directory / (name + '_from_gmsh.vtu.report'),
                         {'vtu.groups_not_serialized'})
        report['negative_controls'] = negative_controls(directory)
        report['conversion_reports_checked'] = 6
        report['status'] = 'passed'
        print('[PASS] independent.mesh.2d_readers files=12 reports=6 negative_controls=5')
    except Exception as error:
        report['error'] = str(error)
        print(f'[FAIL] {error}', file=sys.stderr)
        return 1
    finally:
        if initialized:
            gmsh.finalize()
        report['output_sha256'] = {p.name: hashlib.sha256(p.read_bytes()).hexdigest()
                                   for p in sorted(directory.iterdir()) if p.is_file()}
        (directory / 'result.json').write_text(json.dumps(report, indent=2) + '\n', encoding='utf-8')
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
