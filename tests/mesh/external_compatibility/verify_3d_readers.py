"""Independently read four linear 3D cell families with official Gmsh/VTK.

Synthetic SI solids and analytic volumes are specified here independently of
the C++ producer. No MPMC parser, metric or expected-output manifest is used.
"""
import argparse
from datetime import datetime, timezone
import hashlib
import importlib.metadata
import json
from pathlib import Path
import platform
import subprocess
import sys
import xml.etree.ElementTree as ET

import gmsh
from vtkmodules.vtkCommonCore import vtkVersion
from vtkmodules.vtkFiltersVerdict import vtkCellSizeFilter
from vtkmodules.vtkIOXML import vtkXMLUnstructuredGridReader

import verify_2d_readers as shared
from verify_2d_readers import VerificationError, require, close, check_array, check_report


def case_data(name):
    # Gmsh/VTK linear types, physical vertices, boundary faces, analytic m^3.
    solids = {
        'tetrahedron': (4, 10, [(0,0,0), (2,0,0), (0,3,0), (0,0,4)],
                        [(0,1,2), (0,1,3), (0,2,3), (1,2,3)], 4.),
        'hexahedron': (5, 12, [(0,0,0), (2,0,0), (2,3,0), (0,3,0),
                               (0,0,4), (2,0,4), (2,3,4), (0,3,4)],
                        [(0,1,2,3), (4,5,6,7), (0,1,5,4), (1,2,6,5), (2,3,7,6), (3,0,4,7)], 24.),
        'wedge': (6, 13, [(0,0,0), (2,0,0), (0,3,0), (0,0,4), (2,0,4), (0,3,4)],
                  [(0,1,2), (3,4,5), (0,1,4,3), (1,2,5,4), (2,0,3,5)], 12.),
        'pyramid': (7, 14, [(0,0,0), (2,0,0), (2,3,0), (0,3,0), (1,1.5,4)],
                    [(0,1,2,3), (0,1,4), (1,2,4), (2,3,4), (3,0,4)], 8.),
    }
    gmsh_type, vtk_type, base, faces, volume = solids[name]
    width = len(base)
    sparse = [50,7,90,12,110,19,44,3][:width]
    ids = sparse + [v+1000 for v in sparse]
    coords = base + [(10+2*x, -5+2*y, 2+2*z) for x,y,z in base]
    cells = {31: (tuple(range(width, 2*width)), 8*volume),
             901: (tuple(range(width)), volume)}
    # Handwritten face sets, then the documented generated-ID allocation rule.
    keys = [tuple(sorted(v+block*width for v in face)) for block in (0,1) for face in faces]
    annotations = {tuple(sorted(faces[0])): 701,
                   tuple(sorted(v+width for v in faces[0])): 702}
    generated = 10000
    face_ids = {}
    for key in sorted(keys, key=lambda item: (len(item), item)):
        if key in annotations:
            tag = annotations[key]
        else:
            generated += 1
            tag = generated
        face_ids[tag] = tuple(sorted(ids[v] for v in key))
    return gmsh_type, vtk_type, ids, coords, cells, face_ids


def gmsh_elements(dimension, entity=-1):
    kinds, tags, nodes = gmsh.model.mesh.getElements(dimension, entity)
    widths = {2:3, 3:4} if dimension == 2 else {4:4, 5:8, 6:6, 7:5}
    result = {}
    for kind, block_tags, block_nodes in zip(kinds, tags, nodes, strict=True):
        kind = int(kind)
        require(kind in widths, f'unexpected Gmsh {dimension}D type {kind}')
        width = widths[kind]
        require(len(block_nodes) == width*len(block_tags), 'Gmsh connectivity shape')
        for i, tag in enumerate(block_tags):
            tag = int(tag)
            require(tag not in result, 'duplicate Gmsh element tag')
            result[tag] = kind, tuple(int(v) for v in block_nodes[i*width:(i+1)*width])
    return result


def verify_gmsh(path, name, with_groups):
    kind, _, ids, coords, cells, faces = case_data(name)
    gmsh.clear()
    gmsh.logger.start()
    try:
        gmsh.open(str(path))
        tags, xyz, _ = gmsh.model.mesh.getNodes()
        require(len(xyz) == 3*len(tags), 'Gmsh coordinate count')
        nodes = {int(tag): tuple(float(v) for v in xyz[3*i:3*i+3]) for i, tag in enumerate(tags)}
        require(len(nodes) == len(tags) and set(nodes) == set(ids), 'Gmsh stable node IDs')
        for tag, point in zip(ids, coords, strict=True):
            for a, b in zip(nodes[tag], point, strict=True):
                close(a, b, 'Gmsh point coordinate')
        actual = gmsh_elements(3)
        require(set(actual) == set(cells), 'Gmsh stable volume-cell IDs')
        # Official reference-element quadrature/Jacobian, not MPMC's fixed
        # tetrahedral decomposition. These straight-sided solids integrate exactly.
        integration_points, weights = gmsh.model.mesh.getIntegrationPoints(kind, 'Gauss2')
        require(len(weights) > 0, 'Gmsh integration rule is empty')
        volumes = {}
        for tag, (row, volume) in cells.items():
            require(actual[tag] == (kind, tuple(ids[i] for i in row)), 'Gmsh type/node ordering')
            _, determinants, _ = gmsh.model.mesh.getJacobian(tag, integration_points)
            require(all(value > 0 for value in determinants), 'Gmsh nonpositive Jacobian')
            measured = sum(float(w)*float(d) for w, d in zip(weights, determinants, strict=True))
            close(measured, volume, 'Gmsh analytic volume m^3')
            volumes[str(tag)] = measured
        actual_faces = gmsh_elements(2)
        require(set(actual_faces) == set(faces), 'Gmsh stable surface IDs')
        for tag, row in faces.items():
            require(actual_faces[tag][0] == {3:2, 4:3}[len(row)] and
                    tuple(sorted(actual_faces[tag][1])) == row, 'Gmsh surface type/connectivity')
        expected_groups = {
            (2,11): ('base & inlet', {701}), (2,12): ('translated base', {702}),
            (3,21): ('rock & sand', {31,901}), (3,22): ('selected region', {31}),
        } if with_groups else {}
        require({tuple(map(int, group)) for group in gmsh.model.getPhysicalGroups()} == set(expected_groups),
                'Gmsh physical-group keys')
        for (dim, tag), (label, members) in expected_groups.items():
            require(gmsh.model.getPhysicalName(dim, tag) == label, 'Gmsh physical-group name')
            actual_members = set()
            for entity in gmsh.model.getEntitiesForPhysicalGroup(dim, tag):
                actual_members.update(gmsh_elements(dim, int(entity)))
            require(actual_members == members, 'Gmsh physical-group membership')
        require(len(gmsh.view.getTags()) == 0, 'Gmsh unexpectedly serialized fields')
        require(not any(line.startswith(('Error:', 'Warning:')) for line in gmsh.logger.get()),
                'Gmsh reader warning/error')
        return volumes
    except VerificationError:
        raise
    except Exception as error:
        raise VerificationError(f'Gmsh read failed: {path.name}: {error}') from error
    finally:
        gmsh.logger.stop()


def verify_vtu(path, name, with_fields):
    _, kind, _, coords, cells, _ = case_data(name)
    reader = vtkXMLUnstructuredGridReader()
    errors = []
    reader.AddObserver('ErrorEvent', lambda *_: errors.append('reader error'))
    reader.SetFileName(str(path))
    require(reader.CanReadFile(str(path)) == 1, 'VTK cannot recognize file')
    reader.Update()
    require(not errors and reader.GetErrorCode() == 0, 'VTK parse/pipeline error')
    grid = reader.GetOutput()
    require(grid.GetNumberOfPoints() == len(coords) and grid.GetNumberOfCells() == len(cells), 'VTK counts')
    for i, point in enumerate(coords):
        for a, b in zip(grid.GetPoint(i), point, strict=True):
            close(a, b, 'VTK point coordinate/order')
    pd, cd = grid.GetPointData(), grid.GetCellData()
    identity = cd.GetArray('mpmc_global_cell_id')
    require(identity is not None and identity.IsA('vtkUnsignedLongLongArray'), 'VTK UInt64 cell IDs')
    require(identity.GetNumberOfComponents() == 1 and identity.GetNumberOfTuples() == len(cells), 'VTK ID shape')
    tags = [int(identity.GetValue(i)) for i in range(len(cells))]
    require(len(set(tags)) == len(tags) and set(tags) == set(cells), 'VTK cell ID values')
    for i, tag in enumerate(tags):
        require(grid.GetCellType(i) == kind, 'VTK cell type')
        cell = grid.GetCell(i)
        row = tuple(cell.GetPointId(j) for j in range(cell.GetNumberOfPoints()))
        require(row == cells[tag][0], 'VTK cell node ordering')
    size = vtkCellSizeFilter()
    size.AddObserver('ErrorEvent', lambda *_: errors.append('volume filter error'))
    size.SetInputData(grid)
    size.ComputeVertexCountOff()
    size.ComputeLengthOff()
    size.ComputeAreaOff()
    size.ComputeVolumeOn()
    size.Update()
    require(not errors and size.GetErrorCode() == 0, 'VTK volume filter failed')
    volumes = size.GetOutput().GetCellData().GetArray('Volume')
    check_array(volumes, 1, [(cells[tag][1],) for tag in tags])
    require({pd.GetArrayName(i) for i in range(pd.GetNumberOfArrays())} ==
            ({'temperature','position'} if with_fields else set()), 'VTK point field names')
    require({cd.GetArrayName(i) for i in range(cd.GetNumberOfArrays())} ==
            ({'mpmc_global_cell_id','marker','cell_pair'} if with_fields else {'mpmc_global_cell_id'}),
            'VTK cell field names')
    if with_fields:
        check_array(pd.GetArray('temperature'), 1, [(300.+i,) for i in range(len(coords))])
        check_array(pd.GetArray('position'), 3, coords)
        markers = [3.1 if tag == 31 else 90.1 for tag in tags]
        check_array(cd.GetArray('marker'), 1, [(v,) for v in markers])
        check_array(cd.GetArray('cell_pair'), 2, [(v,-v) for v in markers])
    return {str(tag): volumes.GetValue(i) for i, tag in enumerate(tags)}


def negative_controls(directory):
    output = directory / 'negative-controls'
    output.mkdir()
    controls = []
    for label in ('wrong_type', 'wrong_connectivity', 'wrong_field_association', 'missing_cell_id'):
        tree = ET.parse(directory / 'tetrahedron.vtu')
        arrays = {node.get('Name'): node for node in tree.findall('.//DataArray')}
        if label == 'wrong_type':
            arrays['types'].text = '9 9'  # Valid VTK_QUAD width, wrong dimension/type.
        elif label == 'wrong_connectivity':
            values = arrays['connectivity'].text.split()
            values[0], values[1] = values[1], values[0]
            arrays['connectivity'].text = ' '.join(values)
        elif label == 'wrong_field_association':
            arrays['marker'].text = ' '.join(reversed(arrays['marker'].text.split()))
        else:
            arrays['mpmc_global_cell_id'].set('Name', 'missing_cell_id')
        path = output / (label + '.vtu')
        tree.write(path, encoding='utf-8', xml_declaration=True)
        controls.append((label, path, verify_vtu))
    # Keep group names/tags/counts but swap the two surface memberships.
    lines = (directory / 'tetrahedron.msh').read_text(encoding='utf-8').splitlines()
    start = lines.index('$Entities') + 1
    points, curves, surfaces, _ = map(int, lines[start].split())
    changed = 0
    for i in range(start+1+points+curves, start+1+points+curves+surfaces):
        row = lines[i].split()
        if row[7] == '1' and row[8] in ('11','12'):
            row[8] = '12' if row[8] == '11' else '11'
            lines[i] = ' '.join(row)
            changed += 1
    require(changed == 2, 'Gmsh membership mutation did not apply twice')
    path = output / 'wrong_surface_membership.msh'
    path.write_text('\n'.join(lines)+'\n', encoding='utf-8')
    controls.append(('wrong_surface_membership', path, verify_gmsh))
    text = (directory / 'tetrahedron.vtu').read_text(encoding='utf-8')
    require('</UnstructuredGrid>' in text, 'XML mutation did not apply')
    path = output / 'broken_xml.vtu'
    path.write_text(text.replace('</UnstructuredGrid>', '</WrongGrid>', 1), encoding='utf-8')
    controls.append(('broken_xml', path, verify_vtu))
    rejected = {}
    for label, path, verify in controls:
        try:
            verify(path, 'tetrahedron', True)
        except VerificationError as error:
            rejected[label] = str(error)
            print(f'[PASS] independent.3d.negative.{label}')
        else:
            raise VerificationError(f'corrupt export escaped the oracle: {label}')
    return rejected


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--producer', type=Path, required=True)
    parser.add_argument('--output-dir', type=Path, required=True, help='new directory outside checkout')
    args = parser.parse_args()
    directory = args.output_dir.resolve()
    directory.mkdir(parents=True, exist_ok=False)
    producer = args.producer.resolve(strict=True)
    report = {'started_utc': datetime.now(timezone.utc).isoformat(),
              'platform': platform.platform(), 'python': sys.version,
              'gmsh': gmsh.__version__, 'vtk': vtkVersion.GetVTKVersion(),
              'packages': {name: importlib.metadata.version(name) for name in ('gmsh','vtk','numpy')},
              'producer_sha256': hashlib.sha256(producer.read_bytes()).hexdigest(),
              'oracle_sha256': hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
              'shared_oracle_sha256': hashlib.sha256(Path(shared.__file__).read_bytes()).hexdigest(),
              'passed_files': [], 'volumes_m3': {}, 'status': 'failed'}
    initialized = False
    try:
        run = subprocess.run([str(producer), '--emit-3d', '.'], cwd=directory,
                             capture_output=True, text=True, timeout=60, check=True)
        report['producer_stdout'] = run.stdout
        gmsh.initialize(['mpmc-independent-reader-3d'], readConfigFiles=False)
        initialized = True
        gmsh.option.setNumber('General.Terminal', 0)
        for name in ('tetrahedron','hexahedron','wedge','pyramid'):
            for suffix, verify, payload in (
                ('.msh', verify_gmsh, True), ('.vtu', verify_vtu, True),
                ('_from_vtu.msh', verify_gmsh, False), ('_from_gmsh.vtu', verify_vtu, False),
            ):
                filename = name + suffix
                report['volumes_m3'][filename] = verify(directory / filename, name, payload)
                report['passed_files'].append(filename)
                print(f'[PASS] independent.3d.read.{filename}')
            check_report(directory / (name + '_from_vtu.msh.report'), {'gmsh.fields_not_serialized'})
            check_report(directory / (name + '_from_gmsh.vtu.report'),
                         {'vtu.groups_not_serialized', 'vtu.face_tags_not_serialized',
                          'vtu.vertex_ids_remapped', 'vtu.face_ids_remapped'})
        report['negative_control_rejections'] = negative_controls(directory)
        report['negative_controls'] = len(report['negative_control_rejections'])
        report['conversion_reports_checked'] = 8
        report['status'] = 'passed'
        print('[PASS] independent.mesh.3d_readers files=16 reports=8 negative_controls=6')
    except Exception as error:
        report['error'] = str(error)
        print(f'[FAIL] {error}', file=sys.stderr)
        return 1
    finally:
        if initialized:
            gmsh.finalize()
        report['output_sha256'] = {p.relative_to(directory).as_posix(): hashlib.sha256(p.read_bytes()).hexdigest()
                                   for p in sorted(directory.rglob('*')) if p.is_file()}
        (directory / 'result.json').write_text(json.dumps(report, indent=2)+'\n', encoding='utf-8')
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
