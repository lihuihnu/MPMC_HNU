"""Official writers -> MPMC import/canonical export -> official readers.

Only synthetic linear ASCII fixtures; no MPMC-generated input or golden data.
The .import files expose actual pre-export state, checked against reader data.
"""
import argparse
from datetime import datetime, timezone
import hashlib
import json
import math
from pathlib import Path
import platform
import subprocess
import sys

import gmsh

import verify_3d_readers as reader
from verify_3d_chains import incidence, check_fields, check_report, write_vtu as write_vtu_common
import verify_2d_readers as oracle
from verify_2d_readers import VerificationError, require, close


VTK_KIND = {2:5, 3:9}
GROUPS = {(1,11): ('bottom', {701}), (1,12): ('upper side', {702}),
          (2,21): ('rock & sand', {31,901}), (2,22): ('selected region', {31})}


def fixture(name, disconnected):
    require(name == 'mixed', 'unknown 2D chain fixture')
    ids, points, cells = oracle.case_data(name)
    rows = {tag:list(value[1]) for tag,value in cells.items()}
    if disconnected:
        for old in (1,2):
            rows[31] = [len(points) if v == old else v for v in rows[31]]
            points.append(points[old])
            ids.append(ids[old]+1000)
    return ({tag:value[0] for tag,value in cells.items()}, ids, points, rows,
            {tag:value[2] for tag,value in cells.items()})


def write_gmsh(path, data):
    kinds, ids, points, rows, _ = data
    gmsh.clear()
    gmsh.model.add('independent-mixed-input')
    for entity in (10,20):
        gmsh.model.addDiscreteEntity(2, entity)
    # Deliberately reverse insertion order and use sparse, nonmonotone tags.
    order = list(reversed(range(len(points))))
    gmsh.model.mesh.addNodes(2, 10, [ids[i] for i in order], [x for i in order for x in points[i]])
    for entity, tag in ((10,901), (20,31)):
        gmsh.model.mesh.addElementsByType(entity, kinds[tag], [tag], [ids[i] for i in rows[tag]])
    # Only two exterior edges are explicit; the importer must build the rest.
    for entity, tag, row in ((41,701,rows[901][-2:]), (42,702,rows[31][:2])):
        gmsh.model.addDiscreteEntity(1, entity)
        gmsh.model.mesh.addElementsByType(entity, 1, [tag], [ids[i] for i in row])
    for (dim, tag), (label, _) in GROUPS.items():
        entities = {(1,11):[41], (1,12):[42], (2,21):[10,20], (2,22):[20]}[(dim,tag)]
        gmsh.model.addPhysicalGroup(dim, entities, tag, label)
    gmsh.option.setNumber('Mesh.MshFileVersion', 4.1)
    gmsh.option.setNumber('Mesh.Binary', 0)
    gmsh.option.setNumber('Mesh.SaveAll', 1)
    gmsh.write(str(path))
    return {i: tag for i, tag in enumerate(ids)}


def write_vtu(path, data):
    return write_vtu_common(path, data, VTK_KIND)


def read_official(path):
    if path.suffix == '.vtu':
        grid = reader.read_vtu(path)
        nodes = {i+1: grid.GetPoint(i) for i in range(grid.GetNumberOfPoints())}
        ids = grid.GetCellData().GetArray('mpmc_global_cell_id')
        require(ids is not None and ids.IsA('vtkUnsignedLongLongArray') and
                ids.GetNumberOfComponents() == 1 and ids.GetNumberOfTuples() == grid.GetNumberOfCells(),
                'chain VTU UInt64 cell identity')
        from vtkmodules.vtkFiltersVerdict import vtkCellSizeFilter
        measure = vtkCellSizeFilter()
        measure.SetInputData(grid)
        measure.SetComputeArea(True)
        measure.SetComputeVolume(False)
        measure.Update()
        areas_array = measure.GetOutput().GetCellData().GetArray('Area')
        require(areas_array is not None, 'official VTK area array')
        cells, areas, cell_faces = {}, {}, {}
        tags = []
        for i in range(grid.GetNumberOfCells()):
            tag = int(ids.GetValue(i))
            require(tag not in cells, 'duplicate chain VTU cell ID')
            cell = grid.GetCell(i)
            kind = {v:k for k,v in VTK_KIND.items()}.get(grid.GetCellType(i))
            require(kind is not None, 'unsupported chain VTU type')
            cells[tag] = kind, tuple(cell.GetPointId(j)+1 for j in range(cell.GetNumberOfPoints()))
            cell_faces[tag] = []
            for j in range(cell.GetNumberOfEdges()):
                face = cell.GetEdge(j)
                cell_faces[tag].append(tuple(face.GetPointId(k)+1 for k in range(face.GetNumberOfPoints())))
            areas[tag] = float(areas_array.GetValue(i))
            tags.append(tag)
        fields = {}
        for location, arrays, keys in (('point', grid.GetPointData(), list(nodes)),
                                        ('cell', grid.GetCellData(), tags)):
            for i in range(arrays.GetNumberOfArrays()):
                array = arrays.GetArray(i)
                require(array is not None, 'non-numeric chain field')
                if array.GetName() == 'mpmc_global_cell_id':
                    continue
                require(array.GetNumberOfTuples() == len(keys), 'chain field tuple count')
                fields[location, array.GetName()] = {key: array.GetTuple(j) for j,key in enumerate(keys)}
        return dict(nodes=nodes, cells=cells, areas=areas, cell_faces=cell_faces,
                    faces={}, groups={}, fields=fields)
    gmsh.clear()
    gmsh.logger.start()
    try:
        gmsh.open(str(path))
        tags, xyz, _ = gmsh.model.mesh.getNodes()
        nodes = {int(tag): tuple(float(x) for x in xyz[3*i:3*i+3]) for i,tag in enumerate(tags)}
        require(len(nodes) == len(tags), 'duplicate chain Gmsh node ID')
        cells = oracle.gmsh_elements(2)
        areas, cell_faces = {}, {}
        for tag, (kind, row) in cells.items():
            points, weights = gmsh.model.mesh.getIntegrationPoints(kind, 'Gauss2')
            _, determinants, _ = gmsh.model.mesh.getJacobian(tag, points)
            require(len(weights) > 0 and all(d > 0 for d in determinants), 'chain Gmsh Jacobian')
            areas[tag] = sum(float(w)*float(d) for w,d in zip(weights, determinants, strict=True))
            cell_faces[tag] = [(row[i],row[(i+1)%len(row)]) for i in range(len(row))]
        faces = {tag: tuple(row) for tag,(_,row) in oracle.gmsh_elements(1).items()}
        groups = {}
        for dim, tag in gmsh.model.getPhysicalGroups():
            dim, tag = int(dim), int(tag)
            members = set()
            for entity in gmsh.model.getEntitiesForPhysicalGroup(dim, tag):
                members.update(oracle.gmsh_elements(dim, int(entity)))
            groups[dim,tag] = gmsh.model.getPhysicalName(dim,tag), members
        require(not any(line.startswith(('Warning:', 'Error:')) for line in gmsh.logger.get()),
                'chain Gmsh reader warning/error')
        require(len(gmsh.view.getTags()) == 0, 'unexpected Gmsh field views')
        return dict(nodes=nodes, cells=cells, areas=areas, cell_faces=cell_faces,
                    faces=faces, groups=groups, fields={})
    finally:
        gmsh.logger.stop()


def check_geometry(snapshot, name, data, mapping, disconnected):
    kinds, _, points, rows, areas = data
    require(set(snapshot['nodes']) == set(mapping.values()), 'chain vertex identities/count')
    for i, identity in mapping.items():
        for actual, expected in zip(snapshot['nodes'][identity], points[i], strict=True):
            close(actual, expected, 'chain coordinate m')
    require(set(snapshot['cells']) == set(rows), 'chain cell identities')
    expected_cells = oracle.case_data(name)[2]
    centers = {}
    for tag, row in rows.items():
        kind, actual = snapshot['cells'][tag]
        expected = tuple(mapping[i] for i in row)
        require(kind == kinds[tag] and any(actual == expected[i:]+expected[:i]
                for i in range(len(row))), 'chain cell type/cyclic winding')
        edges = {tuple(sorted((expected[i],expected[(i+1)%len(row)]))) for i in range(len(row))}
        require({tuple(sorted(edge)) for edge in snapshot['cell_faces'][tag]} == edges,
                'chain cell edge connectivity')
        vertices = [snapshot['nodes'][v] for v in actual]
        centers[tag] = oracle.analytic_geometry(vertices, expected_cells[tag])
        close(snapshot['areas'][tag], areas[tag], 'chain analytic cell area m^2')
        if 'centroids' in snapshot:
            for a,b in zip(snapshot['centroids'][tag],centers[tag],strict=True):
                close(a,b,'imported area centroid m')
    close(sum(snapshot['areas'].values()), 6., 'chain total area m^2')
    faces = incidence(snapshot)
    internal = sum(len(owners) == 2 for owners in faces.values())
    require(internal == (0 if disconnected else 1), 'chain shared-edge incidence')
    require(len(faces) == (7 if disconnected else 6), 'chain unique edge count')
    interface_result = None
    if not disconnected:
        key = tuple(sorted((mapping[1],mapping[2])))
        require(faces[key] == {31,901}, 'chain shared-edge identities')
        a,b = [snapshot['nodes'][v] for v in key]
        length = math.hypot(b[0]-a[0],b[1]-a[1])
        close(length,math.sqrt(5),'chain shared-edge length m')
        normals = {}
        for tag in (31,901):
            nx,ny = (b[1]-a[1])/length, (a[0]-b[0])/length
            center = centers[tag]
            if nx*((a[0]+b[0])/2-center[0])+ny*((a[1]+b[1])/2-center[1]) < 0:
                nx,ny = -nx,-ny
            sign = 1 if tag == 901 else -1
            close(nx,sign*2/math.sqrt(5),'chain outward normal x')
            close(ny,sign/math.sqrt(5),'chain outward normal y')
            normals[tag] = (nx,ny)
        interface_result = dict(length_m=length,owner_normals=normals)
    return dict(unique_edges=len(faces),interior_edges=internal,boundary_edges=len(faces)-internal,
                areas_m2=snapshot['areas'],centroids_m=centers,shared_edge=interface_result)


def read_import(path, source):
    snapshot = dict(nodes={}, cells={}, areas={}, centroids={}, cell_faces={}, faces={}, groups={}, fields={})
    face_details, cell_face_ids = {}, {}
    for line in path.read_text(encoding='utf-8').splitlines():
        words = line.split()
        kind, tag = map(int, words[:2])
        target = snapshot[{0:'nodes',3:'cells',2:'faces'}[kind]]
        require(tag not in target, 'duplicate imported entity identity')
        if kind == 0:
            require(len(words) == 4, 'imported vertex record shape')
            target[tag] = (*map(float,words[2:]),0.)
            continue
        width = int(words[2])
        row = tuple(map(int,words[3:3+width]))
        tail = words[3+width:]
        if kind == 3:
            target[tag] = {3:2,4:3}[width], row
            snapshot['areas'][tag] = float(tail[0])
            snapshot['centroids'][tag] = tuple(map(float,tail[1:3]))
            count = int(tail[3])
            require(len(tail) == 4+count, 'imported cell record shape')
            cell_face_ids[tag] = tuple(map(int,tail[4:]))
        else:
            require(width == 2, 'imported edge width')
            target[tag] = row
            count = int(tail[0])
            require(len(tail) == 9+count, 'imported edge record shape')
            face_details[tag] = dict(owners=set(map(int,tail[1:1+count])), boundary=int(tail[1+count]),
                physical=int(tail[2+count]), length=float(tail[3+count]), owner=int(tail[4+count]),
                normal=tuple(map(float,tail[5+count:7+count])),center=tuple(map(float,tail[7+count:])))
    for tag, ids in cell_face_ids.items():
        snapshot['cell_faces'][tag] = [snapshot['faces'][i] for i in ids]
    faces = incidence(snapshot)
    require(len(faces) == len(snapshot['faces']), 'imported duplicate/orphan edges')
    for tag, row in snapshot['faces'].items():
        detail = face_details[tag]
        owners = faces[tuple(sorted(row))]
        require(detail['owners'] == owners, 'imported cell-edge / edge-cell reciprocity')
        require(detail['boundary'] == (int(len(owners) == 1) if source == 'gmsh' else -1), 'imported boundary classification')
        require(detail['owner'] in owners, 'imported geometry owner')
        a,b = [snapshot['nodes'][v] for v in row]
        length = math.hypot(b[0]-a[0],b[1]-a[1])
        close(detail['length'],length,'imported edge length m')
        center = ((a[0]+b[0])/2,(a[1]+b[1])/2)
        for actual,expected in zip(detail['center'],center,strict=True):
            close(actual,expected,'imported edge centroid m')
        nx,ny = (b[1]-a[1])/length,(a[0]-b[0])/length
        cell_center = snapshot['centroids'][detail['owner']]
        if nx*(center[0]-cell_center[0])+ny*(center[1]-cell_center[1]) < 0:
            nx,ny = -nx,-ny
        for actual,expected in zip(detail['normal'],(nx,ny),strict=True):
            close(actual,expected,'imported owner outward normal')
    return snapshot, face_details


def run_case(directory, producer, name, source, disconnected):
    data = fixture(name, disconnected)
    suffix = '.msh' if source == 'gmsh' else '.vtu'
    stem = name+'_'+source+('_disconnected' if disconnected else '')
    path = directory/(stem+'_input'+suffix)
    mapping = (write_gmsh if source == 'gmsh' else write_vtu)(path,data)
    original = read_official(path)
    check_geometry(original,name,data,mapping,disconnected)
    check_fields(original,data,mapping,source == 'vtu')
    require(original['groups'] == (GROUPS if source == 'gmsh' else {}), 'official input physical groups')
    if source == 'gmsh':
        require(set(original['faces']) == {701,702}, 'official input must omit generated faces')
    run = subprocess.run([str(producer),'--convert-2d',source,path.name,stem], cwd=directory,
                         capture_output=True,text=True,timeout=60)
    require(run.returncode == 0, f'MPMC conversion failed: {run.stdout}\n{run.stderr}')
    imported, details = read_import(directory/(stem+'.import'),source)
    imported_result = check_geometry(imported,name,data,mapping,disconnected)
    require(all(detail['physical'] == ({701:11,702:12}.get(tag,0) if source == 'gmsh' else 0)
                for tag,detail in details.items()), 'imported boundary physical tags')
    for tag,row in original['faces'].items():
        require(tag in imported['faces'] and set(row) == set(imported['faces'][tag]), 'imported explicit face IDs')
    results = {}
    for target in ('gmsh','vtu'):
        output = directory/(stem+('.msh' if target == 'gmsh' else '.vtu'))
        actual = read_official(output)
        lookup = {tag:i+1 for i,tag in enumerate(imported['nodes'])}
        target_mapping = mapping if target == 'gmsh' else {i:lookup[tag] for i,tag in mapping.items()}
        results[target] = check_geometry(actual,name,data,target_mapping,disconnected)
        check_fields(actual,data,target_mapping,source == target == 'vtu')
        require(actual['groups'] == (GROUPS if source == target == 'gmsh' else {}), 'output physical groups')
        if target == 'gmsh':
            require(set(actual['faces']) == set(imported['faces']), 'exported imported/generated face identities')
            for tag,row in imported['faces'].items():
                require(set(actual['faces'][tag]) == set(row), 'exported face identity binding')
        codes = set()
        if source == 'gmsh' and target == 'vtu':
            codes = {'vtu.groups_not_serialized','vtu.face_tags_not_serialized',
                     'vtu.vertex_ids_remapped','vtu.face_ids_remapped'}
        if source == 'vtu' and target == 'gmsh':
            codes = {'gmsh.fields_not_serialized'}
        if target == 'vtu':
            # VTU carries neither stable point nor face IDs. Derive the IDs
            # that its documented 1..N import convention assigns to the actual
            # reader connectivity, then compare entity bindings, not ID sets.
            keys = sorted(incidence(actual), key=lambda key:(len(key),key))
            reimported_faces = {key:i+1 for i,key in enumerate(keys)}
            vertex_loss = any(tag != index for tag,index in lookup.items())
            face_loss = any(tag != reimported_faces[tuple(sorted(lookup[v] for v in row))]
                            for tag,row in imported['faces'].items())
            require(vertex_loss == ('vtu.vertex_ids_remapped' in codes), 'actual VTU vertex ID loss differs')
            require(face_loss == ('vtu.face_ids_remapped' in codes), 'actual VTU face ID loss differs')
            results[target]['identity_losses'] = {'vertices':vertex_loss,'faces':face_loss}
        check_report(Path(str(output)+'.report'),codes)
        if disconnected:
            try:
                check_geometry(actual,name,data,target_mapping,False)
            except VerificationError as error:
                require(str(error) == 'chain shared-edge incidence', 'control failed for unrelated reason')
                results[target]['connected_contract_rejection'] = str(error)
            else:
                raise VerificationError('disconnected input was silently welded')
    print(f'[PASS] independent.2d.chain.{stem}')
    return {'input':path.name, 'imported':imported_result, 'outputs':results, 'converter_stdout':run.stdout}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--producer',type=Path,required=True)
    parser.add_argument('--output-dir',type=Path,required=True)
    args = parser.parse_args()
    directory = args.output_dir.resolve()
    directory.mkdir(parents=True,exist_ok=False)
    producer = args.producer.resolve(strict=True)
    report = dict(status='failed',started_utc=datetime.now(timezone.utc).isoformat(),
                  platform=platform.platform(),python=sys.version,gmsh=gmsh.__version__,
                  vtk=oracle.vtkVersion.GetVTKVersion(),
                  producer_sha256=hashlib.sha256(producer.read_bytes()).hexdigest(),
                  scripts_sha256={Path(p).name:hashlib.sha256(Path(p).read_bytes()).hexdigest()
                                  for p in (__file__,oracle.__file__,reader.__file__,write_vtu_common.__code__.co_filename)},cases={})
    gmsh.initialize(['mpmc-independent-chain-2d'],readConfigFiles=False)
    gmsh.option.setNumber('General.Terminal',0)
    try:
        for name in ('mixed',):
            for source in ('gmsh','vtu'):
                for disconnected in (False,True):
                    key=name+'_'+source+('_disconnected' if disconnected else '')
                    report['cases'][key]=run_case(directory,producer,name,source,disconnected)
        report.update(status='passed',inputs=2,chains=4,conversion_reports=4,
                      disconnected_inputs=2,disconnected_chains=4,disconnected_reports=4,import_snapshots=4)
        print('[PASS] independent.mesh.2d_chains inputs=2 chains=4 reports=4 disconnected_chains=4')
    except Exception as error:
        report['error']=str(error)
        print(f'[FAIL] {error}',file=sys.stderr)
        return 1
    finally:
        gmsh.finalize()
        report['output_sha256']={p.relative_to(directory).as_posix():hashlib.sha256(p.read_bytes()).hexdigest()
                                 for p in sorted(directory.rglob('*')) if p.is_file()}
        (directory/'result.json').write_text(json.dumps(report,indent=2)+'\n',encoding='utf-8')
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
