"""Official VTK group input -> MPMC import snapshot -> canonical VTU -> VTK.

Synthetic metadata tests include overlapping groups, empty groups, UTF-8 names,
full-width stable identities and same tags on different entity kinds.
"""
import argparse
import copy
from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path
import platform
import subprocess
import sys

from vtkmodules.vtkCommonCore import vtkVersion
from vtkmodules.vtkIOXML import vtkXMLUnstructuredGridWriter
import verify_2d_chains as two
import verify_3d_chains as three
from verify_3d_readers import read_vtu
from vtu_groups import read_groups, write_groups


def require(value, message):
    if not value:
        raise AssertionError(message)


def write(grid, path):
    writer = vtkXMLUnstructuredGridWriter()
    errors = []
    writer.AddObserver('ErrorEvent', lambda *_: errors.append('VTK writer error'))
    writer.SetFileName(str(path))
    writer.SetInputData(grid)
    writer.SetDataModeToAscii()
    writer.SetCompressorTypeToNone()
    require(writer.Write() == 1 and not errors and writer.GetErrorCode() == 0, 'VTK group write failure')


def convert(producer, dim, path, stem):
    return subprocess.run([str(producer), '--convert-vtu-groups', str(dim), path.name, stem],
                          cwd=path.parent, text=True, capture_output=True, timeout=60)


def case(directory, producer, dim, name):
    module = two if dim == 2 else three
    data = module.fixture(name, False)
    kinds, ids, points, rows, measures = data
    ids = [0, 2**64-1] + [2**53+i for i in range(len(ids)-2)]
    data = kinds, ids, points, rows, measures
    faces = three.face_identity_fixture(data, dim)
    stem = f'{dim}d_{name}'
    path = directory / (stem + '_input.vtu')
    mapping = three.write_vtu(path, data, two.VTK_KIND if dim == 2 else three.VTK_KIND,
                              stable_ids=True, face_ids=faces)
    groups = [dict(location=2, dimension=dim-1, tag=17, name='入口 & <wall> "A"', members=[0, 2**64-1]),
              dict(location=2, dimension=dim-1, tag=19, name='overlapping interface', members=[0]),
              dict(location=3, dimension=dim, tag=17, name='rock', members=[31, 901]),
              dict(location=3, dimension=dim, tag=23, name='selected material', members=[31]),
              dict(location=0, dimension=0, tag=2**32-1, name='node set', members=[ids[0], ids[1]]),
              dict(location=1, dimension=1, tag=29, name='empty edges', members=[]),
              dict(location=3, dimension=dim, tag=31, name='', members=[])]
    grid = read_vtu(path)
    write_groups(grid, groups)
    write(grid, path)
    require(read_groups(read_vtu(path)) == groups, 'official input groups')
    run = convert(producer, dim, path, stem)
    require(run.returncode == 0, 'group conversion failed: ' + run.stdout + run.stderr)
    require(json.loads((directory/(stem+'.groups.json')).read_text(encoding='utf-8')) == groups,
            'imported group snapshot mismatch')
    output = directory/(stem+'.vtu')
    require(read_groups(read_vtu(output)) == groups, 'exported group binding mismatch')
    require((directory/(stem+'.report')).read_text().splitlines() == ['lossless'], 'group lossless report')
    geometry = module.check_geometry(module.read_official(output), name, data, mapping, False)
    controls = []
    for label in ('name_binding', 'face_membership', 'cell_membership', 'missing_member', 'schema_version'):
        modified = copy.deepcopy(groups)
        if label == 'name_binding':
            modified[0]['name'], modified[1]['name'] = modified[1]['name'], modified[0]['name']
        elif label == 'face_membership':
            modified[0]['members'], modified[1]['members'] = modified[1]['members'], modified[0]['members']
        elif label == 'cell_membership':
            modified[3]['members'] = [901]
        elif label == 'missing_member':
            unknown = next(value for value in range(1, 1000) if value not in ids)
            modified[4]['members'] = [unknown]
        bad = read_vtu(output)
        write_groups(bad, modified)
        if label == 'schema_version':
            bad.GetFieldData().GetArray('mpmc_group_schema_version').SetValue(0, 2)
        bad_path = directory/(stem+'_'+label+'.vtu')
        write(bad, bad_path)
        loaded = read_vtu(bad_path)  # Structural XML/VTK validity is required.
        if label in ('missing_member', 'schema_version'):
            rejected = convert(producer, dim, bad_path, stem+'_'+label+'_rejected')
            require(rejected.returncode != 0, 'MPMC accepted invalid group ' + label)
            controls.append(dict(case=label, diagnostic=rejected.stderr))
        else:
            require(read_groups(loaded) != groups, 'binding oracle missed ' + label)
            controls.append(dict(case=label, diagnostic='group semantic binding mismatch'))
    return dict(groups=groups, geometry=geometry, negative_controls=controls)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--producer', required=True, type=Path)
    parser.add_argument('--output-dir', required=True, type=Path)
    args = parser.parse_args()
    directory = args.output_dir.resolve()
    directory.mkdir(parents=True, exist_ok=False)
    producer = args.producer.resolve(strict=True)
    report = dict(status='failed', started_utc=datetime.now(timezone.utc).isoformat(),
                  python=sys.version, platform=platform.platform(), vtk=vtkVersion.GetVTKVersion(),
                  producer_sha256=hashlib.sha256(producer.read_bytes()).hexdigest(),
                  oracle_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(), cases={})
    try:
        for dim, name in ((2, 'mixed'), (3, 'hexa_pyramid'), (3, 'tetra_wedge')):
            report['cases'][f'{dim}d_{name}'] = case(directory, producer, dim, name)
        report.update(status='passed', chains=3, conversion_reports=3, negative_controls=15)
        print('[PASS] independent.mesh.vtu_group_chains chains=3 reports=3 negative_controls=15')
    except Exception as error:
        report['error'] = str(error)
        print('[FAIL] ' + str(error), file=sys.stderr)
        return 1
    finally:
        report['output_sha256'] = {p.name: hashlib.sha256(p.read_bytes()).hexdigest()
                                   for p in sorted(directory.iterdir()) if p.is_file()}
        (directory/'result.json').write_text(json.dumps(report, indent=2, ensure_ascii=False)+'\n', encoding='utf-8')
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
