"""Opt-in audit of existing MRST folder exports; never modify source data.

The authoritative new producer is export_mpmc_mesh.m. This migration/audit
reader consumes the user's historical numeric folder format, records hashes,
and checks every HDF5 dataset independently with h5py. No data is vendored.
"""
import argparse
import hashlib
import json
import subprocess
import time
from pathlib import Path

import h5py
import numpy as np
import psutil


def compare_hdf5(a, b):
    with h5py.File(a) as left, h5py.File(b) as right:
        names = []
        left.visititems(lambda name, obj: names.append(name) if isinstance(obj, h5py.Dataset) else None)
        other = []
        right.visititems(lambda name, obj: other.append(name) if isinstance(obj, h5py.Dataset) else None)
        assert names == other
        for name in names:
            x, y = left[name], right[name]
            assert x.shape == y.shape and x.dtype == y.dtype, name
            for start in range(0, x.size, 1048576):
                assert np.array_equal(x[start:start+1048576], y[start:start+1048576]), name
    return len(names)


def migrate(source, output):
    consumed = {}

    def read(name):
        path = source / name / 'data.csv'
        digest = hashlib.sha256()
        with path.open('rb') as stream:
            for block in iter(lambda: stream.read(1024*1024), b''):
                digest.update(block)
        consumed[name] = digest.hexdigest()
        return np.loadtxt(path, delimiter=',', ndmin=1)

    def integers(values, subtract=0):
        assert np.all(np.isfinite(values)) and np.all(values == np.floor(values))
        assert np.all((values >= subtract) & (values <= 2**53))
        return values.astype('uint64') - np.uint64(subtract)

    dim = int(read('griddim')[0])
    coords = read('nodes/coords')
    if dim == 2:
        coords = np.column_stack([coords, np.zeros(len(coords))])
    nf, nc = int(read('faces/num')[0]), int(read('cells/num')[0])
    fo, fn = integers(read('faces/nodePos'), 1), integers(read('faces/nodes'), 1)
    co = integers(read('cells/facePos'), 1)
    cell_faces = read('cells/faces')
    cf = integers(cell_faces[:, 0] if cell_faces.ndim == 2 else cell_faces, 1)
    neighbors = integers(read('faces/neighbors'))
    owner = np.repeat(np.arange(1, nc+1, dtype='uint64'), np.diff(co).astype('int64'))
    first, second = neighbors[cf, 0] == owner, neighbors[cf, 1] == owner
    assert np.all(first ^ second)
    signs = first.astype('int64')*2-1
    with h5py.File(output, 'x') as f:
        def put(name, value):
            f.create_dataset(name, data=np.asarray(value).reshape(-1))

        def text(name, value):
            put(name, np.frombuffer(value.encode('utf-8'), dtype='uint8'))

        def attribute(name, location, unit, value):
            value = np.asarray(value)
            p = 'arrays/' + name + '/'
            put(p+'values', value)
            put(p+'components', np.array([value.shape[1] if value.ndim == 2 else 1], dtype='uint64'))
            text(p+'location', location)
            text(p+'unit', unit)

        f.create_group('arrays')
        put('schema_version', np.array([1], dtype='uint64'))
        put('dimension', np.array([dim], dtype='uint64'))
        text('length_unit', 'm'); text('z_convention', 'depth')
        put('points', coords.astype('float64'))
        for name, size in [('node_ids', len(coords)), ('face_ids', nf), ('cell_ids', nc)]:
            put(name, np.arange(size, dtype='uint64'))
        for name, data in [('face_offsets', fo), ('face_nodes', fn), ('cell_offsets', co), ('cell_faces', cf), ('signs', signs)]:
            put(name, data)
        if cell_faces.ndim == 2:
            attribute('cell_face_tags', 'incidence', '1', cell_faces[:, 1].astype('int64'))
        specs = [
            ('nodes/global', 'mrst_node_global', 'node', '1', 'uint64'),
            ('faces/global', 'mrst_face_global', 'face', '1', 'uint64'),
            ('cells/global', 'mrst_cell_global', 'cell', '1', 'uint64'),
            ('cells/indexMap', 'mrst_index_map', 'cell', '1', 'uint64'),
            ('cartDims', 'cart_dims', 'metadata', '1', 'uint64'),
            ('faces/tag', 'face_tags', 'face', '1', 'int64'),
            ('faces/areas', 'reference_face_areas', 'face', 'm' if dim == 2 else 'm2', 'float64'),
            ('faces/normals', 'reference_face_normals', 'face', 'm' if dim == 2 else 'm2', 'float64'),
            ('faces/centroids', 'reference_face_centroids', 'face', 'm', 'float64'),
            ('cells/volumes', 'reference_cell_volumes', 'cell', 'm2' if dim == 2 else 'm3', 'float64'),
            ('cells/centroids', 'reference_cell_centroids', 'cell', 'm', 'float64'),
        ]
        for path, name, location, unit, dtype in specs:
            if (source/path/'data.csv').exists():
                value = read(path)
                if value.ndim == 2 and value.shape[1] == 2 and ('centroids' in name or 'normals' in name):
                    value = np.column_stack([value, np.zeros(len(value))])
                if dtype == 'uint64': value = integers(value)
                if name.endswith('_global') and value.shape[0] != {'node':len(coords),'face':nf,'cell':nc}[location]:
                    location = 'metadata'
                attribute(name, location, unit, value.astype(dtype))
        source_files = [p.relative_to(source).as_posix() for p in source.rglob('data.csv')]
        ignored = [p for p in source_files if p[:-9] not in consumed]
        text('provenance', json.dumps({'source':str(source), 'sha256':consumed,
             'unmigrated_source_files':ignored, 'producer':'mrst_folder_audit.py'}, ensure_ascii=False))
    return {'cells':nc, 'faces':nf, 'nodes':len(coords), 'source_hashes':consumed,
            'unmigrated_source_files':ignored}


def measured_run(command, folder):
    started = time.monotonic()
    log = folder/'convert.log'
    peak = 0
    with log.open('w', encoding='utf-8') as stream:
        process = subprocess.Popen(command, cwd=folder, stdout=stream, stderr=subprocess.STDOUT)
        monitored = psutil.Process(process.pid)
        while process.poll() is None:
            try: peak = max(peak, monitored.memory_info().rss)
            except psutil.NoSuchProcess: pass
            time.sleep(0.1)
        code = process.wait()
    return {'exit_code':code, 'seconds':time.monotonic()-started, 'peak_rss_bytes':peak,
            'log':log.read_text(encoding='utf-8', errors='replace')}


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--study-root', type=Path, required=True)
    parser.add_argument('--converter', type=Path, required=True)
    parser.add_argument('--output-dir', type=Path, required=True)
    parser.add_argument('--maximum-cells', type=int, default=3000000)
    args = parser.parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)
    results = []
    sources = sorted(p.parents[2] for p in args.study_root.rglob('G/cells/num/data.csv'))
    for i, source in enumerate(sources):
        count = int(np.loadtxt(source/'cells/num/data.csv'))
        if count > args.maximum_cells: continue
        folder = args.output_dir/f'grid-{i:02d}'
        folder.mkdir()
        item = {'source':str(source)}
        print('START', i, count, str(source), flush=True)
        try:
            item.update(migrate(source, folder/'input.h5'))
            item.update(measured_run([str(args.converter.resolve()), 'input.h5', 'output.h5', '--geometry', '--graph'], folder))
            if (folder/'output.h5').exists():
                item['equal_hdf5_datasets'] = compare_hdf5(folder/'input.h5', folder/'output.h5')
        except Exception as error:
            item['error'] = repr(error)
        results.append(item)
        (args.output_dir/'results.json').write_text(json.dumps(results, ensure_ascii=False, indent=2), encoding='utf-8')
        print('DONE', i, json.dumps({k:v for k,v in item.items() if k not in ['source_hashes','source']}, ensure_ascii=False), flush=True)


if __name__ == '__main__':
    main()
