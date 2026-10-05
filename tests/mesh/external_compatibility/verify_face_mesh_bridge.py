"""Independent VTK -> C++ -> HDF5 -> C++ -> VTK/HDF5 contract tests.

Synthetic analytical prisms/polygons, not reservoir physics validation. Checks
typed IDs, shared faces, geometry, NNC graph, every HDF5 array, and rejection of
contradictory/unsafe inputs. Requires pinned vtk, numpy, h5py from requirements.
"""
import argparse
import hashlib
import json
import platform
import shutil
import subprocess
from pathlib import Path
import xml.etree.ElementTree as ET

import h5py
import numpy as np
import vtk
from vtk.util.numpy_support import numpy_to_vtk, vtk_to_numpy
from mrst_folder_audit import compare_hdf5


def produce(path, dimension):
    xy = np.array([[0,0],[1,0],[1,1],[.5,1.5],[0,1]], dtype=float)
    if dimension == 3:
        points = np.vstack([np.column_stack([xy,np.zeros(5)]), np.column_stack([xy,np.ones(5)]),
                            [[2,0,0],[2,1,0],[2,0,1],[2,1,1]]])
    else:
        points = np.vstack([np.column_stack([xy,np.zeros(5)]), [[2,0,0],[2,1,0]]])
    grid = vtk.vtkUnstructuredGrid(); pts = vtk.vtkPoints()
    pts.SetData(numpy_to_vtk(points, deep=True)); grid.SetPoints(pts)
    if dimension == 3:
        faces = [list(reversed(range(5))), list(range(5,10))]
        faces += [[i,(i+1)%5,(i+1)%5+5,i+5] for i in range(5)]
        stream = vtk.vtkIdList(); stream.InsertNextId(len(faces))
        for face in faces:
            stream.InsertNextId(len(face))
            for node in face: stream.InsertNextId(node)
        grid.InsertNextCell(vtk.VTK_POLYHEDRON, stream)
        grid.InsertNextCell(vtk.VTK_HEXAHEDRON, 8, [1,10,11,2,6,12,13,7])
    else:
        grid.InsertNextCell(vtk.VTK_POLYGON, 5, [0,1,2,3,4])
        grid.InsertNextCell(vtk.VTK_QUAD, 4, [1,5,6,2])
    node_ids = np.arange(len(points),dtype='uint64')+np.uint64(2**53+101)
    node_ids[-1] = np.iinfo('uint64').max
    for data,name,values in [(grid.GetPointData(),'mpmc_global_vertex_id',node_ids),
                             (grid.GetCellData(),'mpmc_global_cell_id',np.array([2**53+7,99],dtype='uint64')),
                             (grid.GetCellData(),'region',np.array([-3,8],dtype='int64'))]:
        a=numpy_to_vtk(values,deep=True); a.SetName(name); data.AddArray(a)
    writer=vtk.vtkXMLUnstructuredGridWriter(); writer.SetInputData(grid)
    writer.SetFileName(str(path)); writer.SetDataModeToAscii(); writer.SetCompressor(None)
    assert writer.Write()==1
    return grid, points


def run(exe, directory, source, target, dimension, expected=0, extra=()):
    command=[str(exe),source,target,'--geometry','--graph']
    if dimension==2: command+=['--2d']
    command+=list(extra)
    p=subprocess.run(command,cwd=directory,capture_output=True,text=True,timeout=120)
    assert p.returncode==expected, (command,p.returncode,p.stdout,p.stderr)
    return p.stdout+p.stderr


def read_vtk(path):
    reader=vtk.vtkXMLUnstructuredGridReader(); reader.SetFileName(str(path)); reader.Update()
    result=reader.GetOutput(); assert result.GetNumberOfCells()>0
    return result


def signature(grid, cell):
    c=grid.GetCell(cell)
    if c.GetCellDimension()==2:
        return {frozenset([c.GetPointId(i),c.GetPointId((i+1)%c.GetNumberOfPoints())]) for i in range(c.GetNumberOfPoints())}
    return {frozenset(c.GetFace(j).GetPointId(i) for i in range(c.GetFace(j).GetNumberOfPoints())) for j in range(c.GetNumberOfFaces())}


def main():
    parser=argparse.ArgumentParser(); parser.add_argument('--converter',type=Path,required=True)
    parser.add_argument('--output-dir',type=Path,required=True); args=parser.parse_args()
    exe=args.converter.resolve(); args.output_dir.mkdir(parents=True,exist_ok=True)
    repo=Path(__file__).resolve().parents[3]
    evidence={'vtk':vtk.vtkVersion.GetVTKVersion(),'h5py':h5py.__version__,
              'platform':platform.platform(),'python':platform.python_version(),
              'converter_sha256':hashlib.sha256(exe.read_bytes()).hexdigest(),
              'oracle_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
              'head':subprocess.check_output(['git','rev-parse','HEAD'],cwd=repo,text=True).strip(),
              'dirty':bool(subprocess.check_output(['git','status','--porcelain'],cwd=repo)),
              'cases':[]}
    for dim in [2,3]:
        directory=args.output_dir/f'{dim}d'; directory.mkdir()
        original,points=produce(directory/'external.vtu',dim)
        log=run(exe,directory,'external.vtu','mesh.h5',dim)
        assert 'volume=2.25' in log or abs(float(log.split('volume=')[1].split()[0])-2.25)<1e-12
        log+=run(exe,directory,'mesh.h5','export.vtu',dim)
        returned=read_vtk(directory/'export.vtu')
        assert np.array_equal(vtk_to_numpy(returned.GetPoints().GetData()),points)
        for c in range(2): assert signature(returned,c)==signature(original,c)
        shared=signature(returned,0)&signature(returned,1); assert len(shared)==1
        for location,name in [('GetPointData','mpmc_global_vertex_id'),('GetCellData','mpmc_global_cell_id'),('GetCellData','region')]:
            assert np.array_equal(vtk_to_numpy(getattr(returned,location)().GetArray(name)),vtk_to_numpy(getattr(original,location)().GetArray(name)))
        size=vtk.vtkCellSizeFilter(); size.SetInputData(returned);size.Update()
        measure=vtk_to_numpy(size.GetOutput().GetCellData().GetArray('Area' if dim==2 else 'Volume'))
        np.testing.assert_allclose(measure,[1.25,1],rtol=1e-13,atol=1e-13)
        log+=run(exe,directory,'export.vtu','returned.h5',dim)
        equal=compare_hdf5(directory/'mesh.h5',directory/'returned.h5')
        # VTK must also be able to rewrite all schema arrays without changing
        # the contract; metadata does not rely on custom XML attributes.
        writer=vtk.vtkXMLUnstructuredGridWriter(); writer.SetInputData(returned); writer.SetFileName(str(directory/'vtk-rewrite.vtu'))
        writer.SetDataModeToAscii();writer.SetCompressor(None);assert writer.Write()==1
        log+=run(exe,directory,'vtk-rewrite.vtu','vtk-rewrite.h5',dim)
        compare_hdf5(directory/'mesh.h5',directory/'vtk-rewrite.h5')
        # No overwrite; bad schema/version/type/sign/extent/reference geometry.
        negatives=0
        for name,mutate in [
            ('version',lambda h:h['schema_version'].__setitem__(0,2)),
            ('offset',lambda h:h['face_offsets'].__setitem__(1,2**63)),
            ('sign',lambda h:h['signs'].__setitem__(0,0)),
            ('id',lambda h:h['node_ids'].__setitem__(1,h['node_ids'][0])),
            ('unknown',lambda h:h.create_dataset('uncontracted',data=[1])),
        ]:
            file=directory/(name+'.h5');shutil.copyfile(directory/'mesh.h5',file)
            with h5py.File(file,'r+') as h: mutate(h)
            run(exe,directory,file.name,name+'-out.h5',dim,expected=1);negatives+=1
        file=directory/'bad-reference.h5';shutil.copyfile(directory/'mesh.h5',file)
        with h5py.File(file,'r+') as h:
            g=h.create_group('arrays/reference_cell_volumes')
            g['values']=np.array([999.,999.]);g['components']=np.array([1],dtype='uint64')
            for name,value in [('location','cell'),('unit','m2' if dim==2 else 'm3')]:g[name]=np.frombuffer(value.encode(),dtype='uint8')
        run(exe,directory,file.name,'bad-reference-out.h5',dim,expected=2);negatives+=1
        file=directory/'bad-reference-association.h5';shutil.copyfile(directory/'mesh.h5',file)
        with h5py.File(file,'r+') as h:
            g=h.create_group('arrays/reference_cell_volumes')
            g['values']=np.array([1.25,1.]);g['components']=np.array([1],dtype='uint64')
            for name,value in [('location','metadata'),('unit','m2' if dim==2 else 'm3')]:
                g[name]=np.frombuffer(value.encode(),dtype='uint8')
        run(exe,directory,file.name,'bad-reference-association-out.h5',dim,expected=1);negatives+=1
        # Contradictory canonical incidence must not override the standard cells.
        tree=ET.parse(directory/'export.vtu');a=tree.find(".//FieldData/DataArray[@Name='mpmc_fm_signs']")
        tokens=a.text.split();tokens[0]=str(-int(tokens[0]));a.text=' '.join(tokens);tree.write(directory/'bad.vtu')
        run(exe,directory,'bad.vtu','bad-out.h5',dim,expected=1);negatives+=1
        (directory/'validation.log').write_text(log,encoding='utf-8')
        evidence['cases'].append({'dimension':dim,'cells':2,'shared_faces':1,'equal_datasets':equal,'negative_controls':negatives})
    evidence['files_sha256']={p.relative_to(args.output_dir).as_posix():hashlib.sha256(p.read_bytes()).hexdigest()
                              for p in args.output_dir.rglob('*') if p.is_file()}
    (args.output_dir/'evidence.json').write_text(json.dumps(evidence,indent=2),encoding='utf-8')
    print('[PASS] independent.face_mesh.hdf5_vtu dimensions=2 negative_controls=16')


if __name__=='__main__':main()
