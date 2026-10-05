function [G, rock, metadata] = import_mpmc_mesh(filename)
%IMPORT_MPMC_MESH Reconstruct an MRST grid from MPMC HDF5 schema v1.
%   [G, ROCK, META] = IMPORT_MPMC_MESH(FILE) restores topology, references,
%   indexMap and typed attributes. Pass META to EXPORT_MPMC_MESH on rewrite.
%   Transport does not certify computational geometry. Run computeGeometry
%   explicitly when desired and compare against the transported references.
%   See also export_mpmc_mesh, h5read.
root=h5info(filename,'/');
expected={'cell_faces','cell_ids','cell_offsets','dimension','face_ids','face_nodes', ...
    'face_offsets','length_unit','node_ids','points','provenance','schema_version','signs','z_convention'};
assert(isequal(sort({root.Datasets.Name}),sort(expected)) && ...
    isequal({root.Groups.Name},{'/arrays'}) && isempty(root.Links), 'mpmc:Schema','Unknown/missing schema entries.');
assert(isequal(read('/schema_version'),uint64(1)), 'mpmc:Version','Unsupported schema.');
rawDimension=read('/dimension');
assert(isa(rawDimension,'uint64'), 'mpmc:Dimension','Dimension must be UInt64.');
dim = double(rawDimension);
assert(isscalar(dim) && any(dim==[2 3]), 'mpmc:Dimension','Invalid dimension.');
assert(strcmp(txt('/length_unit'),'m'), 'mpmc:Unit','Coordinates must be metres.');
xyz = reshape(read('/points'),3,[])';
assert(isa(xyz,'double') && all(isfinite(xyz(:))), 'mpmc:Coordinates','Invalid coordinates.');
if dim==2, assert(all(xyz(:,3)==0), 'mpmc:Surface','Embedded surface unsupported.'); end
metadata.node_ids = read('/node_ids'); metadata.face_ids = read('/face_ids'); metadata.cell_ids = read('/cell_ids');
np = numel(metadata.node_ids); nf = numel(metadata.face_ids); nc = numel(metadata.cell_ids);
assert(size(xyz,1)==np, 'mpmc:Shape','Coordinate count mismatch.');
for name = {'node_ids','face_ids','cell_ids'}
    ids = metadata.(name{1});
    assert(isa(ids,'uint64') && numel(unique(ids))==numel(ids), 'mpmc:ID','Invalid stable IDs.');
end
fo=refs('/face_offsets'); fn=refs('/face_nodes'); co=refs('/cell_offsets'); cf=refs('/cell_faces');
assert(numel(fo)==nf+1 && fo(1)==1 && fo(end)==numel(fn)+1 && all(diff(fo)>=0) && all(fn<=np), 'mpmc:CSR','Invalid face connectivity.');
assert(numel(co)==nc+1 && co(1)==1 && co(end)==numel(cf)+1 && all(diff(co)>=0) && all(cf<=nf), 'mpmc:CSR','Invalid cell connectivity.');
signs=read('/signs');
assert(isa(signs,'int64') && numel(signs)==numel(cf) && all(abs(signs)==1), 'mpmc:Signs','Invalid incidence signs.');
owner=repelem((1:nc)',diff(co)); side=1+(signs<0);
linear=cf+(double(side)-1)*nf;
assert(numel(unique(linear))==numel(linear), 'mpmc:Incidence','Duplicate face owner.');
neighbors=zeros(nf,2); neighbors(linear)=owner;
assert(all(any(neighbors>0,2)), 'mpmc:Incidence','Orphan face.');
G=struct('griddim',dim,'type',{{'import_mpmc_mesh'}});
G.nodes=struct('num',np,'coords',xyz(:,1:dim));
G.faces=struct('num',nf,'nodePos',fo,'nodes',fn,'neighbors',neighbors);
G.cells=struct('num',nc,'facePos',co,'faces',cf);
rock=struct(); metadata.arrays=struct();
metadata.z_convention=txt('/z_convention'); metadata.provenance=txt('/provenance');
assert(any(strcmp(metadata.z_convention,{'depth','elevation'})), 'mpmc:Z','Unknown Z convention.');
info=h5info(filename,'/arrays');
assert(isempty(info.Datasets) && isempty(info.Links), 'mpmc:Schema','Unexpected arrays entry.');
for i=1:numel(info.Groups)
    entry=info.Groups(i); p=entry.Name; parts=strsplit(p,'/'); name=parts{end};
    assert(isequal(sort({entry.Datasets.Name}),{'components','location','unit','values'}) && ...
        isempty(entry.Groups) && isempty(entry.Links), 'mpmc:Schema','Unexpected attribute entry.');
    raw=read([p '/components']); assert(isa(raw,'uint64') && all(raw<=flintmax), 'mpmc:Shape','Invalid component type.');
    n=double(raw); value=read([p '/values']);
    assert(isscalar(n) && n>0 && mod(numel(value),n)==0, 'mpmc:Shape','Invalid attribute shape.');
    value=reshape(value,n,[])';
    a=struct('location',txt([p '/location']),'unit',txt([p '/unit']),'values',value);
    metadata.arrays.(name)=a;
end
validate_mpmc_arrays(metadata.arrays,dim,np,nf,nc,numel(cf));
for field=fieldnames(metadata.arrays)'
    name=field{1}; value=metadata.arrays.(name).values;
    if startsWith(name,'rock_'), rock.(name(6:end))=value;
    elseif strcmp(name,'mrst_index_map'), G.cells.indexMap=exact_double(value);
    elseif strcmp(name,'mrst_node_global'), G.nodes.global=exact_double(value);
    elseif strcmp(name,'mrst_face_global'), G.faces.global=exact_double(value);
    elseif strcmp(name,'mrst_cell_global'), G.cells.global=exact_double(value);
    elseif strcmp(name,'cart_dims'), G.cartDims=exact_double(value(:)');
    elseif strcmp(name,'cell_face_tags'), G.cells.faces(:,2)=exact_double(value);
    elseif strcmp(name,'face_tags'), G.faces.tag=value;
    elseif strcmp(name,'nnc_cells'), G.nnc.cells=exact_double(value)+1;
    elseif startsWith(name,'nnc_'), G.nnc.(name(5:end))=value;
    elseif any(strcmp(name,{'reference_face_areas','reference_face_normals','reference_face_centroids', ...
            'reference_cell_volumes','reference_cell_centroids'}))
        key=name(11:end); pos=strfind(key,'_'); kind=key(1:pos(1)-1); field=key(pos(1)+1:end);
        if strcmp(kind,'face'), kind='faces'; else, kind='cells'; end
        if size(value,2)==3, value=value(:,1:dim); end
        G.(kind).(field)=value;
    end
end
    function v=read(p)
        descriptor=h5info(filename,p);
        assert(numel(descriptor.Dataspace.Size)==1, 'mpmc:Shape','Dataset must be rank one.');
        v=h5read(filename,p); v=v(:);
    end
    function v=txt(p)
        v=read(p); assert(isa(v,'uint8'), 'mpmc:Text','Text must be UTF8 bytes.');
        v=native2unicode(v','UTF-8');
    end
    function v=exact_double(v)
        assert(all(abs(double(v(:)))<=flintmax), 'mpmc:Precision','MRST double indexing would lose precision.');
        v=double(v);
    end
    function v=refs(p)
        v=read(p);
        assert(isa(v,'uint64') && all(v<uint64(flintmax)), 'mpmc:Index','Invalid/lossy local reference.');
        v=double(v)+1;
    end
end
