function validate_mpmc_arrays(arrays, dim, np, nf, nc, ni)
%VALIDATE_MPMC_ARRAYS Validate typed attributes of the version 1 bridge.
nn = 0;
if isfield(arrays,'nnc_cells')
    a=arrays.nnc_cells; v=a.values;
    assert(strcmp(a.location,'nnc') && isa(v,'uint64') && size(v,2)==2 && ...
        all(v(:)<nc) && all(v(:,1)~=v(:,2)), 'mpmc:NNC','Invalid NNC endpoints.');
    nn=size(v,1);
end
expected=struct('node',np,'face',nf,'cell',nc,'incidence',ni,'nnc',nn);
for field=fieldnames(arrays)'
    name=field{1}; a=arrays.(name); v=a.values;
    assert(~isempty(regexp(name,'^[a-zA-Z0-9_]+$','once')), 'mpmc:Name','Invalid array name.');
    assert(ismatrix(v) && size(v,2)>0 && isreal(v) && all(isfinite(v(:))) && ...
        any(strcmp(class(v),{'double','int64','uint64'})), 'mpmc:Array','Invalid typed array.');
    assert(ischar(a.unit) && ~isempty(a.unit), 'mpmc:Unit','Attribute unit is required.');
    if isfield(expected,a.location)
        assert(size(v,1)==expected.(a.location), 'mpmc:Shape','Attribute row count mismatch.');
    else
        assert(strcmp(a.location,'metadata'), 'mpmc:Location','Unknown attribute location.');
    end
    [reserved,location,components,unit]=reference_geometry_contract(name,dim);
    if reserved
        assert(strcmp(a.location,location) && isa(v,'double') && ...
            size(v,2)==components && strcmp(a.unit,unit), ...
            'mpmc:ReferenceGeometry','Invalid reserved reference-geometry array contract.');
    end
end
if isfield(arrays,'mrst_index_map')
    a=arrays.mrst_index_map; v=a.values;
    assert(strcmp(a.location,'cell') && isa(v,'uint64') && size(v,2)==1 && ...
        all(v>=1) && numel(unique(v))==nc, 'mpmc:IndexMap','Invalid active-to-logical mapping.');
end
if isfield(arrays,'cart_dims')
    a=arrays.cart_dims; v=a.values;
    assert(strcmp(a.location,'metadata') && isa(v,'uint64') && ...
        isequal(size(v),[dim 1]) && all(v>0), 'mpmc:CartDims','Invalid logical dimensions.');
    count=uint64(1);
    for k=1:dim
        assert(v(k)<=idivide(intmax('uint64'),count), 'mpmc:CartDims','Logical size overflow.');
        count=count*v(k);
    end
    if isfield(arrays,'mrst_index_map')
        assert(all(arrays.mrst_index_map.values<=count), 'mpmc:IndexMap','indexMap exceeds logical grid.');
    end
end
end

function [reserved,location,components,unit]=reference_geometry_contract(name,dim)
reserved=true; location=''; components=0; unit='';
if dim==2
    face_measure='m'; cell_measure='m2';
else
    face_measure='m2'; cell_measure='m3';
end
switch name
    case 'reference_face_areas'
        location='face'; components=1; unit=face_measure;
    case 'reference_face_normals'
        location='face'; components=3; unit=face_measure;
    case 'reference_face_centroids'
        location='face'; components=3; unit='m';
    case 'reference_cell_volumes'
        location='cell'; components=1; unit=cell_measure;
    case 'reference_cell_centroids'
        location='cell'; components=3; unit='m';
    otherwise
        reserved=false;
end
end
