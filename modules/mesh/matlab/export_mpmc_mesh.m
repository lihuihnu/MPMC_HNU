function report = export_mpmc_mesh(filename, G, rock, metadata)
%EXPORT_MPMC_MESH Write an MRST ordinary FV grid to MPMC HDF5 schema v1.
%   REPORT = EXPORT_MPMC_MESH(FILE, G, ROCK, METADATA) preserves explicit
%   faces, incidence, source indexMap, typed fields and geometry references.
%   Coordinates must be metres; ROCK.perm must be square metres. No geometric
%   repair, axis flip, cell splitting, or NNC-to-face substitution is made.
%   METADATA from IMPORT_MPMC_MESH preserves stable IDs and additional arrays.
%   New files only. Validate computational geometry using the C++ mesh tool.
%   See also import_mpmc_mesh, h5create, h5write.
if nargin < 3, rock = struct(); end
if nargin < 4, metadata = struct(); end
assert(~isfile(filename), 'mpmc:Exists', 'Refusing to overwrite an existing mesh.');
assert(isfield(G, 'griddim') && any(G.griddim == [2 3]), 'mpmc:Dimension', 'Expected ordinary 2D or 3D grid.');
dim = G.griddim;
xyz = double(G.nodes.coords);
assert(size(xyz, 2) == dim && all(isfinite(xyz(:))), 'mpmc:Coordinates', 'Expected finite Cartesian coordinates in metres.');
if dim == 2, xyz(:, 3) = 0; end
np = size(xyz, 1); nf = double(G.faces.num); nc = double(G.cells.num);
assert(np == G.nodes.num, 'mpmc:Shape', 'Node count mismatch.');
fo = indices(G.faces.nodePos, true); fn = indices(G.faces.nodes, true);
co = indices(G.cells.facePos, true); cf = indices(G.cells.faces(:, 1), true);
assert(numel(fo) == nf+1 && fo(1) == 0 && fo(end) == numel(fn) && all(fo(2:end) >= fo(1:end-1)), 'mpmc:CSR', 'Invalid face CSR.');
assert(numel(co) == nc+1 && co(1) == 0 && co(end) == numel(cf) && all(co(2:end) >= co(1:end-1)), 'mpmc:CSR', 'Invalid cell CSR.');
assert(all(fn < np) && all(cf < nf), 'mpmc:Index', 'Reference out of range.');
N = double(G.faces.neighbors);
assert(isequal(size(N), [nf 2]) && all(isfinite(N(:))) && all(N(:) == fix(N(:))) && all(N(:) >= 0 & N(:) <= nc), 'mpmc:Neighbors', 'Invalid neighbors.');
assert(all(any(N > 0, 2)) && all(N(:,1) ~= N(:,2)), 'mpmc:Neighbors', 'Orphan/self-neighbor face.');
owner = repelem((1:nc)', double(diff(co)));
first = N(double(cf)+1, 1) == owner;
second = N(double(cf)+1, 2) == owner;
assert(all(xor(first, second)), 'mpmc:Incidence', 'Cell faces contradict neighbors.');
signs = int64(first)*2-1;
counts = accumarray(double(cf)+1, 1, [nf 1]);
assert(isequal(counts, sum(N > 0, 2)), 'mpmc:Incidence', 'Face incidence count mismatch.');

nodeIds = ids(metadata, 'node_ids', np);
faceIds = ids(metadata, 'face_ids', nf);
cellIds = ids(metadata, 'cell_ids', nc);
arrays = struct();
if isfield(metadata, 'arrays'), arrays = metadata.arrays; end
if isfield(G.cells, 'indexMap')
    add('mrst_index_map', 'cell', '1', indices(G.cells.indexMap, false));
end
if isfield(G, 'cartDims'), add('cart_dims', 'metadata', '1', indices(G.cartDims(:), false)); end
if size(G.cells.faces, 2) > 1
    assert(size(G.cells.faces,2) == 2, 'mpmc:Tags', 'Unsupported cell-face columns.');
    add('cell_face_tags', 'incidence', '1', integer64(G.cells.faces(:,2)));
end
for kind = {'node','face','cell'}
    singular=kind{1}; plural=[singular 's'];
    if isfield(G.(plural),'global')
        location=singular; value=G.(plural).global;
        if size(value,1)~=G.(plural).num, location='metadata'; end
        add(['mrst_' singular '_global'],location,'1',indices(value,false));
    end
end
if isfield(G.faces, 'tag'), add('face_tags', 'face', '1', integer64(G.faces.tag)); end
if isfield(G.faces, 'areas'), add('reference_face_areas', 'face', measure(dim-1), G.faces.areas); end
if isfield(G.faces, 'normals'), add('reference_face_normals', 'face', measure(dim-1), pad3(G.faces.normals)); end
if isfield(G.faces, 'centroids'), add('reference_face_centroids', 'face', 'm', pad3(G.faces.centroids)); end
if isfield(G.cells, 'volumes'), add('reference_cell_volumes', 'cell', measure(dim), G.cells.volumes); end
if isfield(G.cells, 'centroids'), add('reference_cell_centroids', 'cell', 'm', pad3(G.cells.centroids)); end
for field = fieldnames(rock)'
    name = field{1}; value = rock.(name);
    assert(isnumeric(value) || islogical(value), 'mpmc:Rock', 'Rock fields must be numeric.');
    unit = 'unspecified';
    if any(strcmp(name, {'poro','ntg','satnum','SATNUM','region'})), unit = '1'; end
    if strcmp(name, 'perm'), unit = 'm2'; end
    if islogical(value), value = uint64(value); end
    add(['rock_' name], 'cell', unit, value);
end
if isfield(G, 'nnc')
    assert(isfield(G.nnc, 'cells'), 'mpmc:NNC', 'NNC cells are required.');
    add('nnc_cells', 'nnc', '1', indices(G.nnc.cells, true));
    for field = setdiff(fieldnames(G.nnc), {'cells'})'
        name = field{1}; value = G.nnc.(name);
        assert(isnumeric(value), 'mpmc:NNC', 'Non-numeric NNC data needs an explicit contract.');
        add(['nnc_' name], 'nnc', 'unspecified', value);
    end
end
% Unknown structures are not silently discarded. Numeric application metadata
% belongs in metadata.arrays with an explicit association and unit.
check_fields(G, {'nodes','faces','cells','griddim','cartDims','type','nnc'});
check_fields(G.nodes, {'coords','num','global'});
check_fields(G.faces, {'num','nodePos','nodes','neighbors','areas','normals','centroids','tag','global'});
check_fields(G.cells, {'num','facePos','faces','indexMap','volumes','centroids','global'});
z = 'depth'; if isfield(metadata,'z_convention'), z = metadata.z_convention; end
assert(any(strcmp(z, {'depth','elevation'})), 'mpmc:Z', 'Invalid Z convention.');
source = struct('producer', 'export_mpmc_mesh', 'matlab', version, 'geometry_policy', 'source_reference');
if isfield(G,'type'), source.mrst_type = G.type; end
provenance = jsonencode(source);
if isfield(metadata,'provenance'), provenance = metadata.provenance; end

% Finish validation before creating any output file.
for field=fieldnames(arrays)'
    name=field{1}; value=arrays.(name).values;
    assert(isnumeric(value) || islogical(value), 'mpmc:Array','Expected numeric attribute.');
    if isfloat(value), value=double(value);
    elseif startsWith(class(value),'uint') || islogical(value), value=uint64(value);
    else, value=int64(value); end
    arrays.(name).values=value;
end
validate_mpmc_arrays(arrays,dim,np,nf,nc,numel(cf));

put('/schema_version', uint64(1)); put('/dimension', uint64(dim));
bytes('/length_unit', 'm'); bytes('/z_convention', z); bytes('/provenance', provenance);
put('/points', reshape(xyz', [], 1));
put('/node_ids', nodeIds); put('/face_ids', faceIds); put('/cell_ids', cellIds);
put('/face_offsets', fo); put('/face_nodes', fn);
put('/cell_offsets', co); put('/cell_faces', cf); put('/signs', signs);
% Always create the arrays group, including for a geometry-only input.
fid = H5F.open(filename, 'H5F_ACC_RDWR', 'H5P_DEFAULT');
cleanup = onCleanup(@() H5F.close(fid));
gid = H5G.create(fid, '/arrays', 'H5P_DEFAULT', 'H5P_DEFAULT', 'H5P_DEFAULT'); H5G.close(gid);
clear cleanup
for field = fieldnames(arrays)'
    name = field{1}; a = arrays.(name);
    value = a.values;
    prefix = ['/arrays/' name];
    put([prefix '/values'], reshape(value', [], 1));
    put([prefix '/components'], uint64(size(value,2)));
    bytes([prefix '/location'], a.location); bytes([prefix '/unit'], a.unit);
end
report = struct('schema_version',1,'nodes',np,'faces',nf,'cells',nc, ...
    'transport_written',true,'computational_geometry_validated',false);

    function add(name, location, unit, values)
        arrays.(name) = struct('location',location,'unit',unit,'values',values);
    end
    function bytes(path, value)
        put(path, uint8(unicode2native(char(value),'UTF-8'))');
    end
    function put(path, value)
        % A rank-one on-disk array avoids MATLAB/HDF5 dimension reversal.
        if isempty(value)
            h5create(filename,path,Inf,'Datatype',class(value),'ChunkSize',1024);
        else
            h5create(filename,path,numel(value),'Datatype',class(value));
            block = 1048576;
            for start = 1:block:numel(value)
                count = min(block,numel(value)-start+1);
                h5write(filename,path,value(start:start+count-1),start,count);
            end
        end
    end
end
function value = indices(value, subtract)
assert(isnumeric(value) && isreal(value) && all(isfinite(value(:))) && all(value(:)>=1) && all(value(:)==fix(value(:))), 'mpmc:Index','Expected positive integer indices.');
if isfloat(value), assert(all(value(:)<=flintmax), 'mpmc:Precision','Floating-point ID exceeds exact range.'); end
value = uint64(value);
if subtract, value = value-1; end
end
function v = ids(meta,name,count)
v = uint64((0:count-1)');
if isfield(meta,name), v = meta.(name); end
assert(isa(v,'uint64') && numel(v)==count && numel(unique(v))==count, 'mpmc:ID','Stable IDs must be unique UInt64.');
v = v(:);
end
function v = integer64(v)
assert(all(isfinite(v(:))) && all(v(:)==fix(v(:))) && all(abs(double(v(:)))<=flintmax), 'mpmc:Tag','Invalid integer tags.');
v = int64(v);
end
function v = pad3(v)
if size(v,2)==2, v(:,3)=0; end
assert(size(v,2)==3, 'mpmc:Shape','Expected 2/3 coordinate components.');
end
function u = measure(d)
u = 'm'; if d>1, u = ['m' num2str(d)]; end
end
function check_fields(s, allowed)
extra = setdiff(fieldnames(s), allowed);
assert(isempty(extra), 'mpmc:Unsupported','Uncontracted MRST fields: %s',strjoin(extra,', '));
end
