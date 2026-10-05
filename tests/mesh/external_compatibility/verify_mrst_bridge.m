function verify_mrst_bridge(outputDirectory)
%VERIFY_MRST_BRIDGE Local MRST/MATLAB producer and reconstruction checks.
%   Requires MRST startup and modules/mesh/matlab on the path. Output must be
%   a new directory. The C++/VTK test consumes these files in the local audit.
assert(~isfolder(outputDirectory),'Use a fresh evidence directory.');
mkdir(outputDirectory);
grids={computeGeometry(cartGrid([4 3],[4 3])), ...
    computeGeometry(cartGrid([3 2 2],[3 2 2])), ...
    computeGeometry(pebi(triangleGrid([0 0;1 0;1 1;0 1;.3 .4;.7 .6])))};
for k=1:numel(grids)
    G=grids{k}; rock=struct('poro',repmat(.2,G.cells.num,1), ...
        'perm',repmat([1e-12 2e-12 3e-12],G.cells.num,1));
    meta=struct();
    meta.node_ids=uint64((1:G.nodes.num)')+bitshift(uint64(1),53);
    meta.node_ids(end)=intmax('uint64');
    meta.arrays.region=struct('location','cell','unit','1','values',int64(-ones(G.cells.num,1)));
    if k==2
        G.nnc=struct('cells',[1 G.cells.num],'T',1e-10);
    else
        G.nnc=struct('cells',zeros(0,2),'T',zeros(0,1));
    end
    file=fullfile(outputDirectory,sprintf('mrst-%d.h5',k));
    export_mpmc_mesh(file,G,rock,meta);
    refused=false;
    try, export_mpmc_mesh(file,G,rock,meta); catch, refused=true; end
    assert(refused,'Existing files must not be overwritten.');
    [R,r,metadata]=import_mpmc_mesh(file);
    assert(isequal(G.nodes.coords,R.nodes.coords) && isequal(G.faces.nodes,R.faces.nodes));
    assert(isequal(G.faces.neighbors,R.faces.neighbors) && isequal(G.cells.faces,R.cells.faces));
    assert(isequal(G.nnc,R.nnc) && isequal(rock,r) && isequal(metadata.node_ids,meta.node_ids));
    recomputed=computeGeometry(R);
    assert(max(abs(recomputed.cells.volumes-G.cells.volumes))<1e-12);
    export_mpmc_mesh(fullfile(outputDirectory,sprintf('mrst-%d-return.h5',k)),R,r,metadata);
    bad=metadata; bad.arrays.region.values=bad.arrays.region.values(1:end-1);
    rejected=false; target=fullfile(outputDirectory,sprintf('invalid-%d.h5',k));
    try, export_mpmc_mesh(target,G,rock,bad); catch, rejected=true; end
    assert(rejected && ~isfile(target),'Invalid metadata must fail before writing.');
    base=metadata.arrays.reference_cell_volumes;
    bad=metadata; bad.arrays.reference_cell_volumes.location='metadata';
    target=fullfile(outputDirectory,sprintf('invalid-reference-location-%d.h5',k)); rejected=false;
    try, export_mpmc_mesh(target,G,rock,bad); catch, rejected=true; end
    assert(rejected && ~isfile(target),'Reference geometry location must fail before writing.');
    bad=metadata; bad.arrays.reference_cell_volumes.values=[base.values base.values];
    target=fullfile(outputDirectory,sprintf('invalid-reference-components-%d.h5',k)); rejected=false;
    try, export_mpmc_mesh(target,G,rock,bad); catch, rejected=true; end
    assert(rejected && ~isfile(target),'Reference geometry components must fail before writing.');
    bad=metadata; bad.arrays.reference_cell_volumes.unit='invalid';
    target=fullfile(outputDirectory,sprintf('invalid-reference-unit-%d.h5',k)); rejected=false;
    try, export_mpmc_mesh(target,G,rock,bad); catch, rejected=true; end
    assert(rejected && ~isfile(target),'Reference geometry unit must fail before writing.');
    bad=metadata; bad.arrays.reference_cell_volumes.values=int64(zeros(size(base.values)));
    target=fullfile(outputDirectory,sprintf('invalid-reference-type-%d.h5',k)); rejected=false;
    try, export_mpmc_mesh(target,G,rock,bad); catch, rejected=true; end
    assert(rejected && ~isfile(target),'Reference geometry type must fail before writing.');
end
fprintf('[PASS] MATLAB MRST bridge: 2D/3D Cartesian, irregular PEBI, typed IDs, NNC, no-overwrite, invalid fields/reference contracts\n');
end
