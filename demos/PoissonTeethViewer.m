%% preparation
close all;
path(pathdef);
addpath(path,genpath([pwd '/utils/']));

%% .off mesh path
GroupName = 'Galago';
PoissonTeethPath = ['./meshes/PoissonTeeth/' GroupName '/off/'];

%% read in all meshes
numRandMeshes = 10;

MeshList = cell(numRandMeshes,1);

cback = 0;
for j=11:20
    MeshList{j} = Mesh('off', [PoissonTeethPath GroupName '_' sprintf('%02d', j) '.off']);
    for cc=1:cback
        fprintf('\b');
    end
    cback = fprintf(['%4d/' num2str(numRandMeshes) ' done.\n'], j);
end

drawMeshList(MeshList, struct('DisplayLayout',[2,5],'linkCamera','on','DisplayOrient','Horizontal'));
