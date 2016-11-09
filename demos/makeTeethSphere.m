clear vas
close all
path(pathdef);
addpath(path,genpath([pwd '/utils/']));

meshPath = '~/Work/MATLAB/DATA/PNAS/meshes/';
% meshPath = '~/Work/MATLAB/DATA/HDM/meshes/';
% meshPath = '~/Downloads/TheSix/';
% outputPath = '~/Dropbox/sphTeethMinSurf/';

dir_struct = dir(meshPath);
names = cell(length(dir_struct)-2,1);

for j=3:length(dir_struct)
    meshName = strtok(dir_struct(j).name, '.');
    names{j-2} = meshName;
end

for j = 1:length(names)
    outputMeshName = [outputPath names{j} '_close.off'];
    if ~exist(outputMeshName, 'file')
        G = Mesh('off', [meshPath names{j} '.off']);
        sTooth = G.MakeToothSphere(1500, 4000, 0.98);
        sTooth.Write(outputMeshName, 'off');
    end
end
