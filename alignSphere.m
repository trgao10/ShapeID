%% preparation
close all;
path(pathdef);
addpath(path,genpath([pwd '/utils/']));
mesh_path = '~/Downloads/SphereTeeth/';

Names = {'a15', 'a16', 'b19', 'h08', 'H10', 'H19'};

%%% read in meshes and observer landmarks
MeshList = cell(length(Names), 1);
for j=1:length(MeshList)
    MeshList{j} = Mesh('off', [mesh_path Names{j} '_sphere.off']);
    MeshList{j}.V = MeshList{j}.V - repmat(mean(MeshList{j}.V, 2), 1, MeshList{j}.nV);
    
    fid = fopen([mesh_path Names{j} '_landmarks.txt']);
    
    tline = fgets(fid);
    while ischar(tline) && isempty(strfind(tline, 'vertex'))
        tline = fgets(fid);
    end
    MeshList{j}.Aux.ObLmkInd = fscanf(fid, '%i');    
    fclose(fid);
end

%%% align every shape to the first shape
for j=2:length(MeshList)
    [U,~,V] = svd(MeshList{j}.V(:,MeshList{j}.Aux.ObLmkInd)*MeshList{1}.V(:,MeshList{1}.Aux.ObLmkInd)');
    R = V*U';
    MeshList{j}.V = R * MeshList{j}.V;
    if (det(R) < 0)
        MeshList{j}.F = MeshList{j}.F([1, 3, 2], :);
    end
end

for j=1:length(MeshList)
    MeshList{j}.Write([Names{j} '_sphere_aligned.off'], 'off', []);
end

