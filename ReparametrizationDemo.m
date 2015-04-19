%% preparation
% close all;
path(pathdef);
addpath(path,genpath([pwd '/utils/']));

%% setup parameter
Names = {'j07', 'j08', 'j09'};
GroupSize = length(Names);

%% setup paths
base_path = [pwd '/'];
data_path = '../DATA/PNAS/';
spreadsheet_path = [data_path 'ClassificationTable.xlsx'];
sample_path = '../cPdist/samples/PNAS/';
LandmarksPath = [data_path 'landmarks_teeth.mat'];
MeshesPath = [data_path 'meshes/'];
MeshSuffix = '_sas.off';

%% load taxa codes
taxa_file = [data_path 'teeth_taxa_table.mat'];
taxa_code = load(taxa_file);
taxa_code = taxa_code.taxa_code;

%% load meshes
MeshList = cell(GroupSize,1);
for j=1:GroupSize
    load([sample_path taxa_code{strcmpi(taxa_code, Names{j})} '.mat']);
    G.Aux.ObLmk = GetLandmarks(G.Aux.name, LandmarksPath, [MeshesPath G.Aux.name MeshSuffix]);
    MeshList{j} = G;
end

originalMeshes = drawMeshList(MeshList, []);
set(originalMeshes, 'Name', 'Orignal Meshes');

%% reparametrize all meshes

%%%% Step 1: Get all conformal parametrizations
R = cell(GroupSize,1);
R{1} = eye(3);

FlatMeshList = cell(GroupSize,1);
FlatMeshOrientation = zeros(GroupSize,1); %% '1' means parametrization flipped
FlatMeshList{1} = Mesh('VF', MeshList{1}.Aux.UniformizationV, MeshList{1}.F);
for j=2:GroupSize
    FlatVertices = MeshList{j}.Aux.UniformizationV;

    %%% Procrustes matching observer landmarks
    [U,~,V] = svd(FlatVertices(:,MeshList{j}.Aux.ObLmk)*FlatMeshList{1}.V(:,MeshList{1}.Aux.ObLmk)');
    R{j} = V*U';
    
    if (det(R{j}(1:2,1:2)) < 0)
        FlatMeshOrientation(j) = 1;
        disp(['Mesh ' num2str(j) ' flipped!']);
    end

    FlatVertices = R{j}*FlatVertices;
    FlatMeshList{j} = Mesh('VF', FlatVertices, MeshList{j}.F);
end

%%% draw Procrustes aligned meshes
drawMeshList(FlatMeshList, []);

%%%% Step 2: TPS 2~end meshes to the first mesh
DeformedMeshList = cell(GroupSize,1);
DeformedMeshList{1} = FlatMeshList{1};
TPS_FEATURESN = DISCtoPLANE(FlatMeshList{1}.V(1:2, MeshList{1}.Aux.ObLmk)','d2p');
for j=2:GroupSize
    TPS_FEATURESM = DISCtoPLANE(FlatMeshList{j}.V(1:2,MeshList{j}.Aux.ObLmk)','d2p');
    tP = DISCtoPLANE(FlatMeshList{j}.V(1:2,:)','d2p');
    [ftps] = TEETH_calc_tps(TPS_FEATURESM,TPS_FEATURESN-TPS_FEATURESM);
    pt = tP + TEETH_eval_tps(ftps,tP);
    DeformedMeshList{j} = Mesh('VF', DISCtoPLANE(pt,'p2d')', FlatMeshList{j}.F);
    DeformedMeshList{j}.V(:,isnan(compl(DeformedMeshList{j}.V))) = 1;
end

%%% draw TPS deformed meshes
drawMeshList(DeformedMeshList, []);

%%%% Step 3: Get standard unit disk mesh
%%%%% generate mesh from distmesh
n = 300;
[X,Y] = meshgrid(-1:(1/n):1);
V = griddata(DeformedMeshList{1}.V(1,:),DeformedMeshList{1}.V(2,:),MeshList{1}.Aux.Conf,X,Y);
% V = 1./V;
V(isnan(V))=0;
fd = @(p) sqrt(sum(p.^2,2))-0.95^2; %%% choose a smaller disk
[UDVs, UDFs] = distmesh2d(fd, @huniform, 0.03, [-1,-1;1,1], MeshList{1}.Aux.UniformizationV(1:2,MeshList{1}.Aux.ObLmk)');
% [UDVs, UDFs] = distmesh2d(fd, @hmatrix, 0.03, [-1,-1;1,1], [], X, Y, V);
% [UDVs, UDFs] = distmesh2d(fd, @hmatrix, 0.03, [-1,-1;1,1],...
%     MeshList{1}.Aux.UniformizationV(1:2,MeshList{1}.Aux.ObLmk)', X, Y, V);
% [UDVs, UDFs] = distmesh2d(fd, @huniform, 0.03, [-1,-1;1,1], []);
%%%%% generate domain mesh from first mesh
% UDVs = DeformedMeshList{1}.V(1:2,:)';
% UDVs(sum(UDVs.^2,2)>0.95^2,:) = [];
% DT = delaunayTriangulation(UDVs);
% UDFs = DT.ConnectivityList;
UnitDistMesh = Mesh('VF', UDVs', UDFs');
figure;UnitDistMesh.draw();

%%%% Step 4: Re-parametrize all meshes
% BBox = [-1.1, -1.1, 1.1, 1.1; -1.1, 1.1, -1.1, 1.1];
ReparametrizedMeshList = cell(GroupSize,1);
for j=1:GroupSize
    disp(['Barycentric Interpolation for Mesh ' num2str(j) '...']);
    NewMeshVertices = zeros(3, size(UDVs,1));
    
    TR = triangulation(DeformedMeshList{j}.F', DeformedMeshList{j}.V(1:2,:)');
    
    cback = 0;
    for k=1:size(UDVs,1)
        BC = TR.cartesianToBarycentric((1:DeformedMeshList{j}.nF)',repmat(UDVs(k,:),DeformedMeshList{j}.nF,1));
        
        tind = find(all(BC>-1e-10,2));
        if(numel(tind)>1)
            tind = tind(1);
        elseif numel(tind)==0
            warning(['For j = ' num2str(j) ', Point ', num2str(k), ' was not found in any triangle.']);
            continue;
        end
        BC = BC(tind,:);
        
        NewMeshVertices(:, k) = MeshList{j}.V(:,TR.ConnectivityList(tind,:))*BC';
        
        for cc=1:cback
            fprintf('\b');
        end
        cback = fprintf(['%4d/' num2str(size(UDVs,1)) ' done.\n'], k);
    end
    
    ReparametrizedMeshList{j} = Mesh('VF', NewMeshVertices, UDFs');
end

%%% draw reparametrized meshes
reparametrizedMeshes = drawMeshList(ReparametrizedMeshList, []);
set(reparametrizedMeshes, 'Name', 'Reparametrized Meshes');

%% write reparametrized meshes to .off files
for j=1:GroupSize
    ReparametrizedMeshList{j}.Write(['./meshes/reparametrized/' MeshList{j}.Aux.name '.off'], 'off', []);
end

domainMesh = DeformedMeshList{1};
domainMesh.Write('./meshes/reparametrized/domainMesh.off', 'off', []);


