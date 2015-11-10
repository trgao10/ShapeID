%% preparation
close all;
path(pathdef);
addpath(path,genpath([pwd '/utils/']));

%% setup parameter
Names = {'ah16', 'Q19'};
% Weights = [1,1,1,1];
GroupSize = length(Names);

%% setup paths
base_path = [pwd '/'];
data_path = '../DATA/PNAS/';
spreadsheet_path = [data_path 'ClassificationTable.xlsx'];
reparametrized_path = './meshes/reparametrized/';
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
    G = Mesh('off', [reparametrized_path taxa_code{strcmpi(taxa_code, Names{j})} '.off']);
    G.Aux.ObLmk = GetLandmarks(taxa_code{strcmpi(taxa_code, Names{j})},...
        LandmarksPath, [MeshesPath taxa_code{strcmpi(taxa_code, Names{j})} MeshSuffix]);
    MeshList{j} = G;
end

reparametrizedMeshes = drawMeshList(MeshList, []);
set(reparametrizedMeshes, 'Name', 'Reparametrized Meshes');

domainMesh = Mesh('off', [reparametrized_path 'domainMesh.off']);

%%%% compute indices of observer landmarks on the first mesh
[~,ObCoords] = GetLandmarks(taxa_code{strcmpi(taxa_code, Names{1})},...
        LandmarksPath, [MeshesPath taxa_code{strcmpi(taxa_code, Names{1})} MeshSuffix]);
 Vtree = kdtree_build( MeshList{1}.V' );
 ObInds = kdtree_nearest_neighbor(Vtree, ObCoords);

%% reconstruction
%%% focus on checking interpolation between two meshes
WeightsFirst = linspace(0,1,9);
reconMesh = cell(size(WeightsFirst));
for j=1:length(WeightsFirst)
    if (j==1)
        [reconMesh{j},ORT,PSD] = PoissonMeshInterpolation(domainMesh, MeshList, [1-WeightsFirst(j),WeightsFirst(j)]);
    else
        reconMesh{j} = PoissonMeshInterpolation(domainMesh, MeshList, [1-WeightsFirst(j),WeightsFirst(j)], ORT, PSD);
    end
    disp([num2str(j) '/' num2str(length(WeightsFirst)) ' done.']);
end

interpolatedMeshes = drawMeshList(reconMesh, struct('DisplayLayout',[4,5],'linkCamera','on'));
set(gcf,'Name','Interpolated Meshes');

for k=1:length(WeightsFirst)
    reconMesh{k}.Write([reparametrized_path sprintf('%02d', k) '.off'], 'off', []);
end

%%% focus on checking one interpolated mesh
% reconMesh = PoissonMeshInterpolation(domainMesh, MeshList, Weights);
% figure;
% reconMesh.draw();

% %%% randomly generate 15 interpolated meshes
% reconMesh = cell(15,1);
% 
% for j=1:length(reconMesh)
%     Weights = rand(size(Names));
%     Weights = -log(Weights);
%     if (j==1)
%         [reconMesh{j},ORT,PSD] = PoissonMeshInterpolation(domainMesh, MeshList, Weights);
%     else
%         reconMesh{j} = PoissonMeshInterpolation(domainMesh, MeshList, Weights, ORT, PSD);
%     end
%     disp([num2str(j) '/' num2str(length(reconMesh)) ' done.']);
% end
% 
% interpolatedMeshes = drawMeshList(reconMesh, struct('DisplayLayout',[3,5],'linkCamera','on'));
% set(gcf,'Name','Randomly Interpolated Meshes');


