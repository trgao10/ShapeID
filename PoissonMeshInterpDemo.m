%% preparation
close all;
path(pathdef);
addpath(path,genpath([pwd '/utils/']));

%% setup parameter
Names = {'j07', 'j08', 'j09'};
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

%% extract affine transformations
domainMesh = Mesh('off', [reparametrized_path 'domainMesh.off']);

AffineTransformations = zeros(GroupSize,3,3,MeshList{j}.nF);
for j=1:GroupSize
    for k=1:MeshList{j}.nF
        
        AffineTransformations(j,:,:,k) = ;
    end
end




