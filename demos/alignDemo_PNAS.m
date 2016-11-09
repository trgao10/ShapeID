%% preparation
close all;
path(pathdef);
addpath(path,genpath([pwd '/utils/']));

%% setup paths and parameters
% base_path = [pwd '/'];
% data_path = [base_path 'DATA/HDM/'];
data_path = '/home/trgao10/Work/MATLAB/DATA/PNAS/';
spreadsheet_path = [data_path 'ClassificationTable.xlsx'];
sample_path = [data_path 'samples/'];
taxa_file = [data_path 'teeth_taxa_table.mat'];

GroupLevel = 'Order';
GroupNames = {'Euprimates','Primates','Dermoptera','Scandentia','Incertae sedis'};

%% load taxa codes
taxa_code = load(taxa_file);
taxa_code = taxa_code.taxa_code;

%% parse GroupNames
[~,ClTable,~] = xlsread(spreadsheet_path);
Names = {};
NamesByGroup = cell(1,length(GroupNames));
TaxaByGroup = cell(1,length(GroupNames));
for j=1:length(GroupNames)
    NamesJ = ClTable(strcmpi(ClTable(1:end,strcmpi(ClTable(1,:),GroupLevel)),GroupNames{j}),1);
    Names = [Names,NamesJ{:}];
    NamesByGroup{j} = NamesJ;
    TaxaByGroup{j} = cellfun(@(x) find(strcmpi(taxa_code, x)), NamesJ);
end
GroupSize = length(Names);

%% visualization options
options.sample_path = sample_path;
options.DisplayLayout = [8,15];
options.DisplayOrient = 'Horizontal';
options.boundary = 'on';
options.names = 'off';
options.linkCamera = 'on';
options.Shading = 'Smooth';

%% load meshes
TAXAinds = zeros(GroupSize,1);
meshList = cell(1,GroupSize);
for j=1:GroupSize
    TAXAinds(j) = find(strcmpi(taxa_code,Names{j}));
    load([sample_path taxa_code{strcmpi(taxa_code,Names{j})} '.mat']);
    meshList{j} = G;
end
Names = taxa_code(TAXAinds); % match upper/lower cases

%% visualize un-aligned meshes
drawMeshList(meshList, options);
set(gcf, 'Name', 'un-aligned');

%% load rigid motions (extract from MST-based synchronization)
rigid_motions = load([data_path 'cPMSTinitRms.mat']);
R = rigid_motions.R(TAXAinds,TAXAinds);

%% align every meshe to the first one
for j=1:GroupSize
    meshList{j}.V = R{1,j} * meshList{j}.V;
    if det(R{1,j}) < 0
        meshList{j}.F = meshList{j}.F([1 3 2], :);
    end
end

%% visualize aligned meshes
drawMeshList(meshList, options);
set(gcf, 'Name', 'aligned');
