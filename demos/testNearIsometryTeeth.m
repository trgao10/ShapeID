%% preparation
close all;
path(pathdef);
addpath(path,genpath([pwd '/utils/']));

%% setup paths and parameters
% base_path = [pwd '/'];
% data_path = [base_path 'DATA/HDM/'];
data_path = '/home/trgao10/Work/MATLAB/DATA/HDM/';
spreadsheet_path = [data_path 'ClassificationTable.xlsx'];
sample_path = [data_path 'samples/'];
GroupLevel = 'Genus';
GroupNames = {'Alouatta','Ateles','Brachyteles','Callicebus','Saimiri'};

numEvals = 30;

%% load taxa codes
taxa_file = [data_path 'hdm_taxa_table.mat'];
taxa_code = load(taxa_file);
taxa_code = taxa_code.taxa_code;
GroupSize = length(taxa_code);

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

%% visualization options
options.sample_path = sample_path;
options.DisplayLayout = [5,10];
options.DisplayOrient = 'Horizontal';
options.boundary = 'on';
options.names = 'off';
options.linkCamera = 'on';
options.Shading = 'Smooth';

%% load meshes
TAXAinds = zeros(GroupSize,1);
meshList = cell(1,GroupSize);
evals = zeros(GroupSize,numEvals);
cback = 0;
for j=1:GroupSize
    TAXAinds(j) = find(strcmpi(taxa_code,Names{j}));
    load([sample_path taxa_code{strcmpi(taxa_code,Names{j})} '.mat']);
    meshList{j} = G;
    [rIdx,cIdx] = find(G.E);
    bandwidth = mean(sqrt(sum((G.V(:,rIdx)-G.V(:,cIdx)).^2)));
    DiffLap = G.ComputeDiffusionLaplacian(bandwidth);
%     CotLap = G.ComputeCotanLaplacian();
    [U, lambda] = eigs(DiffLap, numEvals, 'LM', struct('isreal',1,'issym',1,'maxit',500,'disp',0));
    evals(j,:) = sort(diag(lambda), 'descend');
    
    for cc=1:cback
        fprintf('\b');
    end
    cback = fprintf(['%4d/' num2str(GroupSize) ' done.\n'],j);
end
Names = taxa_code(TAXAinds); % match upper/lower cases

%% draw box plot
boxplot(evals);
title('HDM Dataset, Spectrum of Diffusion Laplacians');
