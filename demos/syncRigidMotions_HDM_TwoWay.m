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
result_path = '/home/trgao10/Work/MATLAB/ArchivedResults/HDM/cPdist/';
dist_path = [result_path 'cPDistMatrix.mat'];
maps_path = [result_path 'cPMapsMatrix.mat'];
% result_path = '/home/trgao10/Work/MATLAB/ArchivedResults/HDM/cPMST/FeatureFixOff/';
% dist_path = [result_path 'cPMSTDistMatrix.mat'];
% maps_path = [result_path 'cPMSTMapsMatrix.mat'];
GroupLevel = 'Genus';
% ChunkSize = 25;
GroupNames = {'Alouatta','Ateles','Brachyteles','Callicebus','Saimiri'};
DisplayLayout = [5,10];
% GroupNames = {'Alouatta','Callicebus'};
% DisplayLayout = [4,5];

debugFlag = true;

%%%%% gold standard
% BaseEps = 0.04;
% BNN = 4;
%%%%% experiment
BaseEps = 0.04;
% BNN = 6; %%% separate folivore from others
% BNN = 7; %%% seperate folivore from others
BNN = 8; %%% separate out Saimiri (insectivore), promising for detecting more
% BNN = 9; %%% separate folivore from others but with errors, same for larger

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% parameters for SynCut
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
numClusters = 3;
d = 3;
maxIter = 10;
tol = 1e-8;
numKmeans = 200;
bandwidth = 'auto'; %%%% 'auto' | some number
adjType = 'dis';
colorList = {'r','b','k','m'};
hsv = rgb2hsv(winter);
close(gcf);

%% visualization options
options = struct('sample_path', sample_path,...
                 'DisplayLayout', DisplayLayout,...
                 'DisplayOrient', 'Horizontal',...
                 'boundary', 'on', 'names', 'off',...
                 'linkCamera', 'on', 'Shading', 'Smooth');
                 
% options.sample_path = sample_path;
% options.DisplayLayout = [2,10];
% options.DisplayOrient = 'Horizontal';
% options.boundary = 'on';
% options.names = 'off';
% options.linkCamera = 'on';
% options.Shading = 'Smooth';

%% load mapsMatrix
load(maps_path);
mapsMatrix = cPMapsMatrix;
% mapsMatrix = ImprMapsMatrix;

%% load taxa codes
taxa_file = [data_path 'hdm_taxa_table.mat'];
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

%% load meshes
TAXAinds = zeros(GroupSize,1);
meshList = cell(1,GroupSize);
for j=1:GroupSize
    TAXAinds(j) = find(strcmpi(taxa_code,Names{j}));
    load([sample_path taxa_code{strcmpi(taxa_code,Names{j})} '.mat']);
    meshList{j} = G;
end
Names = taxa_code(TAXAinds); % match upper/lower cases

LocalByGroup = TaxaByGroup;
for j=1:length(TaxaByGroup)
    for k=1:length(TaxaByGroup{j})
        TaxaByGroup{j}(k) = find(TAXAinds == TaxaByGroup{j}(k));
    end
end

%% visualize un-aligned meshes
if debugFlag
    drawMeshList(meshList, options);
    set(gcf, 'Name', 'un-aligned');
end

%% build weighted graph from cP distance matrix
clear G
%%%%% build KNN graph
load(dist_path);
BaseDistMatrix = cPDistMatrix(TAXAinds,TAXAinds);
BaseDistMatrix = BaseDistMatrix-diag(diag(BaseDistMatrix));

if BNN > (GroupSize-1)
    warning('BNN larger than allowed by sample size');
    BNN = min(BNN,GroupSize-1);
end

%%% only connect BNN-nearest-neighbors
[sDists,rowNNs] = sort(BaseDistMatrix,2);
sDists = sDists(:,2:min(1+BNN,GroupSize));
rowNNs = rowNNs(:,2:min(1+BNN,GroupSize));
BaseWeights = sparse(repmat((1:GroupSize)',1,BNN),rowNNs,sDists);
BaseWeights = max(BaseWeights, BaseWeights');
for j=1:GroupSize
    sDists(j,:) = BaseWeights(j,rowNNs(j,:));
end
sDists = exp(-sDists.^2/BaseEps);

%%%%% build full graph
% G.adjMat = ones(GroupSize) - eye(GroupSize);

G.adjMat = full(BaseWeights);
% G.adjMat = double(G.adjMat > 0);

adjMask = G.adjMat;
G.V = rand(GroupSize,2);
for j=1:length(TaxaByGroup)
    blockIdx = TaxaByGroup{j};
    adjMask(blockIdx,blockIdx) = 0;
    G.V(blockIdx,1) = G.V(blockIdx,1)+(j-1);
    %%%% shrink cluster to its centroid for better visualization quality
    Vblock = G.V(blockIdx,:);
    VCenter = mean(Vblock,1);
    Vblock = (Vblock-repmat(VCenter,size(Vblock,1),1))*0.8 + ...
        repmat(VCenter,size(Vblock,1),1);
    G.V(blockIdx,:) = Vblock;
end

[G.ccRowIdx,G.ccColIdx] = find(triu(adjMask));

Dvec = sum(G.adjMat);
GL = diag(1./sqrt(Dvec))*G.adjMat*diag(1./sqrt(Dvec));
evals = eig(GL);
evals = 1-evals;
evals = sort(evals);
G.specGap = evals(2);
G.numClusters = numClusters;

%% load rigid motions output from cPDist
R = cell(GroupSize,GroupSize);
cback = 0;
for j=1:GroupSize
    for k = 1:GroupSize
        [~,R{k,j},~] = MapToDist(meshList{j}.V,meshList{k}.V,...
            mapsMatrix{TAXAinds(j),TAXAinds(k)},meshList{j}.Aux.VertArea);
    end
    for cc=1:cback
        fprintf('\b');
    end
    cback = fprintf(['%4d/' num2str(GroupSize) ' done.\n'],j);
end

% rigid_motions = load([data_path 'cPMSTinitRms.mat']);
% R = rigid_motions.R(TAXAinds,TAXAinds);

[wGCL, wGCL_Dvec, wGCL_W] = assembleGCL(G.adjMat, R, d);
if ~issymmetric(wGCL)
    warning('wGCL not symmetric!');
    wGCL = triu(wGCL,1);
    wGCL = wGCL+wGCL'+eye(size(wGCL));
end
if ~issymmetric(wGCL_W)
    warning('wGCL_W not symmetric!');
    wGCL_W = triu(wGCL_W,1);
    wGCL_W = wGCL_W+wGCL_W';
end
[RelaxSolCell, RelaxSolMat] = syncSpecRelax(wGCL, d, wGCL_Dvec);

%%%%% transpose each element to obtain the correct solution
SCSolCell = cellfun(@(x) x', RelaxSolCell, 'UniformOutput', false);

%%%%% Option 1: build RCell from R directly
RCell = cell(size(R));
for j=1:size(RCell,1)
    for k=(j+1):size(RCell,2)
        RCell{j,k} = R{j,k};
        RCell{k,j} = RCell{j,k}';
    end
end
%%%%% Option 2: set RCell toR
% RCell = R;
%%%%% Option 3: set RCell to RelaxSolCell
% RCell = RelaxSolCell;

vertPotCell = R(1,:);
GTTotalFrust = 0;
RSTotalFrust = 0;
SCTotalFrust = 0;
for j=1:size(RCell,1)
    for k=(j+1):size(RCell,2)
        GTTotalFrust = GTTotalFrust+G.adjMat(j,k)*norm(vertPotCell{j}'*vertPotCell{k}-RCell{j,k}, 'fro')^2;
        RSTotalFrust = RSTotalFrust+G.adjMat(j,k)*norm(RelaxSolCell{j}*RelaxSolCell{k}'-RCell{j,k}, 'fro')^2;
        SCTotalFrust = SCTotalFrust+G.adjMat(j,k)*norm(SCSolCell{j}'*SCSolCell{k}-RCell{j,k}, 'fro')^2;
    end
end

fprintf('[GroundTruth] GTTotal = %f\n', GTTotalFrust);
fprintf('[RelaxSol] RSTotal = %f\n', RSTotalFrust);
fprintf('[SCSol] SCTotal = %f\n', SCTotalFrust);

[GTSolPerEdgeFrustVec, GTSolPerEdgeFrustMat] =...
    getPerEdgeFrustFromEdgePot(G.adjMat, RCell, cellfun(@(x) x', vertPotCell, 'UniformOutput', false));
[RSSolPerEdgeFrustVec, RSSolPerEdgeFrustMat] =...
    getPerEdgeFrustFromEdgePot(G.adjMat, RCell, RelaxSolCell);
[SCSolPerEdgeFrustVec, SCSolPerEdgeFrustMat] =...
    getPerEdgeFrustFromEdgePot(G.adjMat, RCell, cellfun(@(x) x', SCSolCell, 'UniformOutput', false));
fprintf('[GroundTruth] verify GTTotal = %f\n', sum(GTSolPerEdgeFrustVec));
fprintf('[RelaxSol] verify RSTotal = %f\n', sum(RSSolPerEdgeFrustVec));
fprintf('[SCSol] verify SCTotal = %f\n', sum(SCSolPerEdgeFrustVec));

% figure;
% subplot(1,3,1);
% imagesc(G.adjMat);
% axis square
% title('weighted adjacency matrix');
% subplot(1,3,2);
% imagesc(GTSolPerEdgeFrustMat);
% axis square
% title('GroundTruth Edgewise Frust');
% subplot(1,3,3);
% imagesc(SCSolPerEdgeFrustMat);
% axis square
% title('RelaxSol Edgewise Frust');

% keyboard

% RMat = zeros(size(RCell)*d);
% for j=1:size(RCell,1)
%     blockRowIdx = ((j-1)*d+1):(j*d);
%     for k=1:size(RCell,2)
%         if ~isempty(RCell{j,k})
%             blockColIdx = ((k-1)*d+1):(k*d);
%             RMat(blockRowIdx,blockColIdx) = RCell{j,k};
%         end
%     end
% end

%% visualize aligned meshes
% visCell = vertPotCell;
% visCell = RelaxSolCell;
visCell = SCSolCell;
if debugFlag
    for j=1:GroupSize
        meshList{j}.V = visCell{j} * meshList{j}.V;
        if det(visCell{j}) < 0
            meshList{j}.F = meshList{j}.F([1 3 2], :);
        end
    end
    drawMeshList(meshList, options);
    set(gcf, 'Name', 'aligned');
end

% keyboard

% %%
% TB = readtable([data_path 'ClassificationTable.xlsx']);
% TBName = TB.Name;
% Genus = TB.Genus;
% Chirality = TB.Chirality;
% 
% TaxaPath = [data_path 'hdm_taxa_table.mat'];
% taxa_code = load(TaxaPath);
% taxa_code = taxa_code.taxa_code;
% IndsTBtoTaxaCode = zeros(length(Names),1);
% for j=1:length(TBName)
%     IndsTBtoTaxaCode(j) = find(strcmpi(taxa_code, TBName{j}));
% end
% 
% %%%%% define the indices to be plotted here
% idx1 = intersect(IndsTBtoTaxaCode(strcmpi(Chirality, 'Left')), TAXAinds);
% idx2 = intersect(IndsTBtoTaxaCode(strcmpi(Chirality, 'Right')), TAXAinds);
% 
% idx1ToSample = zeros(size(idx1));
% for j=1:length(idx1ToSample)
%     idx1ToSample(j) = find(idx1(j) == TAXAinds);
% end
% idx2ToSample = zeros(size(idx2));
% for j=1:length(idx2ToSample)
%     idx2ToSample(j) = find(idx2(j) == TAXAinds);
% end

%% run two-way SynCut
debugFlag = false;

params = struct('debugFlag', debugFlag,...
                'd', d,...
                'numClusters', numClusters,...
                'tol', tol,...
                'maxIter', maxIter,...
                'numKmeans', numKmeans,...
                'bandwidth', bandwidth,...
                'adjType', adjType,...
                'hsv', hsv);
[params.vertPotCell] = deal(vertPotCell);

rslt = SynCut_TwoWay_HDM(G, RCell, params);


