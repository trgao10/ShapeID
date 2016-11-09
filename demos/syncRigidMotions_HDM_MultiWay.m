%% preparation
close all;
path(pathdef);
addpath(path,genpath([pwd '/utils/']));

foliIdx = [TaxaByGroup{1};TaxaByGroup{3}];
frugIdx = [TaxaByGroup{2};TaxaByGroup{4}];
inscIdx = TaxaByGroup{5};
DietLabels = zeros(GroupSize,1);
DietLabels(foliIdx) = 1;
DietLabels(frugIdx) = 2;
DietLabels(inscIdx) = 3;

colorsList = [228,26,28;0,126,204;77,175,74;152,78,163;255,127,0;255,255,51;166,86,40;247,129,191]/255;
markerList = {'+','o','x','^','d'};

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

debugFlag = false;

%%%%% gold standard
% BaseEps = 0.04;
% BNN = 4;
%%%%% experiment
% BaseEps = 0.04;
% BNN = 4;
% BNN = 6; %%% separate folivore from others
BNN = 7; %%% perfect spot for separating three clusters
% BNN = 8; %%% separate out Saimiri (insectivore), promising for detecting more
% BNN = 9; %%% separate folivore from others but with errors, same for larger
% BNN = 40;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% parameters for SynCut
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
numClusters = 3;
d = 3;
maxIter = 10;
tol = 1e-8;
numKmeans = 200;
bandwidth = 'auto'; %%%% ['auto' | some number]
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
% BaseWeights = sparse(repmat((1:GroupSize)',1,BNN),rowNNs,exp(-sDists.^2/BaseEps));
BaseWeights = sparse(repmat((1:GroupSize)',1,BNN),rowNNs,sDists);
BaseWeights = max(BaseWeights, BaseWeights');

tD = BaseDistMatrix;
tD = tD+diag(Inf(GroupSize,1));
% epsilon = BaseEps;
epsilon = mean(min(tD, [] ,2));
W = exp(-tD.^2/epsilon^2);
W = W.*double(BaseWeights>0);
D = sum(W,2);
L = diag(1./sqrt(D))*W*diag(1./sqrt(D));
L = (L+L')/2;

[evecs, evals] = eig(L);
devals = 1-diag(evals);
[devals,evalidx] = sort(devals,'ascend');
evecs = evecs(:,evalidx);
Ydm = diag(1./sqrt(D))*evecs(:,2:(numClusters+5))*sqrt(diag(1-devals(2:(numClusters+5))));

embedPts = Ydm(:,1:3);
three_cluster_idx_dm = kmeans(embedPts,numClusters,'MaxIter',1000,'Replicates',numKmeans,'Display','off');
three_clusterLabel_dm = cell(1,numClusters);
for j=1:numClusters
    three_clusterLabel_dm{j} = find(three_cluster_idx_dm == j);
end
five_cluster_idx_dm = kmeans(embedPts,5,'MaxIter',1000,'Replicates',numKmeans,'Display','off');
five_clusterLabel_dm = cell(1,5);
for j=1:5
    five_clusterLabel_dm{j} = find(five_cluster_idx_dm == j);
end


planeFigure = figure;
subplot(1,2,1);
scatter3(Ydm(foliIdx,1),Ydm(foliIdx,2),Ydm(foliIdx,3),20,'r','filled');
hold on
scatter3(Ydm(frugIdx,1),Ydm(frugIdx,2),Ydm(frugIdx,3),20,'b','filled');
scatter3(Ydm(inscIdx,1),Ydm(inscIdx,2),Ydm(inscIdx,3),20,'g','filled');
title('Diffusion Maps');

% keyboard

%%%%% build full graph
% G.adjMat = W;
G.adjMat = full(BaseWeights);
G.adjMat = double(G.adjMat > 0);
G.adjMat = G.adjMat.*W;

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

%% run multi-way SynCut
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

rslt = SynCut(G, RCell, params);

five_cluster_idx = kmeans(rslt.embedPts,5,'MaxIter',1000,'Replicates',numKmeans,'Display','off');
rslt.five_clusterLabel = cell(1,5);
for j=1:5
    rslt.five_clusterLabel{j} = find(five_cluster_idx == j);
end


%% generate statistics for clustering results
figure(planeFigure);
subplot(1,2,2);
scatter3(rslt.embedPts(foliIdx,1),rslt.embedPts(foliIdx,2),rslt.embedPts(foliIdx,3),20,'r','filled');
hold on
scatter3(rslt.embedPts(frugIdx,1),rslt.embedPts(frugIdx,2),rslt.embedPts(frugIdx,3),20,'b','filled');
scatter3(rslt.embedPts(inscIdx,1),rslt.embedPts(inscIdx,2),rslt.embedPts(inscIdx,3),20,'g','filled');
title('SynCut');
% scatter3(rslt.embedPtsFull(foliIdx,2),rslt.embedPtsFull(foliIdx,3),rslt.embedPtsFull(foliIdx,4),20,'r','filled');
% hold on
% scatter3(rslt.embedPtsFull(frugIdx,2),rslt.embedPtsFull(frugIdx,3),rslt.embedPtsFull(frugIdx,4),20,'b','filled');
% scatter3(rslt.embedPtsFull(inscIdx,2),rslt.embedPtsFull(inscIdx,3),rslt.embedPtsFull(inscIdx,4),20,'g','filled');
% title('SynCut');

% embedPts = rslt.embedPts;
embedPts = tsne(rslt.embedPts(:,1:3),[],rslt.embedPts(:,1:3),10);
% embedPts = tsne(rslt.embedPtsFull(:,2:4),[],rslt.embedPts,10);
% embedPts = tsne(rslt.embedPtsFull(:,2:4),DietLabels,rslt.embedPts,10);

compareFig = figure;
subplot(1,2,2);
plot(embedPts(foliIdx,1), embedPts(foliIdx,2), '+', 'Color', colorsList(1,:));
hold on
plot(embedPts(frugIdx,1), embedPts(frugIdx,2), 'o', 'Color', colorsList(2,:));
plot(embedPts(inscIdx,1), embedPts(inscIdx,2), 'x', 'Color', colorsList(3,:));
axis equal
% axis square
legend('Folivore', 'Frugivore', 'Insectivore', 'Location', 'NorthWest');

speciesFig = figure;
subplot(1,2,2);
for j=1:length(GroupNames)
    plot(embedPts(TaxaByGroup{j},1), embedPts(TaxaByGroup{j},2), markerList{j}, 'Color', colorsList(j,:));
    if j==1
        hold on
    end
end
% hold on
% plot(embedPts(frugIdx,1), embedPts(frugIdx,2), 'o', 'Color', colorsList(2,:));
% plot(embedPts(inscIdx,1), embedPts(inscIdx,2), 'x', 'Color', colorsList(3,:));
axis equal
axis square
legend(GroupNames, 'Location', 'NorthWest');
title('\textbf{(b) SynCut}', 'Interpreter', 'Latex', 'FontSize', 18);

% boldify
% subplot(1,2,1);
% Childs = (get(gca, 'Children'));
% for j=1:length(Childs)
%     try
%         if strcmpi(get(Childs(j),'Type'),'line')
%             set(Childs(j),'LineWidth',1.0);
%         end
%     catch ME
%         continue
%     end
% end
% subplot(1,2,2);
% Childs = (get(gca, 'Children'));
% for j=1:length(Childs)
%     try
%         if strcmpi(get(Childs(j),'Type'),'line')
%             set(Childs(j),'LineWidth',1.0);
%         end
%     catch ME
%         continue
%     end
% end

%% compare with diffusion maps
% tD = BaseDistMatrix;
% tD = tD+diag(Inf(GroupSize,1));
% % epsilon = BaseEps;
% epsilon = mean(min(tD, [] ,2));
% W1 = exp(-tD.^2/epsilon^2);
% D = sum(W1,2);
% L = diag(1./sqrt(D))*W1*diag(1./sqrt(D));
% L = (L+L')/2;
% 
% [evecs, evals] = eig(L);
% devals = 1-diag(evals);
% [devals,evalidx] = sort(devals,'ascend');
% evecs = evecs(:,evalidx);
% Ydm = diag(1./sqrt(D))*evecs(:,1:(numClusters+1));
% embedPtsFull = diag(1./sqrt(D))*evecs;

% embedPts = tsne(Ydm(:,1:3),[],Ydm(:,1:3),10);
% embedPts = embedPts(:,[1 3 2]);
embedPts = tsne(Ydm(:,1:3),[],rslt.embedPts,10);
% figure;
% embedPts = tsne(Ydm(:,2:4),[],rslt.embedPts,10);
% embedPts = tsne(Ydm(:,2:4),DietLabels,rslt.embedPts,10);
% close(gcf);
% 
% % figure;
% % scatter3(Ydm(foliIdx,2),Ydm(foliIdx,3),Ydm(foliIdx,4),20,'r','filled');
% % hold on
% % scatter3(Ydm(frugIdx,2),Ydm(frugIdx,3),Ydm(frugIdx,4),20,'b','filled');
% % scatter3(Ydm(inscIdx,2),Ydm(inscIdx,3),Ydm(inscIdx,4),20,'g','filled');
% 
% figure(planeFigure)
% subplot(1,2,1);
% scatter3(Ydm(foliIdx,2),Ydm(foliIdx,3),Ydm(foliIdx,4),20,'r','filled');
% hold on
% scatter3(Ydm(frugIdx,2),Ydm(frugIdx,3),Ydm(frugIdx,4),20,'b','filled');
% scatter3(Ydm(inscIdx,2),Ydm(inscIdx,3),Ydm(inscIdx,4),20,'g','filled');


figure(compareFig);
subplot(1,2,1);
plot(embedPts(foliIdx,1), embedPts(foliIdx,2), '+', 'Color', colorsList(1,:));
hold on
plot(embedPts(frugIdx,1), embedPts(frugIdx,2), 'o', 'Color', colorsList(2,:));
plot(embedPts(inscIdx,1), embedPts(inscIdx,2), 'x', 'Color', colorsList(3,:));
axis equal
legend('Folivore', 'Frugivore', 'Insectivore', 'Location', 'NorthWest');


figure(speciesFig);
subplot(1,2,1);
for j=1:length(GroupNames)
    plot(embedPts(TaxaByGroup{j},1), embedPts(TaxaByGroup{j},2), markerList{j}, 'Color', colorsList(j,:));
    if j==1
        hold on
    end
end
% hold on
% plot(embedPts(frugIdx,1), embedPts(frugIdx,2), 'o', 'Color', colorsList(2,:));
% plot(embedPts(inscIdx,1), embedPts(inscIdx,2), 'x', 'Color', colorsList(3,:));
axis equal
axis square
legend(GroupNames, 'Location', 'NorthWest');
title('\textbf{(a) Diffusion Maps}', 'Interpreter', 'Latex', 'FontSize', 18);

% [Udm, Ldm] = eig(L);
% dims = sum(diag(Ldm)>0);
% Udm = Udm(:,1:dims);
% Ldm = Ldm(1:dims, 1:dims);
% Ydm = diag(1./sqrt(D))*Udm*sqrt(Ldm);

% Ydm = tsne(Ydm, [], 3, 10, 30);
% DM_dist = pdist(Ydm);
% [Ydm,stress] = mdscale(DM_dist,3,'criterion','metricstress');
% 
% figure('Name','Diffusion Maps');
% for j=1:length(GroupNames)
%     plot3(Ydm(PerGroupDelimit(j,1):PerGroupDelimit(j,2),1),...
%         Ydm(PerGroupDelimit(j,1):PerGroupDelimit(j,2),2),...
%         Ydm(PerGroupDelimit(j,1):PerGroupDelimit(j,2),3),...
%         'Color', colorsList(j,:), 'Marker', '.', 'MarkerSize', 15, 'LineStyle', 'none');
%     if (j == 1)
%         grid on
%         axis equal
%         hold on;
%     end
% end
% legend(GroupNames);

% scatter(embedPts(foliIdx,1), embedPts(foliIdx,2), 20, 'r', 'filled');
% hold on
% scatter(embedPts(frugIdx,1), embedPts(frugIdx,2), 20, 'b', 'filled');
% scatter(embedPts(inscIdx,1), embedPts(inscIdx,2), 20, 'g', 'filled');

% scatter3(embedPts(foliIdx,1), embedPts(foliIdx,2), embedPts(foliIdx,3), 20, 'r', 'filled');
% hold on
% scatter3(embedPts(frugIdx,1), embedPts(frugIdx,2), embedPts(frugIdx,3), 20, 'b', 'filled');
% scatter3(embedPts(inscIdx,1), embedPts(inscIdx,2), embedPts(inscIdx,3), 20, 'g', 'filled');


