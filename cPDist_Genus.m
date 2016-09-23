clear vas
% close all
path(pathdef);
addpath(path,genpath([pwd '/utils/']));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% load disance matrix
%%%% write here how to extract distMatrix from input data

load('data/cPDistMatrix.mat')
distMatrix = cPDistMatrix - diag(diag(cPDistMatrix));
distMatrix = (distMatrix + distMatrix')/2;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% embed A into a Euclidean space
%%%% MDS or diffusion maps, or any other methods

%%%%%% Option 1: Diffusion Maps (or other eigen-method) on similarity matrix
% tD = distMatrix;
% tD = tD+diag(Inf(size(distMatrix, 1),1));
% 
% sD = sort(tD, 2);
% selfTuningCol = sD(:, 6);
% tD = tD./sqrt(selfTuningCol*selfTuningCol');
% 
% epsilon = mean(min(tD, [] ,2));
% W = exp(-tD.^2/epsilon^2);
% D = sum(W,2);
% L = diag(1./sqrt(D))*W*diag(1./sqrt(D));
% L = (L+L')/2;
% [Udm, Ldm] = eigs(L, 4, 'LM', struct('isreal',1,'issym',1,'maxit',100,'disp',0));
% % [Udm, Ldm] = eigs(L, 4, 'LM', struct('isreal',1,'issym',1,'maxit',100,...
% %                   'v0',ones(size(distMatrix,1),1)*0.01,...
% %                   'tol',1e-20,'p',40,'disp',1));
% % dims = sum(diag(Ldm)>0);
% Udm = Udm(:,2:end);
% Ldm = Ldm(2:end, 2:end);
% Ydm = diag(1./sqrt(D))*Udm*sqrt(Ldm);
% % Ydm = Udm*sqrt(Ldm);
% % Ydm = diag(sqrt(D))*Udm*sqrt(Ldm);
% Y = Ydm;

%%%%%%% Option 2: MDS on Distance Matrix
[mdsCoords,stress] = mdscale(distMatrix,3,'criterion','metricstress');
Y = mdsCoords;
% [cmdsCoords,stress] = cmdscale(distMatrix,3);
% Y = cmdsCoords;

%%%% check if embedding makes sense
DataPath = '../DATA/PNAS/';
TB = readtable([DataPath 'ClassificationTable.xlsx']);
Name = TB.Name;
Order = TB.Order;
Family = TB.Family;
Genus = TB.Genus;
IGroup = TB.InformalGroup;
Chirality = TB.Chirality;
Infraorder = TB.Infraorder;
DietaryCategory = TB.DietaryCategory;

TaxaPath = [DataPath 'teeth_taxa_table.mat'];
taxa_code = load(TaxaPath);
taxa_code = taxa_code.taxa_code;
IndsTB = zeros(size(distMatrix, 1),1);
for j=1:size(distMatrix, 1)
    IndsTB(j) = find(strcmpi(taxa_code, Name{j}));
end

%%%%% define the indices to be plotted here
uniqueGenus = unique(Genus);
idx = cell(size(uniqueGenus));
for j=1:length(uniqueGenus)
    idx{j} = IndsTB(strcmpi(Genus, uniqueGenus{j}));
end

groupLabels = cell(1, length(Name));
for j=1:length(uniqueGenus)
    for k=1:length(idx{j})
        groupLabels{idx{j}(k)} = uniqueGenus{j};
    end
end

gscatter(Y(:,1), Y(:,2), groupLabels', colormap(colorcube(1.05*length(uniqueGenus))), [], 25);
l= findobj(gcf,'tag','legend'); set(l,'location','northeastoutside');
title('Continuous Procrustes Distance');
