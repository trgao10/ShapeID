clear vas
% close all
path(pathdef);
addpath(path,genpath([pwd '/utils/']));
% addpath(path,'~/Documents/MATLAB/jsonlab/');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% load disance matrix
%%%% write here how to extract distMatrix from input data

%%%%%% Option 1: cP distances
% load('data/cPDistMatrix.mat')
% distMatrix = cPDistMatrix - diag(diag(cPDistMatrix));
% distMatrix = (distMatrix + distMatrix')/2;

%%%%%% Option 2: HDBM
load('data/PNAS_HBDM_cPMST_FeatureFixOff_BNN5.mat');
% load('data/PNAS_HDBM_cPMST_FeatureFixOn_BNN3.mat')
% load('data/PNAS_HDBM_cPComposedLASTbalance_FeatureFixOn_BNN3.mat')
distMatrix = squareform(pdist(HDBM(:,1:end)));

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
ColorCodeName = TB.ColorCodeName;
Genus2 = TB.Genus2;
Generic = TB.Generic;

Label = Generic;

TaxaPath = [DataPath 'teeth_taxa_table.mat'];
taxa_code = load(TaxaPath);
taxa_code = taxa_code.taxa_code;
IndsTBtoTaxaCode = zeros(size(distMatrix, 1),1);
for j=1:size(distMatrix, 1)
    IndsTBtoTaxaCode(j) = find(strcmpi(taxa_code, Name{j}));
end

%%%%% define the indices to be plotted here
uniqueLabels = unique(Label);
idx = cell(size(uniqueLabels));
for j=1:length(uniqueLabels)
    %%%%% this is index in the taxa table (not taxa_code)
    idx{j} = IndsTBtoTaxaCode(strcmpi(Label, uniqueLabels{j}));
%     idx{j} = IndsTB(strcmpi(Label, uniqueLabels{j}));
end

groupLabels = cell(1, length(Name));
for j=1:length(uniqueLabels)
    for k=1:length(idx{j})
        groupLabels{idx{j}(k)} = uniqueLabels{j};
    end
end

genericLevelColors = [1,0,0;     %% Adapis
                      0,0.5,0;   %% Arctocebus
                      1,0,0;     %% Cantius
                      0,0.5,0.5; %% Cheirogaleidae
                      0,0.5,0.5; %% Cheirogaleidae
                      1,1,0;     %% Cynocephalus
                      1,0,0;     %% Donrussellia
                      0,0,0;     %% Eosimias
                      0,1,0;     %% Galago
                      0,1,1;     %% Indridae
                      0,0,1;     %% Lemuridae
                      0,0,0;     %% Lepilemur
                      0.6,0.4,0; %% Leptacodon
                      0,0.5,0;   %% Loris
                      0,0,0;     %% Megaladapis
                      0,0.5,0;   %% Nycticebus
                      0,0.5,0;   %% Perodicticus
                      0.6,0.4,0; %% Plesiadapoidea
                      0,0,1;     %% Prolemur
                      0.6,0.4,0; %% Purgatorius
                      1,1,0;     %% Scandentia
                      0.5,0,0.5; %% Tarsius
                      1,0,1;     %% Teilhardina
                      ];
colors = genericLevelColors;
Y(:,1) = -Y(:,1);
for j = 1:length(uniqueLabels)
    X1 = Y(idx{j}, 1);
    X2 = Y(idx{j}, 2);
    if j == 1
        hold on
    end
    scatter(X1, X2, 25, colors(j, :), 'filled');
end
legend(uniqueLabels, 'location', 'northeastoutside');
for j = 1:length(uniqueLabels)
    if length(idx{j}) < 2
        continue;
    end
    X1 = Y(idx{j}, 1);
    X2 = Y(idx{j}, 2);
    if length(idx{j}) == 2
        plot(X1, X2, 'Color', colors(j, :));
    else
        cHull = convhull(X1, X2);
        patch(X1(cHull), X2(cHull), colors(j, :), 'FaceAlpha', .6);
        plot(X1(cHull), X2(cHull), 'Color', colors(j, :));
    end
end
title('Horizontal Base Diffusion Distance');
