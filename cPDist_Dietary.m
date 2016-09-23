clear vas
% close all
path(pathdef);
addpath(path,genpath([pwd '/utils/']));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% load disance matrix
%%%% write here how to extract distMatrix from input data

%%%%%% Option 1: cP distances
load('data/cPDistMatrix.mat')
distMatrix = cPDistMatrix - diag(diag(cPDistMatrix));
distMatrix = (distMatrix + distMatrix')/2;

%%%%%% Option 2: HDBM
% load('data/PNAS_HDBM_cPMST_FeatureFixOn_BNN3.mat')
% % load('data/PNAS_HDBM_cPComposedLASTbalance_FeatureFixOn_BNN3.mat')
% distMatrix = squareform(pdist(HDBM(:,1:end)));

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
% idx1 = IndsTB(strcmpi(IGroup, 'Prosimian'));
% idx2 = IndsTB(strcmpi(IGroup, 'Plesiadapiform, "proto-primate"'));
% idx2 = setdiff(1:size(distMatrix, 1), idx1);
% idx2 = IndsTB(strcmpi(IGroup, 'Anthropoid'));
% idx1 = IndsTB(strcmpi(Infraorder, 'Lemuriformes'));
% idx2 = IndsTB(strcmpi(Infraorder, 'Lorisiformes'));
% idx3 = IndsTB(strcmpi(IGroup, 'Plesiadapiform, "proto-primate"'));
idx1 = IndsTB(strcmpi(DietaryCategory, 'Folivore'));
idx2 = IndsTB(strcmpi(DietaryCategory, 'Frugivore'));
idx3 = IndsTB(strcmpi(DietaryCategory, 'Frugivore-Folivore'));
idx4 = IndsTB(strcmpi(DietaryCategory, 'Insectivore'));
idx5 = IndsTB(strcmpi(DietaryCategory, 'Omnivore'));

%%%% 3D Plot
figure('Toolbar','none');
meanDistortion = mean(squareform(distMatrix)./pdist(Y)) * mean(pdist(Y)./squareform(distMatrix));
set(gcf, 'Name', ['Average Distance Distortion = ' num2str(meanDistortion)]);
scatter3(Y(idx1,1),Y(idx1,2),Y(idx1,3),30,[1,0,0],'filled');
axis equal
hold on
scatter3(Y(idx2,1),Y(idx2,2),Y(idx2,3),30,[1,0,0],'filled');
scatter3(Y(idx3,1),Y(idx3,2),Y(idx3,3),30,[1,0,0],'filled');
scatter3(Y(idx4,1),Y(idx4,2),Y(idx4,3),30,'b','filled');
scatter3(Y(idx5,1),Y(idx5,2),Y(idx5,3),30,'g','filled');
legend('Folivore', 'Frugivore', 'Frugivore-Folivore', 'Insectivore', 'Omnivore');
cameratoolbar;
cameratoolbar('SetCoordSys', 'none');
