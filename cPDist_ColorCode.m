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
ColorCodeName = TB.ColorCodeName;

Label = ColorCodeName;

TaxaPath = [DataPath 'teeth_taxa_table.mat'];
taxa_code = load(TaxaPath);
taxa_code = taxa_code.taxa_code;
IndsTB = zeros(size(distMatrix, 1),1);
for j=1:size(distMatrix, 1)
    IndsTB(j) = find(strcmpi(taxa_code, Name{j}));
end

%%%%% define the indices to be plotted here
uniqueLabels = unique(Label);
idx = cell(size(uniqueLabels));
for j=1:length(uniqueLabels)
    %%%%% this is index in the taxa table (not taxa_code)
    idx{j} = IndsTB(strcmpi(Label, uniqueLabels{j}));
end

groupLabels = cell(1, length(Name));
for j=1:length(uniqueLabels)
    for k=1:length(idx{j})
        groupLabels{idx{j}(k)} = uniqueLabels{j};
    end
end

% gscatter(Y(:,1), Y(:,2), groupLabels', colormap(colors), [], 25);
% l= findobj(gcf,'tag','legend'); set(l,'location','northeastoutside');
colors = colormap(colorcube(4+length(uniqueLabels)));
% gscatter(Y(:,1), Y(:,2), groupLabels', colormap(colorcube(1.05*length(uniqueLabels))), [], 25);
% colors =  [1,0,0;0,1,0;0,0,1;1,1,0;1,0,1;0,1,1];
% colors = [colors;colors*0.7;colors*0.5];
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
    if length(idx{j}) < 3
        continue;
    end
    X1 = Y(idx{j}, 1);
    X2 = Y(idx{j}, 2);
    cHull = convhull(X1, X2);
    patch(X1(cHull), X2(cHull), colors(j, :), 'FaceAlpha', .6);
    plot(X1(cHull), X2(cHull), 'Color', colors(j, :));
end
title('Continuous Procrustes Distance');
