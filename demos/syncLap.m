clearvars
close all
path(pathdef)
addpath(path,genpath([pwd '/utils/']));

numMeshes = 8;
d = 3;
DisplayLayout = [2,4];
DataPath = '~/Work/MATLAB/DATA/PNAS/';
SamplePath = [DataPath 'samples/'];
TaxaPath = [DataPath 'teeth_taxa_table.mat'];
LandmarkPath = [DataPath 'landmarks_teeth.mat'];
MeshPath = [DataPath 'meshes/'];
MeshSuffix = '_sas.off';

taxa_code = load(TaxaPath);
taxa_code = taxa_code.taxa_code;

rng(1);
randseq = randperm(length(taxa_code));
TaxaInds = randseq(1:numMeshes);

meshList = cell(numMeshes,1);
for j=1:numMeshes
    load([SamplePath taxa_code{TaxaInds(j)} '.mat']);
    [~,ObCoords] = GetLandmarks(taxa_code{strcmpi(taxa_code, taxa_code{TaxaInds(j)})},...
        LandmarkPath, [MeshPath taxa_code{strcmpi(taxa_code, taxa_code{TaxaInds(j)})} MeshSuffix]);
    Vtree = kdtree_build( G.V' );
    G.Aux.ObInds = kdtree_nearest_neighbor(Vtree, ObCoords);
    meshList{j} = G;
end

% T = zeros(numMeshes);
T = zeros(numMeshes*d);

for j=1:numMeshes
    Pts1 = meshList{j}.V(:, meshList{j}.Aux.ObInds);
    for k=(j+1):numMeshes
        % Procrustes Distance between two meshes
        Pts2 = meshList{k}.V(:, meshList{k}.Aux.ObInds);
        % in principle, no translation is needed since all teeth are
        % already centered at the origin
        [U,~,V] = svd(Pts1*Pts2');
        R = V*U';
        
        rowIdx = ((j-1)*d + 1):(j*d);
        colIdx = ((k-1)*d + 1):(k*d);
        T(rowIdx, colIdx) = R;
%         T(j,k) = det(R);
    end
end

T = T+T'+eye(numMeshes*d);

L = 1-T;
D = diag(sum(L));
eigopt = struct('isreal',1,'issym',1,'maxit',5000,'disp',0);
W = sparse(diag(sqrt(diag(D)))*(eye(numMeshes*d) - inv(D)*L)*diag(1./sqrt(diag(D))));
W = (W+W')/2;
[U, lambda] = eigs(W, numMeshes-1, 'LM', eigopt);
[lambda, ia] = sort(diag(lambda), 'descend');
figure;plot(1:(numMeshes-1), lambda, '+'); % the decay of eigenvalues
title('Eigenvalues of the graph Laplacian');

% idx = kmeans(sign(U(:,1)),2,'MaxIter',1000);
idx = kmeans(reshape(U(:,1), [d,numMeshes])',2,'MaxIter',1000);

disp('class labels of each tooth: ');
disp(idx');
