clearvars;
close all;
path(pathdef);
path(path, genpath('./utils'));

% G = Mesh('off','data/sphere258.off');

samples_path = '~/Work/MATLAB/DATA/PNAS/samples/';
MeshName = 'k11';
load([samples_path MeshName '.mat']);
% G = Mesh('off','/home/trgao10/Work/MATLAB/DATA/HDM/sphTeethMinSurf/AMNH-M-98272_M912_close.off');

adjMat = G.ComputeAdjacencyMatrix();
[rIdx,cIdx] = find(adjMat);
Dists = sqrt(sum((G.V(:,rIdx) - G.V(:,cIdx)).^2));
bandwidth = mean(Dists).^2;
alpha = 0;

[DiffLB,DVec] = G.ComputeDiffusionLaplacian(bandwidth,alpha);

eigopt = struct('isreal',1,'issym',1,'maxit',5000,'disp',0);
[evecs, evals] = eigs(DiffLB, 100, 'SM', eigopt);
evecs = diag(1./sqrt(DVec))*evecs;
devals = diag(evals);
[devals,sortIdx] = sort(devals,'ascend');
evecs = evecs(:,sortIdx);

figure;
dispLayout = [3,5];
h = zeros(1,prod(dispLayout));
for j=1:prod(dispLayout)
    h(j) = subplot(dispLayout(1),dispLayout(2),j);
    G.ViewFunctionOnMesh(evecs(:,j+1),struct('mode','rb'));
end
Link = linkprop(h, {'CameraUpVector', 'CameraPosition', 'CameraTarget', 'CameraViewAngle'});
setappdata(gcf, 'StoreTheLink', Link);
