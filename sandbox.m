clear vas
close all
path(pathdef);
addpath(path,genpath([pwd '/utils/']));

G = Mesh('off','data/BTL/BTL_26102016_113231_in_mm_conncomp_HCLaplacianSmooth6.off');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% extract largest connected component
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% bins = conncomp(graph(G.A));
% uniqueBins = unique(bins);
% counts = zeros(1,length(uniqueBins));
% for j=1:length(counts)
%     counts(j) = sum(bins == uniqueBins(j));
% end
% [maxCount,maxIdx] = max(counts);
% delIdx = 1:G.nV;
% delIdx = delIdx(bins ~= uniqueBins(maxIdx));
% G.DeleteVertex(delIdx);
% 
% G.draw();
% G.Write('BTL_26102016_113211_in_mm_conncomp.off','off');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% HC Laplacian smoothing 6 times (in meshlab)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% check duplicate faces
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[F,Inds] = unique(sort(G.F',2),'rows','first');
if ~isempty(setdiff(1:G.nF,Inds))
    disp('found duplicate faces!');
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% check duplicate vertices
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[V,Inds] = unique(sort(G.V',2),'rows','first');
if ~isempty(setdiff(1:G.nV,Inds))
    disp('found duplicate vertices!');
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% compute Laplace-Beltrami eigenfunctions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
G.Aux.LB = G.ComputeCotanLaplacian();

eigopt = struct('isreal',1,'issym',1,'maxit',5000,'disp',0);
[evecs, evals] = eigs(-G.Aux.LB, 100, 'SM', eigopt);
devals = diag(evals);
[devals,sortIdx] = sort(devals,'ascend');
evecs = evecs(:,sortIdx);

%% view an eigenfunction
figure;
dispLayout = [2,4];
h = zeros(1,prod(dispLayout));
for j=1:prod(dispLayout)
    h(j) = subplot(dispLayout(1),dispLayout(2),j);
    G.ViewFunctionOnMesh(evecs(:,j+1),struct('mode','native'));
end
Link = linkprop(h, {'CameraUpVector', 'CameraPosition', 'CameraTarget', 'CameraViewAngle'});
setappdata(gcf, 'StoreTheLink', Link);

% dims = sum(diag(Ldm)>0);
% Udm = Udm(:,1:dims);
% Ldm = Ldm(1:dims, 1:dims);


% Name = 'a10';
% G = Mesh('off', ['/media/trgao10/Work/MATLAB/DATA/PNAS/meshes_Poisson_sphere/sphere/' Name '_sas.off']);
% G.Centralize('ScaleArea');
% NumLandmark = 16;
% [~,Landmarks] = G.GetLandmarksFromPNAS('/media/trgao10/Work/MATLAB/DATA/PNAS/landmarks_teeth.mat','/media/trgao10/Work/MATLAB/DATA/PNAS/meshes/a10_sas.off',NumLandmark);
% tree = KDTreeSearcher(G.V');
% LandmarkInds = tree.knnsearch(Landmarks);
% Landmarks = G.V(:,LandmarkInds)';

% G = Mesh('off', '/home/trgao10/Downloads/MNHNQu_11046.off');
% VertSampInd = G.GeodesicFarthestPointSampling(5000);
% 
% subsampledVertices = G.V(:,VertSampInd);
% scatter3(subsampledVertices(1,:),subsampledVertices(2,:),subsampledVertices(3,:),20,'g','filled');
% 
% [Vs,Fs] = G.PerformGeodesicDelauneyTriangulation(VertSampInd,[]);
% Gsub = Mesh('VF', Vs, Fs);
% Gsub.DeleteIsolatedVertex();
% 
% Gsub.Write('MNHNQu_11046_sub5000.off','off',[]);

