clear vas
close all
path(pathdef);
addpath(path,genpath([pwd '/utils/']));

% Name = 'a10';
% G = Mesh('off', ['/media/trgao10/Work/MATLAB/DATA/PNAS/meshes_Poisson_sphere/sphere/' Name '_sas.off']);
% G.Centralize('ScaleArea');
% NumLandmark = 16;
% [~,Landmarks] = G.GetLandmarksFromPNAS('/media/trgao10/Work/MATLAB/DATA/PNAS/landmarks_teeth.mat','/media/trgao10/Work/MATLAB/DATA/PNAS/meshes/a10_sas.off',NumLandmark);
% tree = KDTreeSearcher(G.V');
% LandmarkInds = tree.knnsearch(Landmarks);
% Landmarks = G.V(:,LandmarkInds)';

G = Mesh('off', '/home/trgao10/Downloads/MNHNQu_11046.off');
VertSampInd = G.GeodesicFarthestPointSampling(5000);

subsampledVertices = G.V(:,VertSampInd);
scatter3(subsampledVertices(1,:),subsampledVertices(2,:),subsampledVertices(3,:),20,'g','filled');

[Vs,Fs] = G.PerformGeodesicDelauneyTriangulation(VertSampInd,[]);
Gsub = Mesh('VF', Vs, Fs);
Gsub.DeleteIsolatedVertex();

Gsub.Write('MNHNQu_11046_sub5000.off','off',[]);

