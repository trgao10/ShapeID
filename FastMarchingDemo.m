clearvars;
close all;
path(pathdef);
path(path, genpath('./utils'));

samples_path = './samples/Clement/';
MeshName = '11';
load([samples_path MeshName '.mat']);
% G = Mesh('off','./meshes/david-head.off');

%%% pick vertices of interest
%%% by invoking the picking function of the display
% G.draw;
% keyboard

%%% fast marching from initial points of interest
D = G.PerformFastMarching(70);
figure;G.ViewFunctionOnMesh(D,struct('mode','native'));

% Subsamples = G.GeodesicFarthestPointSampling(400);

%%% geodesic extraction demo
% geodesic = G.ComputeGeodesic(D, 1504);
% hold on
% plot3(geodesic(1,:),geodesic(2,:),geodesic(3,:),'k-','LineWidth',2);

%%% anisotropic fast marching
[K,M] = G.ComputeCurvature;
% G.ViewFunctionOnMesh(M,struct('mode','native'));
D = G.PerformFastMarching(70,struct('W',M-min(M)));
figure;G.ViewFunctionOnMesh(D,struct('mode','native'));



