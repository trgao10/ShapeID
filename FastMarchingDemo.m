clearvars;
close all;
path(pathdef);
path(path, genpath('./utils'));

samples_path = '~/Work/MATLAB/DATA/PNAS/samples/';
MeshName = 'k11';
load([samples_path MeshName '.mat']);
[BV
% G = Mesh('off','./meshes/david-head.off');

%%% pick vertices of interest
%%% by invoking the picking function of the display
% G.draw;
% keyboard

%%% fast marching from initial points of interest
D = G.PerformFastMarching([1000,1707,200]);
G.ViewFunctionOnMesh(D,struct('mode','native'));
% 
% % Subsamples = G.GeodesicFarthestPointSampling(400);
% 
% %%% geodesic extraction demo
% % geodesic = G.ComputeGeodesic(D, 1504);
% % hold on
% % plot3(geodesic(1,:),geodesic(2,:),geodesic(3,:),'k-','LineWidth',2);
% 
% %%% anisotropic fast marching
% [K,M] = G.ComputeCurvature;
% % G.ViewFunctionOnMesh(M,struct('mode','native'));
% D = G.PerformFastMarching(1707,struct('W',M-min(M)));
% G.ViewFunctionOnMesh(D,struct('mode','native'));
% 
% close all;
% KM = G.ComputeCPMS();
% DM = KM.ComputeUniformization(struct('method','LSCM','boundary_conditions','disc'));
% 
% 
% 
