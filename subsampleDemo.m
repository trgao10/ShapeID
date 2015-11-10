%% load a fine mesh
load('samples/PNAS/D09.mat');

%% get subsampled mesh
subV = G.GeodesicFarthestPointSampling(400);
[subV,subF] = G.PerformGeodesicDelauneyTriangulation(subV,[]);
Gs = Mesh('VF',subV,subF);

%% compute uniformization for the subsampled mesh
[Gs.Aux.Area,Gs.Aux.Center] = Gs.Centralize('ScaleArea');
[~,TriAreas] = Gs.ComputeSurfaceArea;
Gs.Aux.VertArea = (TriAreas'*Gs.F2V)/3;
D = Gs.ComputeCPMS().ComputeUniformization(struct('method','LSCM','boundary_conditions','disc'));
Gs.Aux.UniformizationV = D.V;

%% estimate conformal factor using area distortion
[~,AG] = Gs.ComputeSurfaceArea;
[~,DG] = D.ComputeSurfaceArea;
Gs.Aux.Conf = (AG./DG)'*D.F2V/3;
%%%%% set the G.Aux.Conf to 0 on the boundary
[Gs.BV,Gs.BE] = Gs.FindBoundaries();
Gs.Aux.Conf(Gs.BV) = 0;

%% keep features on the original fine mesh for the downsampled mesh
TREE = kdtree_build(Gs.V');
Gs.Aux.ConfMaxInds = kdtree_nearest_neighbor(TREE, G.V(:,G.Aux.ConfMaxInds)');

%%% sanity check
Gs.draw();
hold on
scatter3(Gs.V(1,Gs.Aux.ConfMaxInds),Gs.V(2,Gs.Aux.ConfMaxInds),...
    Gs.V(3,Gs.Aux.ConfMaxInds),20,'g','filled');

%% similarly, find nearest neighbors of
%%%% G.Aux.ADMaxInds, G.Aux.GaussMaxInds, G.Aux.GaussMinInds
Gs.Aux.ADMaxInds = kdtree_nearest_neighbor(TREE, G.V(:,G.Aux.ADMaxInds)');
Gs.Aux.GaussMaxInds = kdtree_nearest_neighbor(TREE, G.V(:,G.Aux.GaussMaxInds)');
Gs.Aux.GaussMinInds = kdtree_nearest_neighbor(TREE, G.V(:,G.Aux.GaussMinInds)');

%%% sanity check
Gs.draw();
hold on
scatter3(Gs.V(1,Gs.Aux.ADMaxInds),Gs.V(2,Gs.Aux.ADMaxInds),...
    Gs.V(3,Gs.Aux.ADMaxInds),20,'y','filled');

Gs.draw();
hold on
scatter3(Gs.V(1,Gs.Aux.GaussMaxInds),Gs.V(2,Gs.Aux.GaussMaxInds),...
    Gs.V(3,Gs.Aux.GaussMaxInds),20,'r','filled');
scatter3(Gs.V(1,Gs.Aux.GaussMinInds),Gs.V(2,Gs.Aux.GaussMinInds),...
    Gs.V(3,Gs.Aux.GaussMinInds),20,'b','filled');
