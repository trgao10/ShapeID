%% preparation
% close all;
path(pathdef);
addpath(path,genpath([pwd '/utils/']));

%% setup parameter
Names = {'x02', 'Q19'};
GroupSize = length(Names);

%% setup paths
base_path = [pwd '/'];
data_path = '../DATA/PNAS/';
spreadsheet_path = [data_path 'ClassificationTable.xlsx'];
sample_path = '../cPdist/samples/PNAS/';
LandmarksPath = [data_path 'landmarks_teeth.mat'];
MeshesPath = [data_path 'meshes/'];
MeshSuffix = '_sas.off';

%% load taxa codes
taxa_file = [data_path 'teeth_taxa_table.mat'];
taxa_code = load(taxa_file);
taxa_code = taxa_code.taxa_code;

%% load meshes
MeshList = cell(GroupSize,1);
for j=1:GroupSize
    load([sample_path taxa_code{strcmpi(taxa_code, Names{j})} '.mat']);
    G.Aux.ObLmk = GetLandmarks(G.Aux.name, LandmarksPath, [MeshesPath G.Aux.name MeshSuffix]);
    MeshList{j} = G;
end

originalMeshes = drawMeshList(MeshList, []);
set(originalMeshes, 'Name', 'Orignal Meshes');

%% reparametrize all meshes

%%%% Step 1: Get all conformal parametrizations
R = cell(GroupSize,1);
R{1} = eye(3);

FlatMeshList = cell(GroupSize,1);
FlatMeshOrientation = zeros(GroupSize,1); %% '1' means parametrization flipped
FlatMeshList{1} = Mesh('VF', MeshList{1}.Aux.UniformizationV, MeshList{1}.F);
for j=2:GroupSize
    FlatVertices = MeshList{j}.Aux.UniformizationV;

    %%% Procrustes matching observer landmarks
    [U,~,V] = svd(FlatVertices(:,MeshList{j}.Aux.ObLmk)*FlatMeshList{1}.V(:,MeshList{1}.Aux.ObLmk)');
    R{j} = V*U';
    
    if (det(R{j}(1:2,1:2)) < 0)
        FlatMeshOrientation(j) = 1;
        disp(['Mesh ' num2str(j) ' flipped!']);
    end

    FlatVertices = R{j}*FlatVertices;
    FlatMeshList{j} = Mesh('VF', FlatVertices, MeshList{j}.F);
end

%%% draw Procrustes aligned meshes
drawMeshList(FlatMeshList, []);

%%%% Step 2: TPS 2 to the end meshes to the first mesh
DeformedMeshList = cell(GroupSize,1);
DeformedMeshList{1} = FlatMeshList{1};
TPS_FEATURESN = DISCtoPLANE(FlatMeshList{1}.V(1:2, MeshList{1}.Aux.ObLmk)','d2p');
for j=2:GroupSize
    TPS_FEATURESM = DISCtoPLANE(FlatMeshList{j}.V(1:2,MeshList{j}.Aux.ObLmk)','d2p');
    tP = DISCtoPLANE(FlatMeshList{j}.V(1:2,:)','d2p');
    [ftps] = TEETH_calc_tps(TPS_FEATURESM,TPS_FEATURESN-TPS_FEATURESM);
    pt = tP + TEETH_eval_tps(ftps,tP);
    DeformedMeshList{j} = Mesh('VF', DISCtoPLANE(pt,'p2d')', FlatMeshList{j}.F);
    DeformedMeshList{j}.V(:,isnan(compl(DeformedMeshList{j}.V))) = 1;
end

%%% draw TPS deformed meshes
drawMeshList(DeformedMeshList, []);

%%%% Step 3: Get standard unit disk mesh
% %% GMM estimate high density locations
% % initialMu = DeformedMeshList{1}.V(:,MeshList{1}.Aux.ConfMaxInds);
% samplePts = []';
% for j=1:GroupSize
%     candIdx = find(MeshList{j}.Aux.Conf>mean(MeshList{j}.Aux.Conf));
%     samplePts = [samplePts;DeformedMeshList{j}.V(1:2,candIdx)'];
% end
% 
% samplePts(arrayfun(@(idx) norm(samplePts(idx,:)),1:size(samplePts,1))>0.95,:) = [];
% 
% figure;
% scatter(samplePts(:,1),samplePts(:,2),10,'g','filled');
% axis equal;hold on;
% 
% [idx,C] = kmeans(samplePts,4);
% scatter(C(:,1),C(:,2),20,'b','filled');
% 
% GMModel = fitgmdist(samplePts,4,...
%     'Start',struct('mu',C,'Sigma',repmat(0.01*eye(2),1,1,4)),...
%     'options',statset('MaxIter',400));
% scatter(GMModel.mu(:,1),GMModel.mu(:,2),20,'r','filled');
% ezcontour(@(x1,x2)pdf(GMModel,[x1 x2]),get(gca,{'XLim','YLim'}));
% 
% %%%%% generate mesh from distmesh
% n = 300;
% [X,Y] = meshgrid(-1:(1/n):1);
% V = griddata(DeformedMeshList{1}.V(1,:),DeformedMeshList{1}.V(2,:),MeshList{1}.Aux.Conf,X,Y);
% V(isnan(V))=0;
% fd = @(p) sqrt(sum(p.^2,2))-0.9^2; %%% choose a smaller disk
% Gaussian = @(p,mu,sigma) (exp(-sum((p-repmat(mu,size(p,1),1)).*(sigma\(p-repmat(mu,size(p,1),1))')',2))/2);
% % Gaussian = @(p,mu,sigma) (1/sqrt(2*pi*det(sigma)))*(exp(-sum((p-repmat(mu,size(p,1),1)).*(sigma\(p-repmat(mu,size(p,1),1))')',2))/2);
% fh = @(p) 50*(3.1-10*Gaussian(p,GMModel.mu(1,:),GMModel.Sigma(:,:,1)*20))+...
%     50*(3.1-10*Gaussian(p,GMModel.mu(2,:),GMModel.Sigma(:,:,2)*20))+...
%     50*(3.1-10*Gaussian(p,GMModel.mu(3,:),GMModel.Sigma(:,:,3)*20))+...
%     50*(3.1-10*Gaussian(p,GMModel.mu(4,:),GMModel.Sigma(:,:,4)*20));
% % fh = @(p) 50*(101-Gaussian(p,GMModel.mu(1,:),GMModel.Sigma(:,:,1)));
% [UDVs, UDFs] = distmesh2d(fd, fh, 0.015, [-1,-1;1,1], GMModel.mu);
% UnitDistMesh = Mesh('VF', UDVs', UDFs');
% figure;UnitDistMesh.draw();

%%%%% generate domain mesh from a fine distmesh
% fd = @(p) sqrt(sum(p.^2,2))-0.9^2; %%% choose a smaller disk
% [UDVs, UDFs] = distmesh2d(fd, fh, 0.015, [-1,-1;1,1], GMModel.mu);
% UnitDistMesh = Mesh('VF', UDVs', UDFs');
% figure;UnitDistMesh.draw();

%%%%% generate domain mesh from a random mesh
domainIdx = floor(rand()*GroupSize)+1;
UDVs = DeformedMeshList{domainIdx}.V(1:2,:)';
UDVs(sum(UDVs.^2,2)>0.9^2,:) = [];
DT = delaunayTriangulation(UDVs);
UDFs = DT.ConnectivityList;
UnitDistMesh = Mesh('VF', UDVs', UDFs');
figure;UnitDistMesh.draw();

%%
%%%% Step 4: Re-parametrize all meshes
% BBox = [-1.1, -1.1, 1.1, 1.1; -1.1, 1.1, -1.1, 1.1];
ReparametrizedMeshList = cell(GroupSize,1);
for j=1:GroupSize
    disp(['Barycentric Interpolation for Mesh ' num2str(j) '...']);
    NewMeshVertices = zeros(3, size(UDVs,1));
    
    TR = triangulation(DeformedMeshList{j}.F', DeformedMeshList{j}.V(1:2,:)');
    
    cback = 0;
    for k=1:size(UDVs,1)
        BC = TR.cartesianToBarycentric((1:DeformedMeshList{j}.nF)',repmat(UDVs(k,:),DeformedMeshList{j}.nF,1));
        
        tind = find(all(BC>-1e-10,2));
        if(numel(tind)>1)
            tind = tind(1);
        elseif numel(tind)==0
            warning(['For j = ' num2str(j) ', Point ', num2str(k), ' was not found in any triangle.']);
            continue;
        end
        BC = BC(tind,:);
        
        NewMeshVertices(:, k) = MeshList{j}.V(:,TR.ConnectivityList(tind,:))*BC';
        
        for cc=1:cback
            fprintf('\b');
        end
        cback = fprintf(['%4d/' num2str(size(UDVs,1)) ' done.\n'], k);
    end
    
    ReparametrizedMeshList{j} = Mesh('VF', NewMeshVertices, UDFs');
end

%%% draw reparametrized meshes
reparametrizedMeshes = drawMeshList(ReparametrizedMeshList, []);
set(reparametrizedMeshes, 'Name', 'Reparametrized Meshes');

%% write reparametrized meshes to .off files
touch('./meshes/reparametrized/');
for j=1:GroupSize
    ReparametrizedMeshList{j}.Write(['./meshes/reparametrized/' MeshList{j}.Aux.name '.off'], 'off', []);
end

UnitDistMesh.Write('./meshes/reparametrized/domainMesh.off', 'off', []);
