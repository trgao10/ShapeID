%% preparation
close all;
path(pathdef);
addpath(path,genpath([pwd '/utils/']));

%% setup parameter
Names = {'Q18'};
GroupSize = length(Names);

%% setup paths
base_path = [pwd '/'];
data_path = '../DATA/PNAS/';
spreadsheet_path = [data_path 'ClassificationTable.xlsx'];
reparametrized_path = './meshes/reparametrized/';
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
    G = Mesh('off', [reparametrized_path taxa_code{strcmpi(taxa_code, Names{j})} '.off']);
    G.Aux.ObLmk = GetLandmarks(taxa_code{strcmpi(taxa_code, Names{j})},...
        LandmarksPath, [MeshesPath taxa_code{strcmpi(taxa_code, Names{j})} MeshSuffix]);
    MeshList{j} = G;
end

reparametrizedMeshes = drawMeshList(MeshList, []);
set(reparametrizedMeshes, 'Name', 'Reparametrized Meshes');

%% extract affine transformations
domainMesh = Mesh('off', [reparametrized_path 'domainMesh.off']);

% AffineTransformations = zeros(domainMesh.nF, GroupSize, 3, 3);
AffineTransformations = zeros(domainMesh.nF, GroupSize, 3, 2);
ORT = zeros(domainMesh.nF, GroupSize, 3, 3);
PSD = zeros(domainMesh.nF, GroupSize, 3, 3);

cback = 0;
for j=1:domainMesh.nF
    localF = domainMesh.F(:,j);
    localV = domainMesh.V(:,localF);
    DomainLocalFrame = [localV(1:2,2)-localV(1:2,1), localV(1:2,3)-localV(1:2,1)];
%     DomainLocalFrame = [localV(:,2)-localV(:,1), localV(:,3)-localV(:,1), [0,0,1]'];
    for k=1:GroupSize
        localFk = MeshList{k}.F(:,j);
        localVk = MeshList{k}.V(:,localFk);
        kLocalFrame = [localVk(:,2)-localVk(:,1), localVk(:,3)-localVk(:,1)];
%         kLocalFrame = [localVk(:,2)-localVk(:,1), localVk(:,3)-localVk(:,1), [0,0,0]'];
%         kLocalFrame(:,3) = cross(kLocalFrame(:,1),kLocalFrame(:,2));
%         kLocalFrame(:,3) = kLocalFrame(:,3)/norm(kLocalFrame(:,3));
        
        AffineTransformations(j,k,:,:) = kLocalFrame/DomainLocalFrame;
        [U,S,V] = svd(squeeze(AffineTransformations(j,k,:,:)));
        S = [S,zeros(3,1)];
        makeRot = det(U)*det(V);
        if (makeRot > 0)
            V = [[V;zeros(1,2)],[0,0,1]'];
        else
            V = [[V;zeros(1,2)],[0,0,-1]'];
        end
        ORT(j,k,:,:) = U*V';
        PSD(j,k,:,:) = V*S*V';
    end
    
    for cc=1:cback
        fprintf('\b');
    end
    cback = fprintf(['%4d/' num2str(domainMesh.nF) ' done.\n'], j);
end

%% reconstruction
%%%% interpolate affine transformation using polar decomposition
Weights = [1]; %% to be renormalized later
Weights = Weights/sum(Weights);

InterpAffineTransformations = zeros(domainMesh.nF, 3, 3);

cback = 0;
for j=1:domainMesh.nF
    InterpPSD = zeros(3,3);
    InterpQua = zeros(4,1);
    for k=1:GroupSize
        InterpPSD = InterpPSD + Weights(k)*squeeze(PSD(j,k,:,:));
        InterpQua = InterpQua + Weights(k)*qGetQ(squeeze(ORT(j,k,:,:)));
    end
    InterpAffineTransformations(j, :, :) = qGetR(InterpQua/norm(InterpQua))*InterpPSD;
    
    for cc=1:cback
        fprintf('\b');
    end
    cback = fprintf(['%4d/' num2str(domainMesh.nF) ' done.\n'], j);
end

%%% sanity check
maxD = 0;
for j=1:domainMesh.nF
    D = norm(squeeze(AffineTransformations(j,:,:,:))-squeeze(InterpAffineTransformations(j,:,1:2)));
    if (D > maxD)
        maxD = D;
    end
end
disp(['Max Discrepancy = ' num2str(maxD)]);

%%%% Poisson reconstruction
%%%%%% IDEA: \phi is a map from domainMesh to \mathbb{R}^3, which stands
%%%%%% for the embedding; \nabla\phi is a linear transformation of \phi
%%%%%% (depending on the geometry of domainMesh), and we want to minimize
%%%%%% the discrepancy between \nabla\phi and the interpolated Jacobians

%%% sanity check: what is the gradient of the default embedding?
[Gx,Gy] = domainMesh.ComputeGradientMatrix;
GradX = Gx*MeshList{1}.V';
GradY = Gy*MeshList{1}.V';

maxD = 0;
maxIdx = 1;
count = 0;
for j=1:domainMesh.nF
    GradFace = [GradX(j,:);GradY(j,:)]';
    D = norm(squeeze(AffineTransformations(j,:,:,:))-GradFace);
    if (D > maxD)
        maxD = D;
        maxIdx = j;
    end
end
disp(['Max GradFace Discrepancy = ' num2str(maxD)]);

%%% compute divergence of the interpolated jacobian
[DivX,DivY] = domainMesh.ComputeDivergenceMatrix;
DivRhoX = DivX*squeeze(InterpAffineTransformations(:,1,1))...
    +DivY*squeeze(InterpAffineTransformations(:,1,2));
DivRhoY = DivX*squeeze(InterpAffineTransformations(:,2,1))...
    +DivY*squeeze(InterpAffineTransformations(:,2,2));
DivRhoZ = DivX*squeeze(InterpAffineTransformations(:,3,1))...
    +DivY*squeeze(InterpAffineTransformations(:,3,2));

%%% reconstruct
L = domainMesh.ComputeCotanLaplacian;
reconX = [L;ones(1,domainMesh.nV)]\[DivRhoX;0];
reconY = [L;ones(1,domainMesh.nV)]\[DivRhoY;0];
reconZ = [L;ones(1,domainMesh.nV)]\[DivRhoZ;0];

%%% check reconstructed mesh
reconMesh = Mesh('VF',[reconX';reconY';reconZ'],domainMesh.F);
figure;
reconMesh.draw();
set(gcf,'Name','Reconstructed Mesh');



