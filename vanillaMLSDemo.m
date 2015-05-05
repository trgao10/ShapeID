%% preparation
close all;
path(pathdef);
addpath(path,genpath([pwd '/utils/']));

%% specify mesh
Name = 'j08';

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

%% load mesh (called G)
load([sample_path taxa_code{strcmpi(taxa_code, Name)} '.mat']);

figure;
scatter3(G.V(1,:),G.V(2,:),G.V(3,:),5,'k','filled');
axis equal;
axis off;
hold on;
cameratoolbar;
cameratoolbar('SetCoordSys', 'none');

%%%%% 
%%% pick one vertex, displace along the direction of the normal
% vIdx = 2107;
% vIdx = 3077;
% vIdx = 3089;
% vIdx = 1539;
% vIdx = 269;
% vIdx = 932;
vIdx = 1817;
scatter3(G.V(1,vIdx),G.V(2,vIdx),G.V(3,vIdx),20,'g','filled');
testPt = G.V(:,vIdx)+0.02*G.Nv(:,vIdx);
scatter3(testPt(1),testPt(2),testPt(3),20,'r','filled');

%% estimate initial normal at testPt
tic;
projPt = vanillaMLS(testPt, G.V);
toc;
scatter3(projPt(1),projPt(2),projPt(3),20,'m','filled');

% A = sort(squareform(pdist(G.V'))+diag(Inf(G.nV,1)),2);
% sigma = mean(A(:,5)); %% use the mean distance of each vertex to its 5-th neighbor as bandwidth
% distTestPt2PCloud = pdist2(testPt',G.V');
% sd = sort(distTestPt2PCloud);
% iSigma = sd(5);
% Weights = exp(-distTestPt2PCloud.^2/iSigma^2);
% coeff = pca(G.V(:,Weights>1e-3)','Weights',Weights(Weights>1e-3));
% EstNormal = coeff(:,3);
% 
% weightedCenter = Weights*G.V'/sum(Weights);
% quiver3(weightedCenter(1),weightedCenter(2),weightedCenter(3),coeff(1,1),coeff(2,1),coeff(3,1),0.1,'LineWidth',2);
% quiver3(weightedCenter(1),weightedCenter(2),weightedCenter(3),coeff(1,2),coeff(2,2),coeff(3,2),0.1,'LineWidth',2);
% quiver3(weightedCenter(1),weightedCenter(2),weightedCenter(3),-coeff(1,3),-coeff(2,3),-coeff(3,3),0.1,'LineWidth',2);
% 
% Proj = @(r,t,normal) r+t*normal;
% sF = @(t) ((EstNormal'*(G.V-repmat(Proj(testPt,t,EstNormal),1,G.nV))).^2)*exp(-pdist2(Proj(testPt,t,EstNormal)',G.V').^2/iSigma^2)'/sum(exp(-pdist2(Proj(testPt,t,EstNormal)',G.V').^2/iSigma^2));
% tmin = fminbnd(sF,-2*iSigma,2*iSigma);
% ProjTestPt = Proj(testPt,tmin,EstNormal);
% 
% scatter3(ProjTestPt(1),ProjTestPt(2),ProjTestPt(3),20,'b','filled');
% 
% F = @(q) ((EstNormal'*(G.V-repmat(q,1,G.nV))).^2)*(exp(-pdist2(q',G.V').^2/iSigma^2)'/sum(exp(-pdist2(q',G.V').^2/iSigma^2)));
% [q,fval] = fmincon(F,ProjTestPt,[],[],EstNormal',ProjTestPt'*EstNormal);
% 
% scatter3(q(1),q(2),q(3),20,'m','filled');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% test project for sample vertices
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
numSamples = 100;
projVs = zeros(3,numSamples);

cback = 0;
for j=1:numSamples
    testPt = G.V(:,G.Aux.DensityPnts(j))+0.01*G.Nv(:,G.Aux.DensityPnts(j));
    projVs(:,j) = vanillaMLS(testPt, G.V);
    for cc=1:cback
        fprintf('\b');
    end
    cback = fprintf(['%4d/' num2str(numSamples) ' done.\n'], j);
end

figure;
scatter3(G.V(1,G.Aux.DensityPnts(1:numSamples)),...
    G.V(2,G.Aux.DensityPnts(1:numSamples)),...
    G.V(3,G.Aux.DensityPnts(1:numSamples)),5,'k','filled');
axis equal;
axis off;
hold on;
cameratoolbar;
cameratoolbar('SetCoordSys', 'none');
scatter3(projVs(1,:),projVs(2,:),projVs(3,:),5,'r','filled');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% test project for all vertices
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% projVs = zeros(size(G.V));
% 
% cback = 0;
% for j=1:G.nV
%     testPt = G.V(:,j)+0.01*G.Nv(:,j);
%     projVs(:,j) = vanillaMLS(testPt, G.V);
%     for cc=1:cback
%         fprintf('\b');
%     end
%     cback = fprintf(['%4d/' num2str(G.nV) ' done.\n'], j);
% end
% 
% figure;
% scatter3(G.V(1,:),G.V(2,:),G.V(3,:),5,'k','filled');
% axis equal;
% axis off;
% hold on;
% cameratoolbar;
% cameratoolbar('SetCoordSys', 'none');
% scatter3(projVs(1,:),projVs(2,:),projVs(3,:),5,'r','filled');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% gradient descent to find the initial t
% t = 0;
% stepSize = iSigma/10;
% F = @(t,normal,s) ((normal'*(G.V-repmat(testPt+t*normal,1,G.nV))).^2)*exp(-sum((G.V-repmat(testPt+t*normal,1,G.nV)).^2)/s^2)';
% GradFt = @(t,normal,s) 2*((1+((normal'*(G.V-repmat(testPt+t*normal,1,G.nV))).^2)/s^2).*(normal'*(-G.V+repmat(testPt+t*normal,1,G.nV))))*exp(-sum((G.V-repmat(testPt+t*normal,1,G.nV)).^2)/s^2)';
% oldFValue = F(t,EstNormal,iSigma);
% while (true)
%     gradVec = GradFt(t,EstNormal,iSigma);
%     t = t-stepSize*gradVec/norm(gradVec);
%     newFValue = F(t,EstNormal,iSigma);
%     disp(newFValue);
%     if abs(newFValue-oldFValue)<1e-5
%         break;
%     end
%     oldFValue = newFValue;
% end

% %% shift all vertices by certain amount along normals
% testG = G.V+0.02*G.Nv;
% scatter3(testG(1,:),testG(2,:),testG(3,:),5,'r','filled');

% %% experiment MLS3D
% scale = 3;
% % DETERMINE RADIUS OF SUPPORT OF EVERY NODE
% dmI = scale *0.5* ones(4, G.nV);
% [PHI, DPHIx, DPHIy,DPHIz] = MLS3DShape(1,G.nV,G.V(1,:),G.V(2,:),G.V(3,:),...
%     G.nV,testG(1,:),testG(2,:),testG(3,:),dmI,'CUBIC',0.03);
% Vx = PHI*testG(1,:)';
% Vy = PHI*testG(2,:)';
% Vz = PHI*testG(3,:)';

%%% vanillaMLS
% pt = testPt';
% oldProjPt = testPt';
% pCloud = G.V';
% oldErr = Inf;
% nIt = 0;
% dist2TestPt = pdist2(pt,pCloud);
% sd = sort(dist2TestPt);
% sigma = sd(5);
% while(true)
%     dist2TestPt = pdist2(oldProjPt,pCloud);
%     %     sigma = min(dist2TestPt);
%     %     sd = sort(dist2TestPt);
%     %     sigma = sd(3);
%     Weights = exp(-dist2TestPt.^2/sigma^2);
%     Weights = Weights/sum(Weights);
%     weightedCenter = Weights*pCloud;
%     %     weightedCenter = Weights*pCloud/sum(Weights);
%     
%     %     sW = sort(Weights);
%     %     Weights(Weights<sW(end-10)) = 0;
%     %     weightedCenter = Weights*pCloud/sum(Weights);
%     deMeanPCloud = pCloud-repmat(weightedCenter,size(pCloud,1),1);
%     coeff = pca(deMeanPCloud(Weights>1e-3,:),'Weights',Weights(Weights>1e-3));
%     %     deMeanPCloud = deMeanPCloud(Weights>0,:);
%     %     deMeanWeights = Weights(Weights>0);
%     %     coeff = pca(deMeanPCloud,'Weights',deMeanWeights);
%     
%     projPt = (pt-weightedCenter)*coeff(:,1:2)*coeff(:,1:2)'+weightedCenter;
%     
%     %     scatter3(weightedCenter(1),weightedCenter(2),weightedCenter(3),20,'m','filled');
%     %     scatter3(projPt(1),projPt(2),projPt(3),20,'b','filled');
%     
%     err = norm(projPt-oldProjPt);
%     disp(err);
%     if ((norm(projPt-oldProjPt) < 1e-6) || abs(err-oldErr) < 1e-6)
%         break;
%     end
%     oldProjPt = projPt;
%     oldErr = err;
%     nIt = nIt+1;
% end
% 
% scatter3(weightedCenter(1),weightedCenter(2),weightedCenter(3),20,'m','filled');
% scatter3(projPt(1),projPt(2),projPt(3),20,'b','filled');
% quiver3(weightedCenter(1),weightedCenter(2),weightedCenter(3),coeff(1,1),coeff(2,1),coeff(3,1),0.1,'LineWidth',2);
% quiver3(weightedCenter(1),weightedCenter(2),weightedCenter(3),coeff(1,2),coeff(2,2),coeff(3,2),0.1,'LineWidth',2);
% % Weights = exp(-pdist2(testPt',G.V').^2/min(pdist2(testPt',G.V')));
% % [coeff,score,latent] = pca(G.V','Weights',Weights);
% %
% % projTestPt = coeff(:,1:2)*coeff(:,1:2)'*(testPt-G.V(:,vIdx))+G.V(:,vIdx);
% % projTestPt = vanillaMLS(testPt, G.V);
% 
% % quiver3(G.V(1,vIdx),G.V(2,vIdx),G.V(3,vIdx),coeff(1,1),coeff(2,1),coeff(3,1),0.1,'LineWidth',2);
% % quiver3(G.V(1,vIdx),G.V(2,vIdx),G.V(3,vIdx),coeff(1,2),coeff(2,2),coeff(3,2),0.1,'LineWidth',2);
% % quiver3(G.V(1,vIdx),G.V(2,vIdx),G.V(3,vIdx),coeff(1,3),coeff(2,3),coeff(3,3),0.1,'LineWidth',2);
% % scatter3(projTestPt(1),projTestPt(2),projTestPt(3),20,'b','filled');
% % line([testPt(1),projTestPt(1)],[testPt(2),projTestPt(2)],[testPt(3),projTestPt(3)]);






