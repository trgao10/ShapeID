%% preparation
close all;
path(pathdef);
addpath(path,genpath([pwd '/utils/']));

%% setup parameter
GroupLevel = 'Genus';
numRandMeshes = 10;

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
GroupSize = length(taxa_code);
ChunkSize = 55; %% PNAS

%% load CPD matrix and extract MST
cPD = load([data_path 'PNAS_sym_cPdists.mat']);
cPD = cPD.symmetrized_all_cPdistances;
cPD = sparse(cPD);
cPD = tril(cPD, -1);
[ST, PRED] = graphminspantree(cPD, 'Method', 'Kruskal');

%% set fixed weight
Weights = linspace(0,1,12);
Weights([1,length(Weights)]) = [];

%% collection rigid motions
load([data_path 'cPMSTinitRms.mat']);

%% interpolate each segment
PoissonMSTPath = './meshes/PoissonMST/';
touch(PoissonMSTPath);

for j=1:length(PRED)
    if (PRED(j) == 0)
        continue;
    end
    
    %%%% Step 1: load meshes and align them to each other
    G = load([sample_path taxa_code{j} '.mat']); MeshList{1} = G.G;
    G = load([sample_path taxa_code{PRED(j)} '.mat']); MeshList{2} = G.G;
    MeshList{1}.Aux.ObLmk = GetLandmarks(MeshList{1}.Aux.name, LandmarksPath, [MeshesPath MeshList{1}.Aux.name MeshSuffix]);
    MeshList{2}.Aux.ObLmk = GetLandmarks(MeshList{2}.Aux.name, LandmarksPath, [MeshesPath MeshList{2}.Aux.name MeshSuffix]);
    MeshList{2}.V = R{j,PRED(j)}*MeshList{2}.V;
    
    %%%% Step 2: Procrustes matching observer landmarks
    FlatMeshList = cell(2,1);
    FlatMeshList{1} = Mesh('VF', MeshList{1}.Aux.UniformizationV, MeshList{1}.F);

    FlatVertices = MeshList{2}.Aux.UniformizationV;
    [U,~,V] = svd(FlatVertices(:,MeshList{2}.Aux.ObLmk)*FlatMeshList{1}.V(:,MeshList{1}.Aux.ObLmk)');
    RJ = V*U';
    FlatVertices = RJ*FlatVertices;
    FlatMeshList{2} = Mesh('VF', FlatVertices, MeshList{2}.F);
    
    %%%% Step 3: TPS G2 --> G1
    DeformedMeshList = cell(2,1);
    DeformedMeshList{1} = FlatMeshList{1};
    TPS_FEATURESN = DISCtoPLANE(FlatMeshList{1}.V(1:2, MeshList{1}.Aux.ObLmk)','d2p');
    TPS_FEATURESM = DISCtoPLANE(FlatMeshList{2}.V(1:2, MeshList{2}.Aux.ObLmk)','d2p');
    tP = DISCtoPLANE(FlatMeshList{2}.V(1:2,:)','d2p');
    [ftps] = TEETH_calc_tps(TPS_FEATURESM,TPS_FEATURESN-TPS_FEATURESM);
    pt = tP + TEETH_eval_tps(ftps,tP);
    DeformedMeshList{2} = Mesh('VF', DISCtoPLANE(pt,'p2d')', FlatMeshList{2}.F);
    DeformedMeshList{2}.V(:,isnan(compl(DeformedMeshList{2}.V))) = 1;
    
    %%%% Step 4: generate domain mesh
    domainIdx = floor(rand()*2)+1;
    UDVs = DeformedMeshList{domainIdx}.V(1:2,:)';
    UDVs(sum(UDVs.^2,2)>0.9^2,:) = [];
    DT = delaunayTriangulation(UDVs);
    UDFs = DT.ConnectivityList;
    domainMesh = Mesh('VF', UDVs', UDFs');
    
    %%%% Step 5: re-parametrize all meshes
    ReparametrizedMeshList = cell(2,1);
    for t=1:2
        disp(['Barycentric Interpolation for Mesh ' num2str(t) '...']);
        NewMeshVertices = zeros(3, size(UDVs,1));
        
        TR = triangulation(DeformedMeshList{t}.F', DeformedMeshList{t}.V(1:2,:)');
        
        cback = 0;
        for k=1:size(UDVs,1)
            BC = TR.cartesianToBarycentric((1:DeformedMeshList{t}.nF)',repmat(UDVs(k,:),DeformedMeshList{t}.nF,1));
            
            tind = find(all(BC>-1e-10,2));
            if(numel(tind)>1)
                tind = tind(1);
            elseif numel(tind)==0
                warning(['For t = ' num2str(t) ', Point ', num2str(k), ' was not found in any triangle.']);
                continue;
            end
            BC = BC(tind,:);
            
            NewMeshVertices(:, k) = MeshList{t}.V(:,TR.ConnectivityList(tind,:))*BC';
            
            for cc=1:cback
                fprintf('\b');
            end
            cback = fprintf(['%4d/' num2str(size(UDVs,1)) ' done.\n'], k);
        end
        ReparametrizedMeshList{t} = Mesh('VF', NewMeshVertices, UDFs');
    end
    
    %%%% Step 6: generate intermediate meshes
    reconMeshList = cell(length(Weights),1);
    for k=1:length(Weights)
        if (k==1)
            [reconMeshList{k},ORT,PSD] = PoissonMeshInterpolation(domainMesh, ReparametrizedMeshList, [Weights(k),1-Weights(k)]);
        else
            reconMeshList{k} = PoissonMeshInterpolation(domainMesh, ReparametrizedMeshList, [Weights(k),1-Weights(k)], ORT, PSD);
        end
        disp([num2str(k) '/' num2str(length(Weights)) ' done.']);
    end
    
    drawMeshList(reconMeshList, struct('DisplayLayout', [2,5], 'linkCamera', 'on'));
    
    %%%% Step 7: write interpolated meshes to .off files
    linkPath = [PoissonMSTPath taxa_code{PRED(j)} '_' taxa_code{j} '/'];
    touch(linkPath);
    touch([linkPath 'off']);
    
    for k=1:length(Weights)
        reconMeshList{k}.Write([linkPath 'off/' taxa_code{PRED(j)} '_' taxa_code{j} '_' sprintf('%02d', k) '.off'], 'off', []);
    end
    save([linkPath 'reparMeshes.mat'], 'MeshList', 'ReparametrizedMeshList', 'domainMesh');
    save([linkPath 'reconMeshes.mat'], 'reconMeshList', 'Weights');
    
    pause();
end
    
% %% parse GroupNames
% [~,ClTable,~] = xlsread(spreadsheet_path);
% 
% GroupNames = unique(ClTable(2:end,strcmpi(ClTable(1,:),GroupLevel)));
% GroupNames(strcmpi(GroupNames,'NA')) = [];
% 
% for j=1:length(GroupNames)
%     GrouPath = ['./meshes/PoissonTeeth/' GroupNames{j} '/'];
%     touch(GrouPath);
%     if (exist([GrouPath 'off/'],'dir'))
%         listing = dir([GrouPath 'off/']);
%         if (length(listing) == (numRandMeshes+2))
%             continue;
%         end
%     end
%     
%     Names = ClTable(strcmpi(ClTable(1:end,strcmpi(ClTable(1,:),GroupLevel)),GroupNames{j}),1);
%     GroupSize = length(Names);
%     
%     if (exist([GrouPath 'reparametrizedMeshes.mat'], 'file'))
%         load([GrouPath 'reparametrizedMeshes.mat']);
%     else
%         %%% Step 1: align all teeth with the first tooth in the group
%         MeshList = cell(GroupSize,1);
%         TAXAinds = cellfun(@(x) find(strcmpi(taxa_code,x)), Names, 'UniformOutput', false);
%         for k=1:GroupSize
%             load([sample_path taxa_code{strcmpi(taxa_code, Names{k})} '.mat']);
%             G.Aux.ObLmk = GetLandmarks(G.Aux.name, LandmarksPath, [MeshesPath G.Aux.name MeshSuffix]);
%             if (k>=2)
%                 G.V = R{TAXAinds{1},TAXAinds{k}}*G.V;
%             end
%             MeshList{k} = G;
%         end
%         
%         %%% Step 2: reparametrize all teeth in this group
%         RJ = cell(GroupSize,1);
%         RJ{1} = eye(3);
%         
%         FlatMeshList = cell(GroupSize,1);
%         FlatMeshOrientation = zeros(GroupSize,1); %% '1' means parametrization flipped
%         FlatMeshList{1} = Mesh('VF', MeshList{1}.Aux.UniformizationV, MeshList{1}.F);
%         for k=2:GroupSize
%             FlatVertices = MeshList{k}.Aux.UniformizationV;
%             
%             %%% Procrustes matching observer landmarks
%             [U,~,V] = svd(FlatVertices(:,MeshList{k}.Aux.ObLmk)*FlatMeshList{1}.V(:,MeshList{1}.Aux.ObLmk)');
%             RJ{k} = V*U';
%             
%             if (det(RJ{k}(1:2,1:2)) < 0)
%                 FlatMeshOrientation(k) = 1;
%                 disp(['Mesh ' num2str(k) ' flipped!']);
%             end
%             
%             FlatVertices = RJ{k}*FlatVertices;
%             FlatMeshList{k} = Mesh('VF', FlatVertices, MeshList{k}.F);
%         end
%         
%         %%%% Step 3: TPS (2:end) meshes to the first mesh
%         DeformedMeshList = cell(GroupSize,1);
%         DeformedMeshList{1} = FlatMeshList{1};
%         TPS_FEATURESN = DISCtoPLANE(FlatMeshList{1}.V(1:2, MeshList{1}.Aux.ObLmk)','d2p');
%         for k=2:GroupSize
%             TPS_FEATURESM = DISCtoPLANE(FlatMeshList{k}.V(1:2,MeshList{k}.Aux.ObLmk)','d2p');
%             tP = DISCtoPLANE(FlatMeshList{k}.V(1:2,:)','d2p');
%             [ftps] = TEETH_calc_tps(TPS_FEATURESM,TPS_FEATURESN-TPS_FEATURESM);
%             pt = tP + TEETH_eval_tps(ftps,tP);
%             DeformedMeshList{k} = Mesh('VF', DISCtoPLANE(pt,'p2d')', FlatMeshList{k}.F);
%             DeformedMeshList{k}.V(:,isnan(compl(DeformedMeshList{k}.V))) = 1;
%         end
%         
%         %%%% Step 4: generate domain mesh
%         domainIdx = floor(rand()*GroupSize)+1;
%         UDVs = DeformedMeshList{domainIdx}.V(1:2,:)';
%         UDVs(sum(UDVs.^2,2)>0.9^2,:) = [];
%         DT = delaunayTriangulation(UDVs);
%         UDFs = DT.ConnectivityList;
%         domainMesh = Mesh('VF', UDVs', UDFs');
%         
%         %%%% Step 5: re-parametrize all meshes
%         ReparametrizedMeshList = cell(GroupSize,1);
%         for t=1:GroupSize
%             disp(['Barycentric Interpolation for Mesh ' num2str(t) '...']);
%             NewMeshVertices = zeros(3, size(UDVs,1));
%             
%             TR = triangulation(DeformedMeshList{t}.F', DeformedMeshList{t}.V(1:2,:)');
%             
%             cback = 0;
%             for k=1:size(UDVs,1)
%                 BC = TR.cartesianToBarycentric((1:DeformedMeshList{t}.nF)',repmat(UDVs(k,:),DeformedMeshList{t}.nF,1));
%                 
%                 tind = find(all(BC>-1e-10,2));
%                 if(numel(tind)>1)
%                     tind = tind(1);
%                 elseif numel(tind)==0
%                     warning(['For t = ' num2str(t) ', Point ', num2str(k), ' was not found in any triangle.']);
%                     continue;
%                 end
%                 BC = BC(tind,:);
%                 
%                 NewMeshVertices(:, k) = MeshList{t}.V(:,TR.ConnectivityList(tind,:))*BC';
%                 
%                 for cc=1:cback
%                     fprintf('\b');
%                 end
%                 cback = fprintf(['%4d/' num2str(size(UDVs,1)) ' done.\n'], k);
%             end
%             ReparametrizedMeshList{t} = Mesh('VF', NewMeshVertices, UDFs');
%         end
%         save([GrouPath 'reparametrizedMeshes.mat'],'MeshList','ReparametrizedMeshList','domainMesh');
%     end
%     
%     %%%% Step 6: randomly generate 50 meshes for this species group
%     Weights = cell(numRandMeshes,1);
%     ExtremalMeshIndices = cell(numRandMeshes,1);
%     reconMeshList = cell(numRandMeshes,1);
%     for k=1:numRandMeshes
%         if (length(Names)>4)
%             ExtremalMeshIndices{k} = datasample(1:length(Names),4,'Replace',false);
%             Weights{k} = rand(4,1);
%             Weights{k} = -log(Weights{k});
%         else
%             ExtremalMeshIndices{k} = 1:length(Names);
%             Weights{k} = rand(size(Names));
%             Weights{k} = -log(Weights{k});
%         end
%         fullWeights = zeros(size(Names));
%         fullWeights(ExtremalMeshIndices{k}) = Weights{k};
%         if (k==1)
%             [reconMeshList{k},ORT,PSD] = PoissonMeshInterpolation(domainMesh, ReparametrizedMeshList, fullWeights);
%         else
%             reconMeshList{k} = PoissonMeshInterpolation(domainMesh, ReparametrizedMeshList, fullWeights, ORT, PSD);
%         end
%         disp([num2str(k) '/' num2str(numRandMeshes) ' done.']);
%     end
%     save([GrouPath 'reconMeshes.mat'], 'reconMeshList', 'Weights', 'ExtremalMeshIndices');
%     
% %     interpolatedMeshes = drawMeshList(reconMeshList, struct('DisplayLayout',[2,5],'linkCamera','on','DisplayOrient','Horizontal'));
% %     set(gcf,'Name','Randomly Interpolated Meshes');
% 
%     %%%% Step 7: write randomly interpolated meshes to .off files
%     touch([GrouPath 'off/']);
%     for k=1:numRandMeshes
%         reconMeshList{k}.Write([GrouPath 'off/' GroupNames{j} '_' sprintf('%02d', k) '.off'], 'off', []);
%     end
% end
% 
% 
