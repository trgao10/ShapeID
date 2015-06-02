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

%% collection alignment rigid motions
load([data_path 'cPMSTinitRms.mat']);

%% interpolate each segment
PoissonMSTPath = './meshes/PoissonMSTTeeth/';
touch(PoissonMSTPath);

for j=1:length(PRED)
    if (PRED(j) == 0)
        continue;
    end
    linkPath = [PoissonMSTPath taxa_code{PRED(j)} '_' taxa_code{j} '/'];
    if exist(linkPath, 'dir')
        continue;
    end
    touch(linkPath);
    touch([linkPath 'off']);
    touch([linkPath 'ObLmk']);
    
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
    
    [~,ObCoords] = GetLandmarks(taxa_code{strcmpi(taxa_code, MeshList{1}.Aux.name)},...
        LandmarksPath, [MeshesPath taxa_code{strcmpi(taxa_code, MeshList{1}.Aux.name)} MeshSuffix]);
    Vtree = kdtree_build( ReparametrizedMeshList{1}.V' );
    ObInds = kdtree_nearest_neighbor(Vtree, ObCoords);
    kdtree_delete(Vtree);
    
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
    
%     drawMeshList(reconMeshList, struct('DisplayLayout', [2,5], 'linkCamera', 'on'));
    
    %%%% Step 7: write interpolated meshes to .off files
    for k=1:length(Weights)
        reconMeshList{k}.Write([linkPath 'off/' taxa_code{PRED(j)} '_' taxa_code{j} '_' sprintf('%02d', k) '.off'], 'off', []);
        csvwrite([linkPath 'ObLmk/' taxa_code{PRED(j)} '_' taxa_code{j} '_' sprintf('%02d', k) '.csv'],reconMeshList{k}.V(:,ObInds)');
    end
    save([linkPath 'reparMeshes.mat'], 'MeshList', 'ReparametrizedMeshList', 'domainMesh');
    save([linkPath 'reconMeshes.mat'], 'reconMeshList', 'Weights');
end

