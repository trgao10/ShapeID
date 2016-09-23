%% preparation
close all;
path(pathdef);
addpath(path,genpath([pwd '/utils/']));
samples_path = '/media/trgao10/Work/MATLAB/cPdist/samples/PNAS/';

Name = 'ah03';
numEigVals = 200;
EigRange = 1:200;

load([samples_path Name '.mat']);

eigopt = struct('isreal',1,'issym',1,'maxit',5000,'disp',0);
[B,L] = eigs(G.Aux.LB, numEigVals, 'SM', eigopt);

% G.ViewFunctionOnMesh(V(:,2),[]);

reducedV = G.V*(B(:,EigRange)*B(:,EigRange)');
reducedG = Mesh('VF', reducedV, G.F);

reducedG.draw();


