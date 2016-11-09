clear vas
close all
path(pathdef);
addpath(path,genpath([pwd '/utils/']));
addpath(path,genpath('~/Documents/MATLAB/iso2mesh/'));

dicomImgFolder = './data/AMNH-M-203297_M515-437/';
% dicomImgFolder = '/home/trgao10/Downloads/MCZ-38316/';
dicomImgFiles = dir(dicomImgFolder);
numImgFiles = length(dicomImgFiles)-2;
dicomImgFilePaths = cell(1, numImgFiles);
for i = 1:numImgFiles
    dicomImgFilePaths{i} = [dicomImgFolder dicomImgFiles(i+2).name];
end

%%% get image size
info = dicominfo(dicomImgFilePaths{1});
X = repmat(uint16(0), [info.Rows info.Columns 1 numImgFiles]);
for j = 1:numImgFiles
    X(:,:,1,j) = dicomread(dicomImgFilePaths{j});
end

X = X / max(X(:));  %%% entries of X are either 0 or 127, funny enough....
montage(X,[0,1]);

keyboard

%%% apply iso2mesh to generate volumetric mesh
squeezeX = squeeze(X);
[node,elem,face] = vol2mesh(squeezeX,1:size(squeezeX,1),1:size(squeezeX,2),1:size(squeezeX,3),2,2,1);
plotmesh(node,face);
axis equal;

%%% can also apply iso2mesh to generate surfaces
[V,F,regions,holes] = vol2surf(squeezeX,1:size(squeezeX,1),1:size(squeezeX,2),1:size(squeezeX,3),2,1);
trisurf(F(:,1:3),V(:,1),V(:,2),V(:,3));
axis equal;


