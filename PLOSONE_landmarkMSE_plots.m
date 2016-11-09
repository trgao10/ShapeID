% load('/xtmp/ArchivedResults/PNAS/cPMST/FeatureFixOff/cPMSTlmkMSEMatrix.mat');
% 
% yCoords = nonzeros(triu(lmkMSEMatrix,1));
% 
% figure;
% plot(xCoords, yCoords, 'b.');
% axis equal
% axis square
% hold on
% line([0,1],[0,1],'Color','r');

clearvars;
close all;

load('/xtmp/ArchivedResults/PNAS/cPdist/cPlmkMSEMatrix.mat');
xCoords = nonzeros(triu(cPlmkMSEMatrix,1));

%%%% prepare paths
result_path = '/xtmp/ArchivedResults/PNAS/';

%%%% different improvements
MethodsType = {'cPComposedLAST', 'cPComposedLAST', 'cPComposedLAST',...
    'cPComposedLAST', 'cPLAST', 'cPLAST', 'cPMST', 'cPViterbi',...
    'cPViterbi', 'cPViterbi'};
Methods = {'cPComposedLAST', 'cPComposedLASTbalance', 'cPComposedLASTmean',...
    'cPComposedLASTmedian', 'cPLAST', 'cPLASTbalance', 'cPMST', 'cPViterbi',...
    'cPViterbiAngle0.5', 'cPViterbiAngle0.25'};
FeatureFix = {'FeatureFixOff', 'FeatureFixOn'};

figPath = '~/Desktop/PLOSONE/';
if ~exist(figPath, 'dir')
    mkdir(figPath);
    mkdir([figPath 'fig/']);
    mkdir([figPath 'eps/']);
end

figure('Position', [1000,774,740,724]);
for MethodsIdx = 1:length(Methods)
    for FeatureFixIdx = 1:length(FeatureFix)
        lmkMSEMatrixPath = [result_path Methods{MethodsIdx} filesep FeatureFix{FeatureFixIdx} filesep];
        disp(lmkMSEMatrixPath);
        load([lmkMSEMatrixPath MethodsType{MethodsIdx} 'lmkMSEMatrix']);
        yCoords = nonzeros(triu(lmkMSEMatrix,1));
        clf;
        plot(xCoords, yCoords, 'b.');
        axis equal
        axis square
        hold on
        line([0,1],[0,1],'Color','r','LineWidth',2);
        xlabel('Landmark MSEs for cPDist');
        if strcmpi(FeatureFix(FeatureFixIdx), 'FeatureFixOn')
            ylabel(['Landmark MSEs for ' Methods{MethodsIdx} ' (with feature-fixing)']);
        else
            ylabel(['Landmark MSEs for ' Methods{MethodsIdx} ' (without feature-fixing)']);
        end
%         pause()
        savefig([figPath 'fig/' strrep(Methods{MethodsIdx}, '.', '_') '_' FeatureFix{FeatureFixIdx} '.fig']);
        print([figPath 'eps/' strrep(Methods{MethodsIdx}, '.', '_') '_' FeatureFix{FeatureFixIdx} '.eps'],'-depsc');
    end
end