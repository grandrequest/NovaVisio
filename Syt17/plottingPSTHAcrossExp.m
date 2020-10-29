%import createLoader
import jauimodel.*
import vuidocument.*

%Data and export folder paths
dataFolder16 = "/Users/jnagy2/Documents/20200716";
exportName16 = "20200716RodsPSTH.mat";
dataFolder21 = "/Users/jnagy2/Documents/20200721";
exportName21 = "20200721RodsPSTH.mat";
%% Create Parameters for splitting
% Activate GUI and make the tree
params = {'cell.label', 'protocolSettings.lightAmplitude'};

exportPath16 = fullfile(dataFolder16, exportName16); 

%tree.visualize;
%%

function tree = loadTree(eportPath,dataFolder, params) 
    list = riekesuite.analysis.loadEpochList(exportPath, dataFolder);
    tree = riekesuite.analysis.buildTree(list, params);
end 

function out = getEpochBlock(tree, cellType)

child = 0;
for ii=1:node.children.length
    if node.children(ii).splitValue == 0.01
        epochBlockMatrix = riekesuite.getResponseMatrix(node.children(ii), 'Amp1'); **********************************************************************************************************************************************88********