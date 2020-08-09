function tree = createLoader(exportName, dataFolder)
% jauimodel stuff
loader = edu.washington.rieke.Analysis.getEntityLoader(); 
treeFactory = edu.washington.rieke.Analysis.getEpochTreeFactory();

import jauimodel.*
import vuidocument.*
exportPath = fullfile(dataFolder, exportName); 

list = riekesuite.analysis.loadEpochList(exportPath, dataFolder);

% Activate GUI and make the tree

keywordSplitter1= @(epoch)splitOnCellDate(epoch);
keywordSplitter_java1 = riekesuite.util.SplitValueFunctionAdapter.buildMap(list, keywordSplitter1);

tree = riekesuite.analysis.buildTree(list, {keywordSplitter_java1,'protocolSettings.epochGroup:label','protocolID','cell.label', 'protocolSettings.lightAmplitude'});

end 