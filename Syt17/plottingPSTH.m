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

% Used to epochList - could pass in filename instead
exportPath16 = fullfile(dataFolder16, exportName16); 
exportPath21 = fullfile(dataFolder16, exportName21); 

% Get list of epochs. Can use line below to add other experiments. 
list = riekesuite.analysis.loadEpochList(exportPath16, dataFolder16);

%%
%list21 = riekesuite.analysis.loadEpochList(exportPath21, dataFolder16);

%%

tree = riekesuite.analysis.buildTree(list, params);
tree.visualize;
%%
gui = epochTreeGUI(tree)
% gui = epochTreeGUI(tree);
%%
% Grab point in tree where we split on cell/branch above where you are
% doing analysis. 
node = tree; 
node.visualize
%%
% Function handle that passes in an epoch and grabs the split value
getAmpValues = @(epoch) epoch.splitValue;
amps = unique(arrayfun(getAmpValues, tree.leafNodes.elements));
%%
% this is for 0716
cellInfo = struct();
cellInfo.Cell11 = {'ON S' 'KO'};
cellInfo.Cell2 = {'ON S' 'WT'};
cellInfo.Cell4 = {'ON S' 'WT'};
cellInfo.Cell5 = {'ON S' 'WT'};
cellInfo.Cell7 = {'ON S' 'KO'};
cellInfo.Cell3 = {'OFF T' 'WT'};
cellInfo.Cell6 = {'OFF T' 'WT'};
cellInfo.Cell8  = {'OFF T' 'KO'};
cellInfo.Cell12 = {'OFF T' 'KO'};
cellInfo.Cell1 = {'OFF S' 'WT'};
cellInfo.Cell9  = {'OFF S' 'KO'};
cellInfo.Cell10  = {'OFF S' 'KO'};

cellTypes = {'ON S' 'OFF T' 'OFF S'};


%%
%this is for 0728
% cellInfo = struct();
% 
% cellInfo.Cell2 = {'OFF S' 'WT'};
% cellInfo.Cell4 = {'ON S' 'WT'};
% cellInfo.Cell5 = {'OFF S' 'WT'};
% cellInfo.Cell6 = {'ON S' 'WT'};
% cellInfo.Cell1 = {'ON S' 'WT'};
% 
% 
% cellTypes = {'ON S' 'OFF T' 'OFF S'};

%%

%this is for 0806

% cellInfo = struct();
% 
% cellInfo.Cell7 = {'ON S' 'KO'};
% cellInfo.Cell4 = {'ON S' 'KO'};
% cellInfo.Cell2 = {'ON S' 'KO'};
% cellInfo.Cell1 = {'ON S' 'KO'};
% 
% 
% cellTypes = {'ON S' 'OFF T' 'OFF S'};

%%
clear cellPsth chartData


% Iterate Cell, Amplitude
binSize = 30 * 10;
idx = 1; 
for ii=1:length(cellTypes)
    cellType = cellTypes{ii}; 
    
    for cc = 1:node.children.length
        clear Mean_pre Mean_respBS EpochData SampleEpoch PrePts R sp avgPsth
       
        cell = node.children(cc);
        thisCellInfo = cellInfo.(cell.splitValue);
        cellSplitValue = cellInfo.(cell.splitValue);
        
        if strcmp(thisCellInfo{1}, cellType)         
            for aa = 1:cell.children.length
                psthBlock = struct();
                if cellType == "ON S"
                    avgPsth = JennaPSTH(cell.children(aa), binSize);
                else
                    avgPsth = JennaPSTH(cell.children(aa), binSize, false);
                end
                
                psthBlock.psth = avgPsth;
                psthBlock.cellType = thisCellInfo{1};
                psthBlock.strain = thisCellInfo{2}; 
                psthBlock.lightAmp = num2str(cell.children(aa).splitValue);
                psthBlock.label = cell.splitValue;
                psthBlock.name = char(datestr((cell.children(aa).epochList.firstValue.cell.startDate)','yyyymmdd'));
                chartData{idx} = psthBlock;
                idx = idx + 1;
            end
        end 
    end
end 
%%
colors = containers.Map; 
colors('KO') = 'r';
colors('WT') = 'b';

for ii=1:length(cellTypes)
        
    for jj=1:length(amps)
        fig = figure;
        hold on; 
        grid on;

        cellType = cellTypes{ii}; 
        amp = num2str(amps(jj)); 

        for kk=1:length(chartData) 
            block = chartData{kk}; 

            if strcmp(block.cellType, cellType) && strcmp(amp, block.lightAmp)
                plot((1:length(block.psth)) * (binSize)/10000, block.psth, 'displayName', strcat(block.strain, " ", block.label), 'color', colors(block.strain))
            end
                    Name= strcat(block.name,strcat(block.strain, " ", block.label),'.txt');
        dlmwrite(Name, block.psth', 'delimiter', '\t', 'newline', 'unix');
        end 
        
        title(strcat(cellType, ": Light Intensity ",  amp));
        legend()
        xlabel("Time (seconds)");
        ylabel("Spike Rate (Hz)");
        saveas(fig, strcat(cellType, "_", amp, ".png"));
        hold off 
    end 

end
