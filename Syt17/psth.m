%import createLoader

import jauimodel.*
import vuidocument.*

%Data and export folder paths
dataFolder = "/Volumes/homedir/jnagy2/Data/Syt17 Analysis/20200716";
exportName = "20200716RodsPSTH.mat";

%% Create Parameters for splitting
% Activate GUI and make the tree
params = {'cell.label', 'protocolSettings.lightAmplitude'};

% Used to epochList - could pass in filename instead
exportPath = fullfile(dataFolder, exportName); 

%% Create Tree with paramters

% Get list of epochs. Can use line below to add other experiments. 
list = riekesuite.analysis.loadEpochList(exportPath, dataFolder);
% list.append(otherExportPath, otherDataForlder)

tree = riekesuite.analysis.buildTree(list, params);
tree.visualize
% gui = epochTreeGUI(tree);
%%
% Grab point in tree where we split on cell/branch above where you are
% doing analysis. 
node = tree; 
node.visualize
%%
getAmpValues = @(epoch) epoch.splitValue;
amps = unique(arrayfun(getAmpValues, tree.leafNodes.elements));
%%
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
clear cellPsth chartData

txtname= char(datestr((node.epochList.firstValue.cell.startDate)','yyyymmdd'));

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
                    avgPsth = PSTHsymphony3(cell.children(aa), binSize);
                else
                    avgPsth = PSTHsymphony3(cell.children(aa), binSize, false);
                end
                
                psthBlock.psth = avgPsth;
                psthBlock.cellType = thisCellInfo{1};
                psthBlock.strain = thisCellInfo{2}; 
                psthBlock.lightAmp = num2str(cell.children(aa).splitValue);
                psthBlock.label = cell.splitValue;
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
        end 
        
        title(strcat(cellType, ": Light Intensity ",  amp));
        legend()
        xlabel("Time (seconds)");
        ylabel("Spike Rate (Hz)");
        saveas(fig, strcat(cellType, "_", amp, ".png"));
        hold off 
    end 

end
%%
function avgPsth = PSTHsymphony3(node, binSize, adjustPlease)

if nargin == 3
    adjustPlease = adjustPlease;
else 
    adjustPlease = false;
end
   

% Get response matrix for each cell and light amplitude
epochs = node.epochList;
EpochData = riekesuite.getResponseMatrix(epochs, 'Amp1');

% Find time/index where spike occurs. Stores per trace as a
% matrix. 
R = SpikeDetector(EpochData);
sp = R.sp;

% Get Pretime for spike rate calcu
SampleEpoch = node.epochList.firstValue;
preTime= SampleEpoch.protocolSettings.get('preTime') * 10;
stimTime= SampleEpoch.protocolSettings.get('stimTime') * 10;
tailTime= SampleEpoch.protocolSettings.get('tailTime') * 10;

epochLength = preTime + stimTime + tailTime;
nBins = floor(epochLength / binSize); 

%nBins = length(SampleEpoch)

clear spCount_pre spikeRate_pre spCount_resp spikeRate_resp spikeRate_respBS hists psthCellAmp
for m=1:length(sp)
    epochSpikes = sp{m};

    if ~isempty(epochSpikes)

        %[epochPsth, edges] = histcounts(epochSpikes, nBins); 
            
        lowerBound = 0; 
        upperBound = binSize;
        baseSpikeCount = length(find(epochSpikes<=preTime));
        baseSpikeRate = (baseSpikeCount/preTime) * 10000;

        for i=1:nBins
            spikeCount = sum((epochSpikes > lowerBound) & (epochSpikes <= upperBound)) * 10000 / binSize;  
            if adjustPlease == true  
                spikeCountAdj = spikeCount - baseSpikeRate;
            else
                spikeCountAdj = spikeCount;
            end 
            out(i) = spikeCountAdj; 

            lowerBound = upperBound;
            upperBound = upperBound + binSize; 
        end 

        psthCellAmp(m,:) = out;
    end    
end 

avgPsth = mean(psthCellAmp);
end

