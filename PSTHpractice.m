%yay! 

% import createLoader
import jauimodel.*
import vuidocument.*
import celllabels.*

% exportName16 = "20200716RodsPSTH.mat";
% exportName21 = "20200721RodsPSTH.mat";

exportNames = {...
    "20200716RodsPSTH.mat" "20200721RodsPSTH.mat"...
    "20200723RodsPSTH.mat" "20200728RodsPSTH.mat"...
    "20200924RodsPSTH.mat" "20201007RodsPSTH.mat"...
    "20201022RodsPSTH.mat"};
datafolder = "/Users/jnagy2/Documents/errthing";

% create the full file path
exportPath16 = fullfile(datafolder, exportName16);
exportPath21 = fullfile(datafolder, exportName21);

%%
% create the epoch list
%list = riekesuite.analysis.loadEpochList(exportPath16, datafolder);
%list2 = riekesuite.analysis.loadEpochList(exportPath21, datafolder);
% for ii=1:length(list2.elements)
%     list.append(list2.elements(ii))
% end
%%
list = riekesuite.analysis.loadEpochList();
% Iterate through experiment dates, create list and add to lists. 
for ii=1:length(exportNames)
    exportName = exportNames{ii};
    exportPath = fullfile(datafolder, exportName);
    
    list_temp = riekesuite.analysis.loadEpochList(exportPath, datafolder);
    
    % Append elements from list for that experiment date
    for jj=1:length(list_temp.elements)
       list.append(list_temp.elements(jj)); 
    end
    
end
%%
% Decide parameters for tree before making tree
params = {'protocolSettings.lightAmplitude'};
tree = riekesuite.analysis.buildTree(list, params);
%%
tree.visualize
%%
idx = 1;
binSize = 200;
cellStrains = {'KO' 'WT'};
for ii=1:numel(cellTypes)
    for mm =1: numel(cellStrains)
        for jj=1:tree.children.length
            epochs = tree.children(jj).epochList;

            % Get Pretime for spike rate calcu
            SampleEpoch = tree.children(1).epochList.firstValue;
            preTime= SampleEpoch.protocolSettings.get('preTime') * 10;
            stimTime= SampleEpoch.protocolSettings.get('stimTime') * 10;
            tailTime= SampleEpoch.protocolSettings.get('tailTime') * 10;

            ll = 1;
            clear psth responses;
            traces = riekesuite.getResponseMatrix(epochs, 'Amp1');
            for kk=1:numel(epochs.elements)
               epoch = epochs.elements(kk);
               year = num2str(epoch.startDate(1));
               day = num2str(epoch.startDate(3));
               if epoch.startDate(2) < 10
                  month  = strcat("0", num2str(epoch.startDate(2))) ;
               else
                   month = num2str(epoch.startDate(2));
               end
               celllabel = epoch.cell.label.toCharArray';
               date_key = convertStringsToChars(strcat(year, month, day)); 
               cellType = cellLabels(date_key);

               if ismember(celllabel, cellType.keys())
                   cellType = cellType(celllabel);
                   if strcmp(cellType{1}, cellTypes(ii)) && strcmp(cellType{2}, cellStrains{mm})
                       trace = traces(kk,:);
                       responses(ll,:) = trace;
                       ll = ll + 1;
                       currCellType = cellType{1};
                       currCellStrain = cellType{2};
                   end
               end 
            end
            if exist('responses', 'var')
                psthBlock = struct();
                psthBlock.psth = calcPsth(responses, preTime, stimTime, tailTime, binSize);
                psthBlock.cellType = cellTypes(ii);
                psthBlock.strain = cellStrains{mm};
                psthBlock.lightAmp = num2str(tree.children(jj).splitValue);
                chartData{idx} = psthBlock;
                idx = idx + 1;
            end 
        end
    end
end

%% Plot psth
colors = containers.Map; 
colors('KO') = 'b';
colors('WT') = 'r';
getAmpValues = @(epoch) epoch.splitValue;
amps = unique(arrayfun(getAmpValues, tree.leafNodes.elements));

for ii=1:length(cellTypes)
    for jj=1:length(amps)
        fig = figure;
        hold on; 
        grid on;

            cellType = cellTypes{ii}; 
            amp = num2str(amps(jj)); 
            for kk=1:length(chartData) 
                block = chartData{kk}; 
                block.cellType, block.lightAmp
                if strcmp(block.cellType, cellType) && strcmp(amp, block.lightAmp)
                    plot((1:length(block.psth)) * (binSize)/10000, block.psth, 'displayName', block.strain, 'color', colors(convertStringsToChars(block.strain)))
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
function psth = calcPsth(traces, preTime, stimTime, tailTime, binSize)
R = SpikeDetector(traces);       
sp = R.sp;

epochLength = size(traces,2);
nBins = floor(epochLength / binSize); 

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
                spikeCountAdj = spikeCount - baseSpikeRate;
            out(i) = spikeCountAdj; 

            lowerBound = upperBound;
            upperBound = upperBound + binSize; 
        end 

        psthCellAmp(m,:) = out;
    end    
end 
psth = mean(psthCellAmp);
end