import jauimodel.*
import vuidocument.*
import celllabels.*

% exportName16 = "20200716RodsPSTH.mat";
% exportName21 = "20200721RodsPSTH.mat";
% create the full file path
%exportPath16 = fullfile(datafolder, exportName16);
%exportPath21 = fullfile(datafolder, exportName21);

exportNames = {...
    "20200716RodsPSTH.mat" "20200721RodsPSTH.mat"...
    "20200924RodsPSTH.mat" "20201007RodsPSTH.mat"...
    "20201022RodsPSTH.mat"};
datafolder = "/Users/jnagy2/Documents/errthing";
%%
% create the epoch list
%list = riekesuite.analysis.loadEpochList(exportPath16, datafolder);
%list2 = riekesuite.analysis.loadEpochList(exportPath21, datafolder);
% for ii=1:length(list2.elements)
%     list.append(list2.elements(ii))
% end
%% TODO Turn into function in Nova visio - move everythin there.
% Iterate through experiment dates, create list and add to lists. 
for ii=1:length(exportNames)
    exportName = exportNames{ii};
    exportPath = fullfile(datafolder, exportName);
    if ii == 1
        list = riekesuite.analysis.loadEpochList(exportPath, datafolder);
    else
        list_temp = riekesuite.analysis.loadEpochList(exportPath, datafolder);
        % Append elements from list for that experiment date
        for jj=1:length(list_temp.elements)
           list.append(list_temp.elements(jj)); 
        end
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
% Iterate: celltype, strain, epoch
for ii=1:numel(cellTypes)
    for mm =1: numel(cellStrains)
        for jj=1:tree.children.length
            epochs = tree.children(jj).epochList;
            cellKeys= {};
            ck = 0;

            % Get Pretime for spike rate calcu
            SampleEpoch = epochs.firstValue;
            preTime= SampleEpoch.protocolSettings.get('preTime') * 10;
            stimTime= SampleEpoch.protocolSettings.get('stimTime') * 10;
            tailTime= SampleEpoch.protocolSettings.get('tailTime') * 10;

            ll = 0;
            clear psth responses;
            traces = riekesuite.getResponseMatrix(epochs, 'Amp1');
            for kk=1:numel(epochs.elements)
               epoch = epochs.elements(kk);
               year = num2str(epoch.startDate(1));
               month = num2str(epoch.startDate(2).', '%02d');
               day = num2str(epoch.startDate(3).', '%02d');
               celllabel = epoch.cell.label.toCharArray';
               date_key = convertStringsToChars(strcat(year, month, day));
               
               cellType = cellLabels(date_key);
                   
               if ismember(celllabel, cellType.keys())
                   cellType = cellType(celllabel);
                   if strcmp(cellType{1}, cellTypes(ii)) && strcmp(cellType{2}, cellStrains{mm})
                       ll = ll + 1;
                       cellKey = strcat(date_key, celllabel);
                       if ~ismember(cellKey, cellKeys)
                         ck = ck + 1;
                         cellKeys{ck} = cellKey;
                       end
               
                       trace = traces(kk,:);
                       responses(ll,:) = trace;
                   end
               end 
            end
            if exist('responses', 'var')
                psthBlock = struct();
                if strcmp(cellTypes{ii}, "OFF S") || strcmp(cellTypes{ii}, "OFF T")
                    psthBlock.psth = calcPsth(responses, preTime, binSize, true);
                else
                    psthBlock.psth = calcPsth(responses, preTime, binSize, false);
                end
                psthBlock.cellType = cellTypes(ii);
                psthBlock.strain = cellStrains{mm};
                psthBlock.lightAmp = num2str(tree.children(jj).splitValue);
                psthBlock.numCells = num2str(length(unique(cellKeys)));
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
                    displayName = strcat(block.strain, ' ', '(', block.numCells, ')');
                    plot((1:length(block.psth)) * (binSize)/10000, block.psth, 'displayName', displayName, 'color', colors(convertStringsToChars(block.strain)))
                end
            end

        title(strcat(cellType, ": Light Intensity ",  amp, " with binsize ", num2str(binSize/1000), " seconds"));
        legend()
        xlabel("Time (seconds)");
        ylabel("Spike Rate (Hz)");
        saveas(fig, strcat(cellType, "_", amp, "_bz_", num2str(binSize), ".png"));
        hold off 
    end 
end 
%%
function psth = calcPsth(traces, preTime, binSize, isOffCell)
R = SpikeDetector(traces);       
sp = R.sp;

epochLength = size(traces,2);
nBins = floor(epochLength / binSize); 

psthCellAmp = zeros(length(sp), nBins);
for m=1:length(sp)
    epochSpikes = sp{m};

    if ~isempty(epochSpikes)

        %[epochPsth, edges] = histcounts(epochSpikes, nBins); 
            
        lowerBound = 0; 
        upperBound = binSize;
        baseSpikeCount = length(find(epochSpikes<=preTime));
        baseSpikeRate = (baseSpikeCount/preTime) * 10000;
        
        out = zeros(1, nBins);
        for i=1:nBins
            spikeCount = sum((epochSpikes > lowerBound) & (epochSpikes <= upperBound)) * 10000 / binSize;  
            if isOffCell == true 
                spikeCountAdj = spikeCount;
            else
                spikeCountAdj = spikeCount - baseSpikeRate;
            end
            out(i) = spikeCountAdj; 

            lowerBound = upperBound;
            upperBound = upperBound + binSize; 
        end 

        psthCellAmp(m,:) = out;
    end    
end 
psth = mean(psthCellAmp);
end