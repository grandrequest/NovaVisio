function avgPsth = PSTH(node, binSize, adjustPlease)

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