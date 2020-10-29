function avgPsth = JennaPSTH(EpochData, preTime, binSize, adjustPlease)

if nargin == 4
    adjustPlease = adjustPlease;
else 
    adjustPlease = false;
end

% Find time/index where spike occurs. Stores per trace as a
% matrix. 
R = SpikeDetector(EpochData);
sp = R.sp;

epochLength = size(EpochData, 2); 
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