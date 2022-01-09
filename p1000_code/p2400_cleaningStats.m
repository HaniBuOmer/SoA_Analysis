%% Cleaning statistics
% 05/18/2021 Hani.   Organized.
% 04/22/2021 Makoto. Created.


% initialization
workingDir      = extractBefore(matlab.desktop.editor.getActiveFilename,'p1000_code'); % ousdide "p1000_code" folder
preprocessPath  = append(workingDir, 'p2100_import_preprocess\');

allSets         = dir( append(preprocessPath,'*.set') );  

cleaningStats   = zeros(length(allSets), 3);

for setIdx = 1:length(allSets)
    
   loadName = allSets(setIdx).name;
   EEG = pop_loadset('filename', loadName, 'filepath', preprocessPath);
   [maxVal, maxIdx] = max(EEG.etc.ic_classification.ICLabel.classifications,[],2);
   brainIdx = find(maxIdx==1);
   brainProb = maxVal(brainIdx);
   beforeRejVar = var(EEG.data, 0, 2);
   EEG2 = pop_subcomp(EEG, brainIdx, 0, 1);
   afterRejVar  = var(EEG2.data, 0 ,2);
   varRej = 1-afterRejVar./beforeRejVar;
   cleaningStats(setIdx,:) = [size(EEG.icaweights, 1)-length(brainIdx) mean(varRej) mean(brainProb)];
    
end

%%% cleaningStats:  1) rejected ICs,
%%%                 2) reduced electrodes data variance,
%%%                 3) probability of Brain data
% icRej: 9.9 (2.5, 6-16); varRej: 0.59 (0.16, 0.32-0.89); brainProb: 0.82 (0.10, 0.55-0.98)
meanVal = mean(cleaningStats) 
stdVal  = std(cleaningStats)
minVal  = min(cleaningStats)
maxVal  = max(cleaningStats)
