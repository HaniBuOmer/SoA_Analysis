%% Reject non-Brain ICs and Calculate PSD for Scalp EEG and IC time series

% initialization
workingDir      = extractBefore(matlab.desktop.editor.getActiveFilename,'p1000_code'); % ousdide "p1000_code" folder
preprocessPath  = append(workingDir, 'p2100_import_preprocess\');
PsdOutputPath      = append(workingDir, 'p2200_blockPsd\');

% freqBins        = 513;   % (EEG.srate*2+1)                     

% preprocessed data directory
allSets         = dir( append(preprocessPath,'*.set') );  %

for setIdx = 1:length(allSets)
    
    loadName = allSets(setIdx).name;
    EEG = pop_loadset('filename', loadName, 'filepath', preprocessPath);
    freqBins        = EEG.srate*2+1; 

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% Perform IC rejection. %%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Obtain the most dominant class label and its label probability.
    [~, mostDominantClassLabelVector] = max(EEG.etc.ic_classification.ICLabel.classifications, [], 2);
    mostDominantClassLabelProbVector = zeros(length(mostDominantClassLabelVector),1);
    for icIdx = 1:length(mostDominantClassLabelVector)
        mostDominantClassLabelProbVector(icIdx)  = EEG.etc.ic_classification.ICLabel.classifications(icIdx, mostDominantClassLabelVector(icIdx));
    end
    
    % Identify brain ICs. The order of the classes are {'Brain'  'Muscle'  'Eye'  'Heart'  'Line Noise'  'Channel Noise'  'Other'}.
    brainLabelProbThresh  = 0; % [0-1]
    brainIdx = find((mostDominantClassLabelVector==1 & mostDominantClassLabelProbVector>=brainLabelProbThresh));
    
    % Perform IC rejection.
    EEG = pop_subcomp(EEG, brainIdx, 0, 1);  % pop_subcomp( INEEG, components, plotag, keepcomp);
    
    % Post-process to update ICLabel data structure.
    EEG.etc.ic_classification.ICLabel.classifications = EEG.etc.ic_classification.ICLabel.classifications(brainIdx,:);
    
    % Post-process to update EEG.icaact.
    EEG.icaact = [];   % not necessary; pop_subcomp can do this
    EEG = eeg_checkset(EEG, 'ica');
    
    
    %%%%%%%%%%%%%%%%%%%%%%
    %%% Calculate PSD. %%%
    %%%%%%%%%%%%%%%%%%%%%%
    allEvents = {EEG.event.type}';
    uniqueEvents = unique(allEvents);
    blockOnsetEventIdx    = find(contains(uniqueEvents, '_start'));
    blockOnsetEventLabels = uniqueEvents(blockOnsetEventIdx);
    
    EEGorig = EEG;
    elecSpectraTensor = zeros(EEG.nbchan            , freqBins,length(blockOnsetEventLabels));
    icSpectraTensor   = zeros(size(EEG.icaweights,1), freqBins,length(blockOnsetEventLabels));
    
    for eventIdx = 1:length(blockOnsetEventLabels)
        
        currentBlock = blockOnsetEventLabels(eventIdx);
        EEG = pop_epoch(EEGorig, currentBlock, [0  20], 'newname', currentBlock{1,1}, 'epochinfo', 'yes');
        EEG = eeg_checkset(EEG, 'ica');
        
        [elecSpectra,freqs] = spectopo( EEG.data,  EEG.pnts, EEG.srate, 'freqfac', 4, ...
            'overlap', round(EEG.srate/2), 'freqrange', [1 50], 'plot', 'off');
        [icSpectra, ~ ]     = spectopo(EEG.icaact, EEG.pnts, EEG.srate, 'freqfac', 4, ...
            'overlap', round(EEG.srate/2), 'freqrange', [1 50], 'plot', 'off');
        elecSpectraTensor(:,:,eventIdx) = elecSpectra;
        icSpectraTensor(  :,:,eventIdx) = icSpectra;
    end
    
    EEG = EEGorig;
    EEG.etc.Hani.elecSpectraTensor = elecSpectraTensor;
    EEG.etc.Hani.icSpectraTensor   = icSpectraTensor;
    EEG.etc.Hani.freqs             = freqs;
    EEG.etc.Hani.labelOrder        = blockOnsetEventLabels;
    EEG                            = pop_saveset(EEG, 'filename', loadName, 'filepath', PsdOutputPath);
end

%%