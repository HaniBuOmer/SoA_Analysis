%% Import and preprocess EEG data for all subjects

clear
% initialization
workingDir        = extractBefore(matlab.desktop.editor.getActiveFilename,'p1000_code'); % ousdide "p1000_code" folder
rawdataPath       = append(extractBefore(workingDir, 'p4xxx202111'), 'p0000_rawAll\');
PreprocessOutPath = append(workingDir, 'p2100_import_preprocess\');
if ~isfolder(PreprocessOutPath)
    mkdir(PreprocessOutPath)
end

%%% dipfit_settings parameters  <<< Check their existance
% hdmFile     = 'E:\\eeglab2019_1\\plugins\\dipfit3.4\\standard_BEM\\standard_vol.mat';
% mriFile     = 'E:\\eeglab2019_1\\plugins\\dipfit3.4\\standard_BEM\\standard_mri.mat';
% chanFile    = 'E:\\eeglab2019_1\\plugins\\dipfit3.4\\standard_BEM\\elec\\standard_1005.elc';
hdmFile     = append(workingDir, 'headModelFiles\standard_vol.mat');
mriFile     = append(workingDir, 'headModelFiles\standard_mri.mat');
chanFile    = append(workingDir, 'headModelFiles\standard_1005.elc');

% pop_chanedit parameters
lookupFile  = chanFile;

allVhdr     = dir( append(rawdataPath,'*.vhdr') ); % 'E:\PhD_Hani\EEG_Project\DrMakoto\p0000_rawAll\*.vhdr'

for vhdrIdx = 1:length(allVhdr)
    tic
    % Step 1: Import data.
    currentVhdr = allVhdr(vhdrIdx).name;
    EEG = pop_loadbv(rawdataPath, currentVhdr);
    
    
    % Step 2: Downsample to 256 Hz.
    EEG = pop_resample(EEG, 256);
    
    
    % Step 3: High-pass filter (cutoff 0.5, TBW 0.5)
    EEG = pop_firws(EEG, 'fcutoff', 0.5, 'ftype', 'highpass', 'wtype', 'hamming', 'forder', 1690, 'minphase', 0);
    
    
    % Step 4: Low-pass filter (cutoff 55, TBW 5)
    EEG = pop_firws(EEG, 'fcutoff', 50, 'ftype', 'lowpass', 'wtype', 'hamming', 'forder', 170, 'minphase', 0);
    
    
    % Set channels
    EEG = pop_chanedit(EEG, 'changefield',{18 'labels' 'AFp9'},...
        'lookup',lookupFile,...
        'changefield',{17 'X' '-31'}, ...
        'changefield',{17 'Y' '95.1'},...
        'changefield',{17 'Z' '-20'},...
        'convert',{'cart2all'},...
        'eval','chans = pop_chancenter( chans, [],[]);');
    
    
    % Step 5: Apply clean_rawdata.
    originalEEG = EEG;
    %   EEG = clean_rawdata(EEG, arg_flatline, arg_highpass, arg_channel, arg_noisy, arg_burst, arg_window)
    EEG = clean_rawdata(EEG, -1, -1, -1, -1, 20, -1); % , 'availableRAM_GB', 16);
    
    
    % Step 6: Interpolate all the removed channels.
    EEG = pop_interp(EEG, originalEEG.chanlocs, 'spherical');
    
    
    % Step 7: Re-reference the data to average.
    EEG.nbchan = EEG.nbchan+1;
    EEG.data(end+1,:) = zeros(1, EEG.pnts);
    EEG.chanlocs(1,EEG.nbchan).labels = 'initialReference';
    EEG = pop_reref(EEG, []);
    EEG = pop_select( EEG,'nochannel',{'initialReference'});
    
    
    %  Step 8: Run ICA using calculated data.
    EEG_forICA = pop_resample(EEG, 100);
    EEG_forICA = pop_runica(EEG_forICA, 'extended',1,'interupt','off');
    EEG.icaweights = EEG_forICA.icaweights;
    EEG.icasphere  = EEG_forICA.icasphere;
    EEG = eeg_checkset(EEG, 'ica');
    
    
    %  Step 9: Fit equivalent current dipoles.
    %%% first coregister coordinates    %% manual: [0 0 0 0 0 -1.5708 1 1 1]
    %     if vhdrIdx == 1
    [~,coordinateTransformParameters] = coregister(EEG.chanlocs, chanFile, 'warp', 'auto', 'manual', 'off'); %'E:\eeglab2019_1\plugins\dipfit3.4\standard_BEM/elec/standard_1005.elc', 'warp', 'auto', 'manual', 'off');
    %     end
    EEG = pop_dipfit_settings( EEG, 'hdmfile',hdmFile,...
        'coordformat','MNI',...
        'mrifile',  mriFile,...
        'chanfile', chanFile,...
        'coord_transform',coordinateTransformParameters,... %'coord_transform',[-0.0019958 -12.9477 -3.0656 1.3835e-08 5.5862e-07 -1.5708 1 1 1] ,...
        'chansel',1:EEG.nbchan);
    EEG = pop_multifit(EEG, 1:EEG.nbchan,'threshold', 100, 'dipplot','off','plotopt',{'normlen' 'on'});

    
    % Step 10: Search for and estimate symmetrically constrained bilateral dipoles
    EEG = fitTwoDipoles(EEG, 'LRR', 35);
    
    
    % Step 11: Run ICLabel (Pion-Tonachini et al., 2019)
    EEG = iclabel(EEG, 'default');
    
   
    % Save the dataset
    EEG.setname = sprintf('subj%03d', vhdrIdx);
    EEG.subject = sprintf('subj%03d', vhdrIdx);
    EEG = pop_saveset(EEG, 'filename', sprintf('subj%03d', vhdrIdx), 'filepath', PreprocessOutPath);
    elapsedtime(vhdrIdx) = toc/60
end

%%