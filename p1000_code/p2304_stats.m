%% statistics

% D0, D100, D200, ..., D1000, M, rest.
% labelSortingOrder = [1 3:11 2 12 13];
% labelSortingOrder = [13 1 12 3:11 2];  % rest, D0, M, D100, D200, ..., D1000
    %EEG.etc.Hani.labelOrder(labelSortingOrder)

% Indices of  condition under study
% selCodIdx      = [1 3:11 2 12 13]; % All conditions
selCodIdx      = [1 3:11 2]; % D0, D100, ..., D1000
nCond          = length(selCodIdx);

% initialization
workingDir     = extractBefore(matlab.desktop.editor.getActiveFilename,'p1000_code'); % ousdide "p1000_code" folder
blockPsd_Path  = append(workingDir, 'p2200_blockPsd\');
OutputPath     = append(workingDir, 'p2300_stats\stats\');
spaceplotsPath = append(workingDir, 'plottingTools\spaceplots_v3');

corrTimeVect    = 0:100:1000;  % for corrcoef.

% figure saving options : % file format & resolution
figFormat  = '-djpeg95';  figResolution  = '-r300';

allSets        = dir( append(blockPsd_Path,'*.set') );


EEG = pop_loadset('filename', 'subj011.set', 'filepath', blockPsd_Path);
deltaIdx = find(EEG.etc.Hani.freqs>=  1 & EEG.etc.Hani.freqs<  4);
thetaIdx = find(EEG.etc.Hani.freqs>=  4 & EEG.etc.Hani.freqs<  8);
alphaIdx = find(EEG.etc.Hani.freqs>=  8 & EEG.etc.Hani.freqs< 13);
beta1Idx = find(EEG.etc.Hani.freqs>= 13 & EEG.etc.Hani.freqs< 20);
beta2Idx = find(EEG.etc.Hani.freqs>= 20 & EEG.etc.Hani.freqs< 35);
gammaIdx = find(EEG.etc.Hani.freqs>= 35 & EEG.etc.Hani.freqs< 50);


groupElecDelta = zeros(EEG.nbchan, nCond, length(allSets)); % elec, conditions, subjects.
groupElecTheta = zeros(EEG.nbchan, nCond, length(allSets)); % elec, conditions, subjects.
groupElecAlpha = zeros(EEG.nbchan, nCond, length(allSets)); % elec, conditions, subjects.
groupElecBeta1 = zeros(EEG.nbchan, nCond, length(allSets)); % elec, conditions, subjects.
groupElecBeta2 = zeros(EEG.nbchan, nCond, length(allSets)); % elec, conditions, subjects.
groupElecGamma = zeros(EEG.nbchan, nCond, length(allSets)); % elec, conditions, subjects.

for setIdx = 1:length(allSets)
    
    loadName = allSets(setIdx).name;
    EEG = pop_loadset('filename', loadName, 'filepath', blockPsd_Path);

%     if setIdx == 1
%         deltaIdx = find(EEG.etc.Hani.freqs>=1 & EEG.etc.Hani.freqs<4);
%         thetaIdx = find(EEG.etc.Hani.freqs>=4 & EEG.etc.Hani.freqs<8);
%         alphaIdx = find(EEG.etc.Hani.freqs>=8 & EEG.etc.Hani.freqs<13);
%         beta1Idx = find(EEG.etc.Hani.freqs>=13 & EEG.etc.Hani.freqs<20);
%         beta2Idx = find(EEG.etc.Hani.freqs>=20 & EEG.etc.Hani.freqs<35);
%         gammaIdx = find(EEG.etc.Hani.freqs>=35 & EEG.etc.Hani.freqs<50);
%     end

    % Perform one-way ANOVA for all channels/ICs.
    for elecIdx = 1:EEG.nbchan
        currentElecDelta = squeeze(mean(EEG.etc.Hani.elecSpectraTensor(elecIdx,deltaIdx,selCodIdx),2));
        currentElecTheta = squeeze(mean(EEG.etc.Hani.elecSpectraTensor(elecIdx,thetaIdx,selCodIdx),2));
        currentElecAlpha = squeeze(mean(EEG.etc.Hani.elecSpectraTensor(elecIdx,alphaIdx,selCodIdx),2));
        currentElecBeta1 = squeeze(mean(EEG.etc.Hani.elecSpectraTensor(elecIdx,beta1Idx,selCodIdx),2));
        currentElecBeta2 = squeeze(mean(EEG.etc.Hani.elecSpectraTensor(elecIdx,beta2Idx,selCodIdx),2));
        currentElecGamma = squeeze(mean(EEG.etc.Hani.elecSpectraTensor(elecIdx,gammaIdx,selCodIdx),2));
    
        groupElecDelta(elecIdx,:,setIdx) = currentElecDelta';
        groupElecTheta(elecIdx,:,setIdx) = currentElecTheta';
        groupElecAlpha(elecIdx,:,setIdx) = currentElecAlpha';
        groupElecBeta1(elecIdx,:,setIdx) = currentElecBeta1';
        groupElecBeta2(elecIdx,:,setIdx) = currentElecBeta2';
        groupElecGamma(elecIdx,:,setIdx) = currentElecGamma';
    end
end


%% Plot bar graphs with statistics.

bins = 0: size(currentElecDelta,1)-1;  % no. of conditions under comparison
useYlim = 0;       % use y-limit; 0 = no, 1 = yes
YLims = [-15 5];   % common y-limit - not used in this script
saveFig = 1;

pvalStack1 = [];
statsResults = [];
pvalStack2 = [];
corrcoeff1 = [];
pvalCounter = 0;

%%% *********************    ***********    *********************

freqBandLabel = {'Delta' 'Theta' 'Alpha' 'Beta1' 'Beta2' 'Gamma'};
for freqBandIdx = 1:6
    
    switch freqBandIdx
        case 1
            currentData = groupElecDelta;
        case 2
            currentData = groupElecTheta;
        case 3
            currentData = groupElecAlpha;
        case 4
            currentData = groupElecBeta1;
        case 5
            currentData = groupElecBeta2;
        case 6
            currentData = groupElecGamma;
    end
    
    figure
    for elecIdx = 1:EEG.nbchan
        subplot(4,5,elecIdx)
        elecPower = squeeze(currentData(elecIdx,:,:));
        meanPower = mean(elecPower,2);
        smePower  = std(elecPower,0,2)/sqrt(length(meanPower));
        bar(bins, meanPower, 'FaceColor', [0.66 0.76 1])
        hold on
        errorbar(bins, meanPower, smePower, 'linestyle', 'none', 'color', [0 0 0])
        
        if elecIdx == 1
            xlabel('Interval (x100ms)')
            ylabel('Power (10\timeslog10 uV^2/Hz dB)')
        end
        
        if useYlim == 1
            ylim([-8 5]) % ylim(YLims)
        end
        
        %%% perform one-way ANOVA ......
        [stats, df, pvals] = statcond(num2cell (elecPower,2)', 'paired', 'on');
        statsResults = [statsResults; stats];
        
        %%% perform Pearson correlation
        lmInput = mean(elecPower,2);
        [R,P]   = corrcoef([corrTimeVect' lmInput]);
        
        
        title(sprintf('%s p=%.5f, r=%.2f (p=%.3f)', EEG.chanlocs(elecIdx).labels, pvals, R(2,1), P(2,1)))
        pvalCounter = pvalCounter + 1;
        pvalStack1  = [pvalStack1;pvals];
        pvalStack2  = [pvalStack2;P(2,1)];
        corrcoeff1  = [corrcoeff1; R(2,1)];
        
    end
    
    set(gcf, 'name', freqBandLabel{freqBandIdx}, 'numbertitle', 'off', 'position', [64 1 1537 1123])
    suptitle(freqBandLabel{freqBandIdx})  % sgtitle(figTitles{freqBandIdx})
    
    if saveFig == 1
        addpath(spaceplotsPath)
        %     spaceplots   % to remove whitespace.
        rmpath(spaceplotsPath)
        print( append(OutputPath,freqBandLabel{freqBandIdx}), figFormat, figResolution)
        close
    end
    
end

% Calculate false discovery rate.
correctedPvalThreshold1 = fdr(pvalStack1, 0.01); %
tTreshold = [tinv(correctedPvalThreshold1/2, df(1)) tinv(1-correctedPvalThreshold1/2, df(1))];
correctedPvalThreshold2 = fdr(pvalStack2, 0.01); %

%  ***********    %%%%%%%%%%%    ***********    %%%%%%%%%%%    ***********
%% Plot statistically significant electrodes (red) with FDR-corrected p < 0.01.

freqBandLabel = {'Delta' 'Theta' 'Alpha' 'Beta1' 'Beta2' 'Gamma'};
figure
for freqBandIdx = 1:6
            
    subplot(2,3,freqBandIdx)
    currentstatsResults = statsResults(1+EEG.nbchan*(freqBandIdx-1):EEG.nbchan*freqBandIdx);
    topoplot(currentstatsResults, EEG.chanlocs, 'maplimits', tTreshold)

    title(freqBandLabel{freqBandIdx})
    
    if freqBandIdx == 6
        currentPosition = get(gca,'position');
        colorbar
        set(gca, 'position', currentPosition)
    end
end
set(gcf, 'position', [247   409   918   594])
set(findall(gcf, '-property', 'fontsize'), 'fontsize', 14)

if saveFig == 1
%     addpath(spaceplotsPath)
%     spaceplots   % to remove whitespace.
%     rmpath(spaceplotsPath)
%     spaceplots   % to remove whitespace etc.
    print( append(OutputPath,'FstatsResults'), figFormat, figResolution)
    close
end

%%%               ******************************

% % Significant F-stats  p-val mask FDR p<0.01.

freqBandLabel = {'Delta' 'Theta' 'Alpha' 'Beta1' 'Beta2' 'Gamma'};
figure
for freqBandIdx = 1:6

    subplot(2,3,freqBandIdx)
    currentPvals = pvalStack1(EEG.nbchan*(freqBandIdx-1)+1:EEG.nbchan*freqBandIdx);
    topoplot(currentPvals<correctedPvalThreshold1, EEG.chanlocs)
    title(freqBandLabel{freqBandIdx})
    
    if freqBandIdx == 6
        currentPosition = get(gca,'position');
        colorbar
        set(gca, 'position', currentPosition)
    end
end
set(gcf, 'position', [247   409   918   594])
set(findall(gcf, '-property', 'fontsize'), 'fontsize', 14)

if saveFig == 1
    %     addpath(spaceplotsPath)
    %     spaceplots   % to remove whitespace.
    %     rmpath(spaceplotsPath)
    print( append(OutputPath,'topoMaskF'), figFormat, figResolution)
    close
end

%  ***********    %%%%%%%%%%%    ***********    %%%%%%%%%%%    ***********

%% Plot Peason's statistically significant electrodes (red) with FDR-corrected p < 0.01.

% Peason's corr coeff.
figure
for freqBandsIdx = 1:6
    switch freqBandsIdx
        case 1
            freqBandLabel = 'Delta';
        case 2
            freqBandLabel = 'Theta';
        case 3
            freqBandLabel = 'Alpha';
        case 4
            freqBandLabel = 'Beta1';
        case 5
            freqBandLabel = 'Beta2';
        case 6
            freqBandLabel = 'Gamma';
    end
            
    subplot(2,3,freqBandsIdx)
    currentCoeffs = corrcoeff1(1+EEG.nbchan*(freqBandsIdx-1):EEG.nbchan*freqBandsIdx);
    topoplot(currentCoeffs, EEG.chanlocs, 'maplimits', [-1 1])
    title(freqBandLabel)
    
    if freqBandsIdx == 6
        currentPosition = get(gca,'position');
        colorbar
        set(gca, 'position', currentPosition)
    end
end
set(gcf, 'position', [247   409   918   594])
set(findall(gcf, '-property', 'fontsize'), 'fontsize', 14)

if saveFig == 1
%     addpath(spaceplotsPath)
%     spaceplots   % to remove whitespace.
%     rmpath(spaceplotsPath)
%     spaceplots   % to remove whitespace etc.
    print( append(OutputPath,'topoPearsonR'), figFormat, figResolution)
    close
end

%%%               ******************************

% Significant Pearson's R p-val mask FDR p<0.01.
figure
for freqBandsIdx = 1:6
    switch freqBandsIdx
        case 1
            freqBandLabel = 'Delta';
        case 2
            freqBandLabel = 'Theta';
        case 3
            freqBandLabel = 'Alpha';
        case 4
            freqBandLabel = 'Beta1';
        case 5
            freqBandLabel = 'Beta2';
        case 6
            freqBandLabel = 'Gamma';
    end
            
    subplot(2,3,freqBandsIdx)
    currentPvals = pvalStack2(EEG.nbchan*(freqBandsIdx-1)+1:EEG.nbchan*freqBandsIdx);
    topoplot(currentPvals<correctedPvalThreshold2, EEG.chanlocs) %topoplot(currentPvals<0.0046, EEG.chanlocs)
    
    title(freqBandLabel)
    
    if freqBandsIdx == 6
        currentPosition = get(gca,'position');
        colorbar
        set(gca, 'position', currentPosition)
    end
end
set(gcf, 'position', [247   409   918   594])
set(findall(gcf, '-property', 'fontsize'), 'fontsize', 14)

if saveFig == 1
%     addpath(spaceplotsPath)
%     spaceplots   % to remove whitespace.
%     rmpath(spaceplotsPath)
    print( append(OutputPath,'topoPearsonP'), figFormat, figResolution)
    close
end

%  ***********    %%%%%%%%%%%    ***********    %%%%%%%%%%%    ***********

%% Line plot with with stats.  (Corrected)

% Prepare colors
customColor = lines(6);

% timePlot = 000:100:200;
timePlot = corrTimeVect;

figure

for freqIdx = 1:6
    
    switch freqIdx
        case 1
            elecPower = squeeze(mean(groupElecDelta(:,:,:),1));
            titleString = 'Delta (1-4 Hz)';
        case 2
            elecPower = squeeze(mean(groupElecTheta(:,:,:),1));
            titleString = 'Theta (4-8 Hz)';
        case 3
            elecPower = squeeze(mean(groupElecAlpha(:,:,:),1));
            titleString = 'Alpha (8-13 Hz)';
        case 4
            elecPower = squeeze(mean(groupElecBeta1(:,:,:),1));
            titleString = 'Beta1 (13-20 Hz)';
        case 5
            elecPower = squeeze(mean(groupElecBeta2(:,:,:),1));
            titleString = 'Beta2 (20-35 Hz)';
        case 6
            elecPower = squeeze(mean(groupElecGamma(:,:,:),1));
            titleString = 'Gamma (35-50 Hz)';
    end
    
    % Calculate the grand-mean power.
    meanPower = mean(elecPower,2);
    
    % Calculate Delta power.
    meanPower = meanPower-meanPower(1);   %  why ?? : compared to 0 delay
%     meanPower = meanPower;    % Absolute
    
    stdPower  = std(elecPower);
%     plot(meanPower, 'LineStyle',LineTypes{freqIdx}, 'color', customColor(freqIdx, :), 'marker', '.', 'markersize', 20, 'LineStyleMode');
    nameval = lineprops(freqIdx);  % Line plot properties
%     plot(timePlot, meanPower, nameval{:}, 'color', customColor(freqIdx, :), 'markersize', 20, 'LineWidth', 2);
    plot(meanPower, nameval{:}, 'color', customColor(freqIdx, :), 'markersize', 10, 'LineWidth', 2);
    
    if freqIdx == 1
        hold on
%         xlabel('Delay Interval (ms)')
        xlabel('Experiment Delay Conditions')
        ylabel('Detla Power (10 \times log10 uV^2/Hz)')
    end
end

legend({'Delta' 'Theta' 'Alpha' 'Beta1' 'Beta2' 'Gamma'}, 'location', 'northwest')
% title('\Delta Global EEG power against conditions') %title('Global EEG power against interval length')
title('\Delta Global EEG power against time delay') %title('Global EEG power against interval length')

hold off
conditions = [ "D-0" "D-100" "D-200" "D-300" "D-400" "D-500" "D-600" "D-700" "D-800" "D-900" "D-1000" ];
set(gca,'XTickLabel',conditions, 'xtick', 1:nCond)
axis square
xlim([0 12])
ylim([-0.5 3.25])

set(gcf, 'numbertitle', 'off', 'position', [1 1 800 800])
set(findall(gcf, '-property', 'fontsize'), 'fontsize', 14)

if saveFig == 1
    addpath(spaceplotsPath)
%     spaceplots   % to remove whitespace.
    rmpath(spaceplotsPath)
    print( append(OutputPath,'grandMeanPowers') , figFormat, figResolution)
    close
end

%  ***********    %%%%%%%%%%%    ***********    %%%%%%%%%%%    ***********

function nameval = lineprops(idx)
%   LINEPROPS(idx) will return a 1x6 cell array of name-value pairs that
%   specify a marker, linestyle, and color combination identifies by idx
%   which is a positive integer.  These name-value pairs can be applied
%   to to the inputs of many plotting function using,
%       plot(___, nameval{:});
%
%   If you're only varying linestyle and color, set the LineStyleOrder
%   and ColorOrder properties of the axes, though be aware of some
%   compatibility issues prior to Matlab r2019b (see footnotes [1,2]).
%
%   LINEPROPS('count') returns the number of property combinations.
%   combinations
%
%   Examples
%   ------------
%   figure()
%   hold on
%   for i = 1:10
%       nameval = lineprops(i);
%       plot(0:.1:1,rand(1,11)+i, nameval{:}, 'LineWidth', 1.5)
%   end
%   legend('Location','BestOutside')
%
% lineprops() is available at
% https://www.mathworks.com/matlabcentral/answers/697655#answer_579560
% Author: Adam Danz
% Copyright (c) 2020  All rights reserved

%% Input validation
narginchk(1,1)
if (ischar(idx) || isa(idx,'string')) && all(strcmpi(idx,'count'))
    returnCount = true;
else
    validateattributes(idx,{'numeric'},{'scalar','numel',1})
    returnCount = false;
end

%% line properties
% List a bunch of markers; they will be selected in
% order and then the selection will start again if
% there are more lines than markers.
markers = {'o','+','*','x','s','d','v','^','>','p','v','<','h'};

% List a bunch of colors; like the markers, they
% will be selected circularly.
colors = {'r' 'g' 'b' 'c' 'm' 'k'}; % not use 'w' or 'y' due to visibility on white axes

% Same with line styles
linestyles = {'-','--',':','-.'};

if returnCount
    nameval = lcm(numel(markers),lcm(numel(colors),numel(linestyles)));
else
    first = @(v)v{1};
    prop = @(options, idx)first(circshift(options,-idx+1));
    nameval = {'Marker', prop(markers,idx), 'LineStyle', prop(linestyles,idx), 'Color', prop(colors,idx)};
end
end
