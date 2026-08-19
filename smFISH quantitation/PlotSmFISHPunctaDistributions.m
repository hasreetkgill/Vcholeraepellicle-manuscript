%% PlotSmFISHPunctaDistributions.m
% Author: Hasreet Gill, Bassler Lab, HHMI/Princeton University
% Revised with Claude (Sonnet 5, Anthropic)

% Plotting stage for the smFISH puncta pipeline -- reads the per-
% condition Struct_*.mat files written by QuantifyPunctaCellsMicro-
% colonies.m (the "data" cell-level struct array specifically; the
% companion Struct_Micro_*.mat microcolony-level files are never read
% here) and produces, per channel and per puncta-fit metric:
%   - AllCells and Filtered (WasSplit==false, plus optional area/aspect-
%     ratio range filters) box/violin plots of the metric split by cell
%     class (Green vs. NonGreen, i.e. FromBingr), each with individual
%     replicates jittered on top and a matching CSV of the raw values.
%   - Pooled and per-replicate empirical CDF plots of the same split.
%   - FarRed-only: the ratio of cells above vs. at/below a fixed
%     absolute threshold, averaged across replicates.
%   - Scatter of the metric vs. distance-to-nearest-microcolony
%     (MinDistMicro), plus the same box/violin/threshold-ratio family
%     re-run within each of two distance-binning schemes (a simple
%     <=500/>500 split, and six 100px-wide bins).
%
% This script is shared by BOTH experiment types QuantifyPunctaCells-
% Microcolonies.m supports (see localGetCellLabels below, which branches
% on "Coculture" vs "SpikeIn" in the group name) -- it doesn't need to
% know which CELL_MASK_MODE produced a given Struct_*.mat file, only
% what that condition's group name implies about its cell-class labels.
%
% Uses swarmchart (built-in) for jittered replicate points, but the
% half-violin plots depend on an external/third-party Violin() function
% that must be on the MATLAB path -- it is not a MATLAB built-in.
%
% IMPORTANT: the channel and metric loops below (both in MAKE OUTPUT
% FOLDERS and in MAIN LOOP) are currently hardcoded to `for c = 2` and
% `for m = 2` -- i.e. this run only processes channel index 2 (FarRed)
% and metric index 4 (FixedSigma... actually metricShort{2} = 'GaussArea'
% here, since m indexes metricShort/metricPretty) with the full ranges
% (1:numel(channels), 1:4) left as trailing comments. Restore those
% ranges to run every channel/metric combination instead of just this
% one.

clear; clc; close all

%% ================= USER SETTINGS =================

dataDir = uigetdir(pwd, 'Choose folder containing Struct_*.mat files');
if isequal(dataDir,0)
    error('No folder selected.')
end
cd(dataDir)

plotRoot = fullfile(dataDir, 'Plots');
if ~exist(plotRoot,'dir'); mkdir(plotRoot); end

% "Filtered" population controls: WasSplit==false is always required;
% these two further narrow it to plausible single-cell size/shape.
useAreaFilter = true;
useAspectRatioFilter = true;

areaMin = 500;
areaMax = 1500;
aspectRatioMax = 3;
aspectRatioMin = 1.1;

% Optional: drop cells whose PixelList bounding box comes within
% boundaryDistPx of any image edge. Off by default.
excludeNearBoundary = false;
boundaryDistPx = 300;
imageWidthPx = 2000;
imageHeightPx = 2000;

% Per-metric (indexed like metricShort below) zero/near-zero exclusion
% toggles and thresholds, applied per channel before any plotting.
excludeZeroRed = [false false false false];
excludeZeroFarRed = [false false false false];
excludeZeroRedNorm = [false false false false];

zeroThrRed = [0 0 0 0];
zeroThrFarRed = [0 0 0 0];
zeroThrRedNorm = [0 0 0 0];

repColors = lines(20);
pointSize = 18;
pointAlpha = 0.55;
jitterAmount = 0.05;

% The 4 columns produced per puncta by extractpunctae in
% QuantifyPunctaCellsMicrocolonies.m: peak intensity, Gaussian-integrated
% area, background-subtracted raw sum, and the fixed-sigma peak
% intensity -- SumRed/SumFarRed store one 4-column row per cell (the sum
% over that cell's assigned puncta).
metricShort = {'Peak','GaussArea','BgSub','FixedSigma'};
metricPretty = {'Peak','Gaussian Area','Background-Subtracted Sum','Fixed-Sigma Peak'};

% Absolute FarRed threshold per metric, used only by the
% "ThresholdRatio" plots below (ratio of cells above vs. at/below this
% cutoff). Fill these in later with the values you want.
% [Peak, GaussArea, BgSub, FixedSigma]
ratioThreshold = [0.004 0.03 0.01 0.006];

% 2-bin MinDistMicro plots
% Note: values exactly equal to 500 are included in the <500 bin here
% to avoid dropping cells exactly at the boundary.
distBinNames2 = {'lt_500','gt_500'};
distBinPretty2 = {'MinDist < 500','MinDist > 500'};

% 100-pixel MinDistMicro plots
distBinNames100 = {'0_100','100_200','200_300','300_400','400_500','gt_500'};
distBinPretty100 = {'0 <= MinDist < 100','100 <= MinDist < 200','200 <= MinDist < 300', ...
                    '300 <= MinDist < 400','400 <= MinDist <= 500','MinDist > 500'};

%% ================= FILE DISCOVERY =================

matFiles = dir(fullfile(dataDir, 'Struct_*.mat'));
matFiles = matFiles(~startsWith(string({matFiles.name}), "Struct_Micro_"));
if isempty(matFiles)
    error('No Struct_*.mat files found in selected folder.')
end

%% ================= MAKE OUTPUT FOLDERS =================
% Pre-builds the full Plots/<Channel>/<Metric>/<AllCells|Filtered>/...
% tree so every plot type below can just imwrite/exportgraphics straight
% into it without checking/creating folders inline every time.

channels = {'Red','FarRed','RedNorm'};
analysisSets = {'AllCells','Filtered'};
baseSubSets = {'CDF_CombinedReplicates','CDF_ByReplicate','DistanceScatter'};
distBinSetNames = {'2Bins','100PixBins'};

for c = 2 %1:numel(channels)

    chanDir = fullfile(plotRoot, channels{c});
    if ~exist(chanDir,'dir'); mkdir(chanDir); end

    for m = 2 %1:4

        metricDir = fullfile(chanDir, metricShort{m});
        if ~exist(metricDir,'dir'); mkdir(metricDir); end

        for a = 1:numel(analysisSets)

            analysisDir = fullfile(metricDir, analysisSets{a});
            if ~exist(analysisDir,'dir'); mkdir(analysisDir); end

            for s = 1:numel(baseSubSets)
                thisDir = fullfile(analysisDir, baseSubSets{s});
                if ~exist(thisDir,'dir'); mkdir(thisDir); end
            end

            if strcmp(channels{c},'FarRed')
                thisDir = fullfile(analysisDir, 'ThresholdRatio');
                if ~exist(thisDir,'dir'); mkdir(thisDir); end
            end

            for bs = 1:numel(distBinSetNames)

                thisBinRoot = fullfile(analysisDir, 'DistBins', distBinSetNames{bs});
                if ~exist(thisBinRoot,'dir'); mkdir(thisBinRoot); end

                thisDir = fullfile(thisBinRoot, 'Box');
                if ~exist(thisDir,'dir'); mkdir(thisDir); end

                thisDir = fullfile(thisBinRoot, 'Violin');
                if ~exist(thisDir,'dir'); mkdir(thisDir); end

                if strcmp(channels{c},'FarRed')
                    thisDir = fullfile(thisBinRoot, 'ThresholdRatio');
                    if ~exist(thisDir,'dir'); mkdir(thisDir); end
                end
            end
        end
    end
end

%% ================= HELPERS =================
% Thin anonymous-handle aliases for the local functions defined at the
% bottom of this file (a style choice for shorter call-site names --
% MATLAB scripts can call their own local functions directly, so this
% indirection isn't strictly required, just kept for readability).

getProbeLabels = @(groupName) localGetProbeLabels(groupName);
getCellLabels = @(groupName) localGetCellLabels(groupName);
metricVecFromStruct = @(T,fieldName,idx) localMetricVecFromStruct(T,fieldName,idx);
applyMetricZeroFilter = @(T,metricVals,doExclude,thr) localApplyMetricZeroFilter(T,metricVals,doExclude,thr);
addReplicateSwarm = @(xBase,yVals,repVals,repColors,pointSize,pointAlpha,jitterAmount) ...
    localAddReplicateSwarm(xBase,yVals,repVals,repColors,pointSize,pointAlpha,jitterAmount);
makeSimpleBox = @(xPos,yVals,colorIn) localMakeSimpleBox(xPos,yVals,colorIn);
safeLogLimits = @(vals) localSafeLogLimits(vals);
saveFig = @(figHandle,filePath) localSaveFig(figHandle,filePath);
savePlotCSV = @(T,filePath) localSavePlotCSV(T,filePath);
makePlotData = @(T,groupName,thisChannel,metricName,analysisName,cellClass,plotGroup,xPos,distanceBin) ...
    localMakePlotData(T,groupName,thisChannel,metricName,analysisName,cellClass,plotGroup,xPos,distanceBin);

%% ================= MAIN LOOP =================

for c = 2 %1:numel(channels)

    thisChannel = channels{c};

    if strcmp(thisChannel,'Red')
        sumField = 'SumRed';
        zeroExcludeVec = excludeZeroRed;
        zeroThrVec = zeroThrRed;
    elseif strcmp(thisChannel,'FarRed')
        sumField = 'SumFarRed';
        zeroExcludeVec = excludeZeroFarRed;
        zeroThrVec = zeroThrFarRed;
    else
        % RedNorm is a derived channel (FarRed/Red ratio per metric),
        % not a real fluorescence column -- see metricVals below.
        sumField = 'RedNorm';
        zeroExcludeVec = excludeZeroRedNorm;
        zeroThrVec = zeroThrRedNorm;
    end

    for m = 2 %1:4

        thisRatioThreshold = ratioThreshold(m);

        fprintf('\n=== %s / %s ===\n', thisChannel, metricPretty{m});
        chanMetricDir = fullfile(plotRoot, thisChannel, metricShort{m});

        for f = 1:numel(matFiles)

            fpath = fullfile(matFiles(f).folder, matFiles(f).name);
            S = load(fpath);

            if ~isfield(S,'data')
                warning('Skipping %s: no variable named "data".', matFiles(f).name);
                continue
            end
            if isempty(S.data)
                warning('Skipping %s: data is empty.', matFiles(f).name);
                continue
            end

            [~, baseName, ~] = fileparts(matFiles(f).name);
            groupName = erase(baseName, 'Struct_');

            % RedNorm (vpsL/gyrA ratio) is only meaningful for the
            % probe pairing it was computed for.
            if strcmp(thisChannel,'RedNorm') && ~contains(groupName,'vpsLgyrA')
                continue
            end

            probeInfo = getProbeLabels(groupName);
            cellInfo = getCellLabels(groupName);

            T = struct2table(S.data);

            % These are exactly the fields QuantifyPunctaCellsMicro-
            % colonies.m's cellprops rows carry -- missing any of them
            % means this file didn't come from that pipeline (or from
            % an incompatible older version of it).
            requiredVars = {'Replicate','XYPos','AspectRatio','Area','MinDistMicro','FromBingr','WasSplit','PixelList'};
            missingVars = requiredVars(~ismember(requiredVars, T.Properties.VariableNames));
            if ~isempty(missingVars)
                warning('Skipping %s: missing required variables: %s', ...
                    matFiles(f).name, strjoin(missingVars, ', '));
                continue
            end

            if strcmp(thisChannel,'RedNorm')
                if ~ismember('SumRed', T.Properties.VariableNames) || ~ismember('SumFarRed', T.Properties.VariableNames)
                    warning('Skipping %s: missing SumRed or SumFarRed for RedNorm.', matFiles(f).name);
                    continue
                end
                redVals = metricVecFromStruct(T, 'SumRed', m);
                farRedVals = metricVecFromStruct(T, 'SumFarRed', m);
                metricVals = farRedVals ./ redVals;
            else
                if ~ismember(sumField, T.Properties.VariableNames)
                    warning('Skipping %s: missing %s', matFiles(f).name, sumField);
                    continue
                end
                metricVals = metricVecFromStruct(T, sumField, m);
            end

            T.MetricValue = metricVals;
            T.IsGreen = logical(T.FromBingr);
            T.IsNonGreen = ~T.IsGreen;

            % Treat "no microcolony present" cases (stored as 0) as very far from microcolonies.
            % These will fall into the >500 distance bin.
            T.MinDistMicro(T.MinDistMicro == 0) = 1000;

            T = applyMetricZeroFilter(T, T.MetricValue, zeroExcludeVec(m), zeroThrVec(m));
            if isempty(T)
                warning('Skipping %s after zero filtering: no rows remain.', matFiles(f).name);
                continue
            end

            % Every plot below uses a log y-axis for MetricValue, so any
            % non-positive or non-finite value has to be dropped first;
            % this reports exactly what's being dropped (split by cell
            % class) before actually filtering, as an audit trail rather
            % than a silent removal.
            valsForLog = T.MetricValue;

            negMask = isfinite(valsForLog) & valsForLog < 0;
            zeroMask = isfinite(valsForLog) & valsForLog == 0;
            nanMask = isnan(valsForLog);
            posInfMask = valsForLog == Inf;
            negInfMask = valsForLog == -Inf;

            validLog = isfinite(valsForLog) & valsForLog > 0;
            nRemoved = sum(~validLog);

            if nRemoved > 0

                fprintf('Removing %d values for log plotting from %s | %s / %s:\n', ...
                    nRemoved, matFiles(f).name, thisChannel, metricShort{m});

                fprintf('    Total: %d negative, %d zero, %d NaN, %d +Inf, %d -Inf\n', ...
                    sum(negMask), sum(zeroMask), sum(nanMask), sum(posInfMask), sum(negInfMask));

                classMasks = {T.IsNonGreen, T.IsGreen};
                classLabels = {cellInfo.NonGreenLabel, cellInfo.GreenLabel};

                for ccPrint = 1:numel(classMasks)

                    thisClassMask = classMasks{ccPrint};
                    thisClassLabel = classLabels{ccPrint};

                    fprintf('    %s: %d negative, %d zero, %d NaN, %d +Inf, %d -Inf\n', ...
                        thisClassLabel, ...
                        sum(negMask & thisClassMask), ...
                        sum(zeroMask & thisClassMask), ...
                        sum(nanMask & thisClassMask), ...
                        sum(posInfMask & thisClassMask), ...
                        sum(negInfMask & thisClassMask));
                end

            end

            T = T(validLog,:);

            if isempty(T)
                warning('Skipping %s after removing invalid log-scale values: no rows remain.', matFiles(f).name);
                continue
            end

            if excludeNearBoundary
                keepBoundary = true(height(T),1);
                for bb = 1:height(T)
                    pixList = T.PixelList{bb};
                    xpix = pixList(:,1);
                    ypix = pixList(:,2);
                    keepBoundary(bb) = min(xpix) > boundaryDistPx && ...
                                       min(ypix) > boundaryDistPx && ...
                                       max(xpix) < imageWidthPx - boundaryDistPx && ...
                                       max(ypix) < imageHeightPx - boundaryDistPx;
                end
                T = T(keepBoundary,:);
            end

            if isempty(T)
                warning('Skipping %s after boundary filtering: no rows remain.', matFiles(f).name);
                continue
            end

            % "Filtered" population: drop any cell that had to be
            % watershed-split from a neighbor, then optionally narrow to
            % a plausible single-cell area/aspect-ratio range.
            Tfilt = T;
            Tfilt = Tfilt(Tfilt.WasSplit == false, :);

            if useAreaFilter
                Tfilt = Tfilt(Tfilt.Area >= areaMin & Tfilt.Area <= areaMax, :);
            end
            if useAspectRatioFilter
                Tfilt = Tfilt(Tfilt.AspectRatio >= aspectRatioMin & Tfilt.AspectRatio <= aspectRatioMax, :);
            end

            if strcmp(thisChannel,'RedNorm')
                yLab = sprintf('vpsL/gyrA: %s', metricPretty{m});
            else
                yLab = sprintf('%s: %s', probeInfo.(thisChannel), metricPretty{m});
            end
            condTitle = strrep(groupName, '_', '\_');

            %% ---------- ALL CELLS BASIC PLOTS ----------
            % Box (hand-drawn, see makeSimpleBox) + jittered per-
            % replicate swarm points, a half-violin equivalent, a pooled
            % CDF, and a per-replicate CDF -- all on the unfiltered (but
            % log-valid) population T, split Green vs. NonGreen.
            green = T(T.IsGreen, :);
            nongreen = T(T.IsNonGreen, :);

            if ~isempty(green) || ~isempty(nongreen)

                fig = figure('Visible','off');
                hold on
                plotData = table();

                if ~isempty(nongreen)
                    addReplicateSwarm(1, nongreen.MetricValue, nongreen.Replicate, repColors, pointSize, pointAlpha, jitterAmount);
                    makeSimpleBox(1, nongreen.MetricValue, [0 0 0]);
                    plotData = [plotData; makePlotData(nongreen,groupName,thisChannel,metricShort{m},'AllCells',cellInfo.NonGreenLabel,'NonGreen',1,'')]; %#ok<AGROW>
                end
                if ~isempty(green)
                    addReplicateSwarm(2, green.MetricValue, green.Replicate, repColors, pointSize, pointAlpha, jitterAmount);
                    makeSimpleBox(2, green.MetricValue, [0 0.6 0]);
                    plotData = [plotData; makePlotData(green,groupName,thisChannel,metricShort{m},'AllCells',cellInfo.GreenLabel,'Green',2,'')]; %#ok<AGROW>
                end

                set(gca,'YScale','log')
                ylim(safeLogLimits([nongreen.MetricValue; green.MetricValue]))
                xlim([0.5 2.5])
                set(gca,'XTick',[1 2],'XTickLabel',{cellInfo.NonGreenLabel, cellInfo.GreenLabel})
                ylabel(yLab)
                title(sprintf('%s | all cells', condTitle), 'Interpreter','none')
                box on; grid on
                hold off

                outPath = fullfile(chanMetricDir,'AllCells',sprintf('%s_Box_%s_%s.png', groupName, thisChannel, metricShort{m}));
                savePlotCSV(plotData, localCsvPath(outPath));
                saveFig(fig, outPath)

                fig = figure('Visible','off');
                hold on

                if ~isempty(nongreen)
                    Violin({nongreen.MetricValue}, 1, 'HalfViolin', 'full', 'QuartileStyle', 'boxplot', 'DataStyle', 'scatter', 'ShowMean', false, 'ShowMedian', true, 'ShowNotches', false, 'ViolinColor', {[0 0 0]});
                end
                if ~isempty(green)
                    Violin({green.MetricValue}, 2, 'HalfViolin', 'full', 'QuartileStyle', 'boxplot', 'DataStyle', 'scatter', 'ShowMean', false, 'ShowMedian', true, 'ShowNotches', false, 'ViolinColor', {[0 0.6 0]});
                end

                set(gca,'YScale','log')
                ylim(safeLogLimits([nongreen.MetricValue; green.MetricValue]))
                xlim([0.5 2.5])
                set(gca,'XTick',[1 2],'XTickLabel',{cellInfo.NonGreenLabel, cellInfo.GreenLabel})
                ylabel(yLab)
                title(sprintf('%s | all cells', condTitle), 'Interpreter','none')
                box on; grid on
                hold off

                outPath = fullfile(chanMetricDir,'AllCells',sprintf('%s_Violin_%s_%s.png', groupName, thisChannel, metricShort{m}));
                saveFig(fig, outPath)
            end

            if ~isempty(green) || ~isempty(nongreen)

                fig = figure('Visible','off');
                plotData = table();

                subplot(1,2,1)
                if ~isempty(nongreen)
                    cdfplot(nongreen.MetricValue); hold on
                    plotData = [plotData; makePlotData(nongreen,groupName,thisChannel,metricShort{m},'AllCells',cellInfo.NonGreenLabel,'NonGreen',1,'')]; %#ok<AGROW>
                end
                set(gca,'XScale','log')
                xlabel(yLab); ylabel('CDF')
                title(sprintf('%s | %s', condTitle, cellInfo.NonGreenLabel), 'Interpreter','none')
                box on; grid on

                subplot(1,2,2)
                if ~isempty(green)
                    cdfplot(green.MetricValue); hold on
                    plotData = [plotData; makePlotData(green,groupName,thisChannel,metricShort{m},'AllCells',cellInfo.GreenLabel,'Green',2,'')]; %#ok<AGROW>
                end
                set(gca,'XScale','log')
                xlabel(yLab); ylabel('CDF')
                title(sprintf('%s | %s', condTitle, cellInfo.GreenLabel), 'Interpreter','none')
                box on; grid on

                outPath = fullfile(chanMetricDir,'AllCells','CDF_CombinedReplicates',sprintf('%s_CDFCombined_%s_%s.png', groupName, thisChannel, metricShort{m}));
                savePlotCSV(plotData, localCsvPath(outPath));
                saveFig(fig, outPath)
            end

            if ~isempty(green) || ~isempty(nongreen)

                fig = figure('Visible','off');
                plotData = table();

                subplot(1,2,1); hold on
                repsNG = unique(nongreen.Replicate);
                for rr = 1:numel(repsNG)
                    thisRep = repsNG(rr);
                    y = nongreen.MetricValue(nongreen.Replicate == thisRep);
                    if ~isempty(y)
                        [fCDF,xCDF] = ecdf(y);
                        plot(xCDF, fCDF, 'LineWidth', 1.5, 'Color', repColors(mod(thisRep-1,size(repColors,1))+1,:));
                        tmp = localMakeCDFData(xCDF,fCDF,groupName,thisChannel,metricShort{m},'AllCells',cellInfo.NonGreenLabel,sprintf('NonGreen_Rep%d',thisRep),'');
                        plotData = [plotData; tmp]; %#ok<AGROW>
                    end
                end
                set(gca,'XScale','log')
                xlabel(yLab); ylabel('CDF')
                title(sprintf('%s | %s by replicate', condTitle, cellInfo.NonGreenLabel), 'Interpreter','none')
                box on; grid on

                subplot(1,2,2); hold on
                repsG = unique(green.Replicate);
                for rr = 1:numel(repsG)
                    thisRep = repsG(rr);
                    y = green.MetricValue(green.Replicate == thisRep);
                    if ~isempty(y)
                        [fCDF,xCDF] = ecdf(y);
                        plot(xCDF, fCDF, 'LineWidth', 1.5, 'Color', repColors(mod(thisRep-1,size(repColors,1))+1,:));
                        tmp = localMakeCDFData(xCDF,fCDF,groupName,thisChannel,metricShort{m},'AllCells',cellInfo.GreenLabel,sprintf('Green_Rep%d',thisRep),'');
                        plotData = [plotData; tmp]; %#ok<AGROW>
                    end
                end
                set(gca,'XScale','log')
                xlabel(yLab); ylabel('CDF')
                title(sprintf('%s | %s by replicate', condTitle, cellInfo.GreenLabel), 'Interpreter','none')
                box on; grid on

                outPath = fullfile(chanMetricDir,'AllCells','CDF_ByReplicate',sprintf('%s_CDFByRep_%s_%s.png', groupName, thisChannel, metricShort{m}));
                savePlotCSV(plotData, localCsvPath(outPath));
                saveFig(fig, outPath)
            end

            %% ---------- FILTERED BASIC PLOTS ----------
            % Same four plot types as above, on Tfilt instead of T.
            greenF = Tfilt(Tfilt.IsGreen, :);
            nongreenF = Tfilt(Tfilt.IsNonGreen, :);

            if ~isempty(greenF) || ~isempty(nongreenF)

                fig = figure('Visible','off');
                hold on
                plotData = table();

                if ~isempty(nongreenF)
                    addReplicateSwarm(1, nongreenF.MetricValue, nongreenF.Replicate, repColors, pointSize, pointAlpha, jitterAmount);
                    makeSimpleBox(1, nongreenF.MetricValue, [0 0 0]);
                    plotData = [plotData; makePlotData(nongreenF,groupName,thisChannel,metricShort{m},'Filtered',cellInfo.NonGreenLabel,'NonGreen',1,'')]; %#ok<AGROW>
                end
                if ~isempty(greenF)
                    addReplicateSwarm(2, greenF.MetricValue, greenF.Replicate, repColors, pointSize, pointAlpha, jitterAmount);
                    makeSimpleBox(2, greenF.MetricValue, [0 0.6 0]);
                    plotData = [plotData; makePlotData(greenF,groupName,thisChannel,metricShort{m},'Filtered',cellInfo.GreenLabel,'Green',2,'')]; %#ok<AGROW>
                end

                set(gca,'YScale','log')
                ylim(safeLogLimits([nongreenF.MetricValue; greenF.MetricValue]))
                xlim([0.5 2.5])
                set(gca,'XTick',[1 2],'XTickLabel',{cellInfo.NonGreenLabel, cellInfo.GreenLabel})
                ylabel(yLab)
                title(sprintf('%s | filtered', condTitle), 'Interpreter','none')
                box on; grid on
                hold off

                outPath = fullfile(chanMetricDir,'Filtered',sprintf('%s_BoxFiltered_%s_%s.png', groupName, thisChannel, metricShort{m}));
                savePlotCSV(plotData, localCsvPath(outPath));
                saveFig(fig, outPath)

                fig = figure('Visible','off');
                hold on

                if ~isempty(nongreenF)
                    Violin({nongreenF.MetricValue}, 1, 'HalfViolin', 'full', 'QuartileStyle', 'boxplot', 'DataStyle', 'scatter', 'ShowMean', false, 'ShowMedian', true, 'ShowNotches', false, 'ViolinColor', {[0 0 0]});
                end
                if ~isempty(greenF)
                    Violin({greenF.MetricValue}, 2, 'HalfViolin', 'full', 'QuartileStyle', 'boxplot', 'DataStyle', 'scatter', 'ShowMean', false, 'ShowMedian', true, 'ShowNotches', false, 'ViolinColor', {[0 0.6 0]});
                end

                set(gca,'YScale','log')
                ylim(safeLogLimits([nongreenF.MetricValue; greenF.MetricValue]))
                xlim([0.5 2.5])
                set(gca,'XTick',[1 2],'XTickLabel',{cellInfo.NonGreenLabel, cellInfo.GreenLabel})
                ylabel(yLab)
                title(sprintf('%s | filtered', condTitle), 'Interpreter','none')
                box on; grid on
                hold off

                outPath = fullfile(chanMetricDir,'Filtered',sprintf('%s_ViolinFiltered_%s_%s.png', groupName, thisChannel, metricShort{m}));
                saveFig(fig, outPath)
            end

            if ~isempty(greenF) || ~isempty(nongreenF)

                fig = figure('Visible','off');
                plotData = table();

                subplot(1,2,1)
                if ~isempty(nongreenF)
                    cdfplot(nongreenF.MetricValue); hold on
                    plotData = [plotData; makePlotData(nongreenF,groupName,thisChannel,metricShort{m},'Filtered',cellInfo.NonGreenLabel,'NonGreen',1,'')]; %#ok<AGROW>
                end
                set(gca,'XScale','log')
                xlabel(yLab); ylabel('CDF')
                title(sprintf('%s | filtered %s', condTitle, cellInfo.NonGreenLabel), 'Interpreter','none')
                box on; grid on

                subplot(1,2,2)
                if ~isempty(greenF)
                    cdfplot(greenF.MetricValue); hold on
                    plotData = [plotData; makePlotData(greenF,groupName,thisChannel,metricShort{m},'Filtered',cellInfo.GreenLabel,'Green',2,'')]; %#ok<AGROW>
                end
                set(gca,'XScale','log')
                xlabel(yLab); ylabel('CDF')
                title(sprintf('%s | filtered %s', condTitle, cellInfo.GreenLabel), 'Interpreter','none')
                box on; grid on

                outPath = fullfile(chanMetricDir,'Filtered','CDF_CombinedReplicates',sprintf('%s_CDFFilteredCombined_%s_%s.png', groupName, thisChannel, metricShort{m}));
                savePlotCSV(plotData, localCsvPath(outPath));
                saveFig(fig, outPath)
            end

            if ~isempty(greenF) || ~isempty(nongreenF)

                fig = figure('Visible','off');
                plotData = table();

                subplot(1,2,1); hold on
                repsNG = unique(nongreenF.Replicate);
                for rr = 1:numel(repsNG)
                    thisRep = repsNG(rr);
                    y = nongreenF.MetricValue(nongreenF.Replicate == thisRep);
                    if ~isempty(y)
                        [fCDF,xCDF] = ecdf(y);
                        plot(xCDF, fCDF, 'LineWidth', 1.5, 'Color', repColors(mod(thisRep-1,size(repColors,1))+1,:));
                        tmp = localMakeCDFData(xCDF,fCDF,groupName,thisChannel,metricShort{m},'Filtered',cellInfo.NonGreenLabel,sprintf('NonGreen_Rep%d',thisRep),'');
                        plotData = [plotData; tmp]; %#ok<AGROW>
                    end
                end
                set(gca,'XScale','log')
                xlabel(yLab); ylabel('CDF')
                title(sprintf('%s | filtered %s by replicate', condTitle, cellInfo.NonGreenLabel), 'Interpreter','none')
                box on; grid on

                subplot(1,2,2); hold on
                repsG = unique(greenF.Replicate);
                for rr = 1:numel(repsG)
                    thisRep = repsG(rr);
                    y = greenF.MetricValue(greenF.Replicate == thisRep);
                    if ~isempty(y)
                        [fCDF,xCDF] = ecdf(y);
                        plot(xCDF, fCDF, 'LineWidth', 1.5, 'Color', repColors(mod(thisRep-1,size(repColors,1))+1,:));
                        tmp = localMakeCDFData(xCDF,fCDF,groupName,thisChannel,metricShort{m},'Filtered',cellInfo.GreenLabel,sprintf('Green_Rep%d',thisRep),'');
                        plotData = [plotData; tmp]; %#ok<AGROW>
                    end
                end
                set(gca,'XScale','log')
                xlabel(yLab); ylabel('CDF')
                title(sprintf('%s | filtered %s by replicate', condTitle, cellInfo.GreenLabel), 'Interpreter','none')
                box on; grid on

                outPath = fullfile(chanMetricDir,'Filtered','CDF_ByReplicate',sprintf('%s_CDFFilteredByRep_%s_%s.png', groupName, thisChannel, metricShort{m}));
                savePlotCSV(plotData, localCsvPath(outPath));
                saveFig(fig, outPath)
            end

            %% ---------- THRESHOLD RATIO PLOTS: FARRED ONLY ----------
            % ratio = (# MetricValue > thisRatioThreshold) / (# MetricValue <= thisRatioThreshold)
            % Computed per replicate, then averaged; individual replicate
            % ratios are overlaid as jittered points on the mean bar.

            if strcmp(thisChannel,'FarRed')

                ratioTables = {T, Tfilt};
                ratioNames = {'AllCells','Filtered'};

                for rset = 1:2

                    Tratio = ratioTables{rset};
                    ratioFolder = ratioNames{rset};

                    if isempty(Tratio)
                        continue
                    end

                    ratioData = table();

                    cellGroups = {
                        Tratio(Tratio.IsNonGreen,:), cellInfo.NonGreenLabel, 1, [0 0 0], 'NonGreen';
                        Tratio(Tratio.IsGreen,:),    cellInfo.GreenLabel,    2, [0 0.6 0], 'Green'
                    };

                    fig = figure('Visible','off');
                    hold on

                    for cg = 1:size(cellGroups,1)

                        Tcg = cellGroups{cg,1};
                        classLabel = cellGroups{cg,2};
                        xpos = cellGroups{cg,3};
                        thisColor = cellGroups{cg,4};
                        className = cellGroups{cg,5};

                        if isempty(Tcg)
                            continue
                        end

                        reps = unique(Tcg.Replicate);
                        repRatios = nan(numel(reps),1);

                        for rr = 1:numel(reps)

                            thisRep = reps(rr);
                            Trep = Tcg(Tcg.Replicate == thisRep,:);
                            vals = Trep.MetricValue;

                            nAbove = sum(vals > thisRatioThreshold);
                            nBelow = sum(vals <= thisRatioThreshold);

                            if nBelow == 0
                                ratioVal = NaN;
                            else
                                ratioVal = nAbove / nBelow;
                            end

                            repRatios(rr) = ratioVal;

                            tmp = localMakeRatioTable(groupName, thisChannel, metricShort{m}, ...
                                ratioFolder, classLabel, string(className), thisRep, ...
                                thisRatioThreshold, nAbove, nBelow, ratioVal);

                            ratioData = [ratioData; tmp]; %#ok<AGROW>
                        end

                        validRatios = repRatios(isfinite(repRatios));

                        if isempty(validRatios)
                            continue
                        end

                        bar(xpos, mean(validRatios), ...
                            'FaceColor', thisColor, ...
                            'FaceAlpha', 0.35, ...
                            'EdgeColor', 'none');

                        for rr = 1:numel(reps)

                            if ~isfinite(repRatios(rr))
                                continue
                            end

                            scatter(xpos + jitterAmount*randn, ...
                                    repRatios(rr), ...
                                    60, ...
                                    repColors(mod(reps(rr)-1,size(repColors,1))+1,:), ...
                                    'filled', ...
                                    'MarkerFaceAlpha', 0.9, ...
                                    'MarkerEdgeColor', 'none');
                        end
                    end

                    xlim([0.5 2.5])
                    set(gca,'XTick',[1 2], ...
                             'XTickLabel',{cellInfo.NonGreenLabel, cellInfo.GreenLabel})

                    ylabel(sprintf('Ratio: > %.4g / <= %.4g', ...
                                   thisRatioThreshold, thisRatioThreshold))

                    title(sprintf('%s | %s | threshold ratio', ...
                                  condTitle, ratioFolder), ...
                                  'Interpreter','none')

                    box on
                    grid on
                    hold off

                    outPath = fullfile(chanMetricDir, ratioFolder, 'ThresholdRatio', ...
                        sprintf('%s_ThresholdRatio_%s_%s.png', ...
                        groupName, thisChannel, metricShort{m}));

                    savePlotCSV(ratioData, localCsvPath(outPath));
                    saveFig(fig, outPath)
                end
            end

            %% ---------- DISTANCE PLOTS: ALL + FILTERED ----------
            distanceTables = {T, Tfilt};
            distanceNames = {'AllCells','Filtered'};

            for dset = 1:2

                Tdist = distanceTables{dset};
                distFolder = distanceNames{dset};

                if isempty(Tdist)
                    continue
                end

                greenDAll = Tdist(Tdist.IsGreen,:);
                nongreenDAll = Tdist(Tdist.IsNonGreen,:);

                %% ---------- DISTANCE SCATTER ----------
                % Raw look at metric (log y) vs. MinDistMicro (linear x),
                % colored by replicate, before any distance binning.
                fig = figure('Visible','off');
                plotData = table();

                subplot(1,2,1); hold on
                if ~isempty(nongreenDAll)
                    reps = unique(nongreenDAll.Replicate);
                    for rr = 1:numel(reps)
                        thisRep = reps(rr);
                        idx = nongreenDAll.Replicate == thisRep;
                        scatter(nongreenDAll.MinDistMicro(idx), nongreenDAll.MetricValue(idx), pointSize, ...
                            'MarkerFaceColor', repColors(mod(thisRep-1,size(repColors,1))+1,:), ...
                            'MarkerEdgeColor','none', 'MarkerFaceAlpha', pointAlpha);
                    end
                    plotData = [plotData; makePlotData(nongreenDAll,groupName,thisChannel,metricShort{m},distFolder,cellInfo.NonGreenLabel,'NonGreen',1,'')]; %#ok<AGROW>
                end
                set(gca,'YScale','log')
                ylim(safeLogLimits(nongreenDAll.MetricValue))
                xlabel('MinDistMicro')
                ylabel(yLab)
                title(sprintf('%s | %s | %s', condTitle, distFolder, cellInfo.NonGreenLabel), 'Interpreter','none')
                box on; grid on

                subplot(1,2,2); hold on
                if ~isempty(greenDAll)
                    reps = unique(greenDAll.Replicate);
                    for rr = 1:numel(reps)
                        thisRep = reps(rr);
                        idx = greenDAll.Replicate == thisRep;
                        scatter(greenDAll.MinDistMicro(idx), greenDAll.MetricValue(idx), pointSize, ...
                            'MarkerFaceColor', repColors(mod(thisRep-1,size(repColors,1))+1,:), ...
                            'MarkerEdgeColor','none', 'MarkerFaceAlpha', pointAlpha);
                    end
                    plotData = [plotData; makePlotData(greenDAll,groupName,thisChannel,metricShort{m},distFolder,cellInfo.GreenLabel,'Green',2,'')]; %#ok<AGROW>
                end
                set(gca,'YScale','log')
                ylim(safeLogLimits(greenDAll.MetricValue))
                xlabel('MinDistMicro')
                ylabel(yLab)
                title(sprintf('%s | %s | %s', condTitle, distFolder, cellInfo.GreenLabel), 'Interpreter','none')
                box on; grid on

                outPath = fullfile(chanMetricDir,distFolder,'DistanceScatter',sprintf('%s_DistanceScatter_%s_%s.png', groupName, thisChannel, metricShort{m}));
                savePlotCSV(plotData, localCsvPath(outPath));
                saveFig(fig, outPath)

                %% ---------- DISTANCE BIN PLOTS ----------
                % Repeats the box/violin/(FarRed) threshold-ratio plot
                % types from above, but computed separately within each
                % bin of MinDistMicro, for both binning schemes
                % (distBinSetNames). Each bin gets a NonGreen bar and a
                % Green bar placed 1 unit apart, with 3 units between
                % consecutive bins (xCenters below).
                for bset = 1:numel(distBinSetNames)

                    thisBinSet = distBinSetNames{bset};
                    d = Tdist.MinDistMicro;

                    if strcmp(thisBinSet,'2Bins')
                        distNames = distBinNames2;
                        distPretty = distBinPretty2;
                        distSets = cell(1,numel(distNames));

                        distSets{1} = Tdist(d >= 0 & d <= 500, :);
                        distSets{2} = Tdist(d > 500, :);

                    else
                        distNames = distBinNames100;
                        distPretty = distBinPretty100;
                        distSets = cell(1,numel(distNames));

                        distSets{1} = Tdist(d >= 0 & d < 100, :);
                        distSets{2} = Tdist(d >= 100 & d < 200, :);
                        distSets{3} = Tdist(d >= 200 & d < 300, :);
                        distSets{4} = Tdist(d >= 300 & d < 400, :);
                        distSets{5} = Tdist(d >= 400 & d <= 500, :);
                        distSets{6} = Tdist(d > 500, :);
                    end
                                        nBins = numel(distSets);
                    xCenters = zeros(1,nBins*2);
                    xLabels = cell(1,nBins*2);
                    for db = 1:nBins
                        xCenters((db-1)*2+1) = (db-1)*3 + 1;
                        xCenters((db-1)*2+2) = (db-1)*3 + 2;
                        xLabels{(db-1)*2+1} = ['NG ' distPretty{db}];
                        xLabels{(db-1)*2+2} = ['G ' distPretty{db}];
                    end

                    %% ---------- DISTANCE BIN BOX ----------
                    fig = figure('Visible','off','Position',[100 100 1800 650]);
                    hold on
                    plotData = table();

                    for db = 1:nBins

                        Tdb = distSets{db};
                        ng = Tdb(Tdb.IsNonGreen,:);
                        g = Tdb(Tdb.IsGreen,:);
                        xng = xCenters((db-1)*2+1);
                        xg  = xCenters((db-1)*2+2);

                        if ~isempty(ng)
                            addReplicateSwarm(xng, ng.MetricValue, ng.Replicate, repColors, pointSize, pointAlpha, jitterAmount);
                            makeSimpleBox(xng, ng.MetricValue, [0 0 0]);
                            plotData = [plotData; makePlotData(ng,groupName,thisChannel,metricShort{m},distFolder,cellInfo.NonGreenLabel,['NonGreen_' distNames{db}],xng,distPretty{db})]; %#ok<AGROW>
                        end
                        if ~isempty(g)
                            addReplicateSwarm(xg, g.MetricValue, g.Replicate, repColors, pointSize, pointAlpha, jitterAmount);
                            makeSimpleBox(xg, g.MetricValue, [0 0.6 0]);
                            plotData = [plotData; makePlotData(g,groupName,thisChannel,metricShort{m},distFolder,cellInfo.GreenLabel,['Green_' distNames{db}],xg,distPretty{db})]; %#ok<AGROW>
                        end
                    end

                    allGroupedVals = [];
                    for db = 1:nBins
                        allGroupedVals = [allGroupedVals; distSets{db}.MetricValue]; %#ok<AGROW>
                    end
                    if ~isempty(allGroupedVals)
                        set(gca,'YScale','log')
                        ylim(safeLogLimits(allGroupedVals))
                    end

                    set(gca,'XTick',xCenters,'XTickLabel',xLabels)
                    xtickangle(45)
                    ylabel(yLab)
                    title(sprintf('%s | %s | %s distance bins', condTitle, distFolder, thisBinSet), 'Interpreter','none')
                    box on; grid on
                    hold off

                    outPath = fullfile(chanMetricDir,distFolder,'DistBins',thisBinSet,'Box', ...
                        sprintf('%s_Box_%s_%s_%s.png', groupName, thisBinSet, thisChannel, metricShort{m}));

                    csvPath = fullfile(chanMetricDir,distFolder,'DistBins',thisBinSet, ...
                        sprintf('%s_DistBins_%s_%s_%s.csv', groupName, thisBinSet, thisChannel, metricShort{m}));

                    savePlotCSV(plotData, csvPath);
                    saveFig(fig, outPath)

                    %% ---------- DISTANCE BIN VIOLIN ----------
                    fig = figure('Visible','off','Position',[100 100 1800 650]);
                    hold on

                    for db = 1:nBins

                        Tdb = distSets{db};
                        ng = Tdb(Tdb.IsNonGreen,:);
                        g = Tdb(Tdb.IsGreen,:);
                        xng = xCenters((db-1)*2+1);
                        xg  = xCenters((db-1)*2+2);

                        if ~isempty(ng)
                            Violin({ng.MetricValue}, xng, 'HalfViolin', 'full', 'QuartileStyle', 'boxplot', 'DataStyle', 'scatter', 'ShowMean', false, 'ShowMedian', true, 'ShowNotches', false, 'ViolinColor', {[0 0 0]});
                        end
                        if ~isempty(g)
                            Violin({g.MetricValue}, xg, 'HalfViolin', 'full', 'QuartileStyle', 'boxplot', 'DataStyle', 'scatter', 'ShowMean', false, 'ShowMedian', true, 'ShowNotches', false, 'ViolinColor', {[0 0.6 0]});
                        end
                    end

                    if ~isempty(allGroupedVals)
                        set(gca,'YScale','log')
                        ylim(safeLogLimits(allGroupedVals))
                    end

                    set(gca,'XTick',xCenters,'XTickLabel',xLabels)
                    xtickangle(45)
                    ylabel(yLab)
                    title(sprintf('%s | %s | %s distance bins', condTitle, distFolder, thisBinSet), 'Interpreter','none')
                    box on; grid on
                    hold off

                    outPath = fullfile(chanMetricDir,distFolder,'DistBins',thisBinSet,'Violin', ...
                        sprintf('%s_Violin_%s_%s_%s.png', groupName, thisBinSet, thisChannel, metricShort{m}));

                    saveFig(fig, outPath)

                    %% ---------- DISTANCE BIN THRESHOLD RATIO: FARRED ONLY ----------
                    if strcmp(thisChannel,'FarRed')

                        ratioData = table();
                        fig = figure('Visible','off','Position',[100 100 1800 650]);
                        hold on

                        for db = 1:nBins

                            Tdb = distSets{db};

                            if isempty(Tdb)
                                continue
                            end

                            cellGroups = {
                                Tdb(Tdb.IsNonGreen,:), cellInfo.NonGreenLabel, xCenters((db-1)*2+1), [0 0 0], 'NonGreen';
                                Tdb(Tdb.IsGreen,:),    cellInfo.GreenLabel,    xCenters((db-1)*2+2), [0 0.6 0], 'Green'
                            };

                            for cg = 1:size(cellGroups,1)

                                Tcg = cellGroups{cg,1};
                                classLabel = cellGroups{cg,2};
                                xpos = cellGroups{cg,3};
                                thisColor = cellGroups{cg,4};
                                className = cellGroups{cg,5};

                                if isempty(Tcg)
                                    continue
                                end

                                reps = unique(Tcg.Replicate);
                                repRatios = nan(numel(reps),1);

                                for rr = 1:numel(reps)

                                    thisRep = reps(rr);
                                    Trep = Tcg(Tcg.Replicate == thisRep,:);
                                    vals = Trep.MetricValue;

                                    nAbove = sum(vals > thisRatioThreshold);
                                    nBelow = sum(vals <= thisRatioThreshold);

                                    if nBelow == 0
                                        ratioVal = NaN;
                                    else
                                        ratioVal = nAbove / nBelow;
                                    end

                                    repRatios(rr) = ratioVal;

                                    tmp = localMakeRatioTable(groupName, thisChannel, metricShort{m}, ...
                                        distFolder, classLabel, string(className), thisRep, ...
                                        thisRatioThreshold, nAbove, nBelow, ratioVal);
                                    tmp.DistanceBin = string(distPretty{db});
                                    tmp.DistanceBinName = string(distNames{db});
                                    tmp.XPos = xpos;
                                    tmp.DistBinSet = string(thisBinSet);

                                    ratioData = [ratioData; tmp]; %#ok<AGROW>
                                end

                                validRatios = repRatios(isfinite(repRatios));

                                if isempty(validRatios)
                                    continue
                                end

                                bar(xpos, mean(validRatios), ...
                                    'FaceColor', thisColor, ...
                                    'FaceAlpha', 0.35, ...
                                    'EdgeColor', 'none');

                                for rr = 1:numel(reps)

                                    if ~isfinite(repRatios(rr))
                                        continue
                                    end

                                    scatter(xpos + jitterAmount*randn, ...
                                            repRatios(rr), ...
                                            60, ...
                                            repColors(mod(reps(rr)-1,size(repColors,1))+1,:), ...
                                            'filled', ...
                                            'MarkerFaceAlpha', 0.9, ...
                                            'MarkerEdgeColor', 'none');
                                end
                            end
                        end

                        xlim([0.5 max(xCenters)+0.5])
                        set(gca,'XTick',xCenters,'XTickLabel',xLabels)
                        xtickangle(45)

                        ylabel(sprintf('Ratio: > %.4g / <= %.4g', ...
                                       thisRatioThreshold, thisRatioThreshold))

                        title(sprintf('%s | %s | %s threshold ratio', ...
                                      condTitle, distFolder, thisBinSet), ...
                                      'Interpreter','none')

                        box on
                        grid on
                        hold off

                        outPath = fullfile(chanMetricDir,distFolder,'DistBins',thisBinSet,'ThresholdRatio', ...
                            sprintf('%s_ThresholdRatio_%s_%s_%s.png', groupName, thisBinSet, thisChannel, metricShort{m}));

                        savePlotCSV(ratioData, localCsvPath(outPath));
                        saveFig(fig, outPath)
                    end
                end
            end
        end
    end
end

fprintf('\nDone.\n')

%% ================= LOCAL FUNCTIONS =================

function probeInfo = localGetProbeLabels(groupName)
% Maps a condition's group name to the actual smFISH probe gene names
% imaged in the Red (561) and FarRed (640) channels for that probe set --
% the channels are the same physical detectors across conditions, but
% which gene each one's probe targets changes by experiment.
probeInfo = struct('Red','Red','FarRed','FarRed');
if contains(groupName,'vpsLhapR')
    probeInfo.Red = 'hapR';
    probeInfo.FarRed = 'vpsL';
elseif contains(groupName,'vpsLgyrA')
    probeInfo.Red = 'gyrA';
    probeInfo.FarRed = 'vpsL';
end
end

function cellInfo = localGetCellLabels(groupName)
% Maps a condition's group name to Green/NonGreen strain labels. This is
% what lets one script serve both experiment types QuantifyPunctaCells-
% Microcolonies.m supports: "Coculture" groups always pair WT (NonGreen)
% against a DeltamshA mutant (Green); "SpikeIn" groups are always WT
% NonGreen, with the Green label depending on which spike-in strain code
% ("63" -> DeltamshA, "340" -> WT) is embedded in the group name.
deltaStr = char(916);
cellInfo = struct('GreenLabel','Green','NonGreenLabel','NonGreen');
if contains(groupName,'Coculture')
    cellInfo.GreenLabel = [deltaStr 'mshA'];
    cellInfo.NonGreenLabel = 'WT';
elseif contains(groupName,'SpikeIn')
    cellInfo.NonGreenLabel = 'WT';
    if contains(groupName,'63')
        cellInfo.GreenLabel = [deltaStr 'mshA'];
    elseif contains(groupName,'340')
        cellInfo.GreenLabel = 'WT';
    end
end
end

function metricVals = localMetricVecFromStruct(T, fieldName, idx)
% Pulls out column idx (1-4, matching metricShort) of a SumRed/SumFarRed
% field. Handles two possible storage shapes: a plain N-by-4 numeric
% matrix (the normal case, one row per cell), or a cell array of
% row-vectors -- struct2table falls back to a cell column whenever not
% every cell's SumRed/SumFarRed row has the same width, which happens for
% any cell with literally zero matching puncta (sum of a 0-by-4 slice
% still returns a 1x4 zero row, but sum of a fully empty 0x0 matrix -- no
% puncta detected in that channel/FOV at all -- returns a 1x0 row instead).
raw = T.(fieldName);
if isnumeric(raw)
    if size(raw,2) < idx
        error('%s does not have column %d.', fieldName, idx)
    end
    metricVals = raw(:,idx);
elseif iscell(raw)
    metricVals = nan(height(T),1);
    for ii = 1:height(T)
        v = raw{ii};
        if isempty(v)
            metricVals(ii) = NaN;
        else
            metricVals(ii) = v(idx);
        end
    end
else
    error('Unsupported data type for %s.', fieldName)
end
end

function T = localApplyMetricZeroFilter(T, metricVals, doExclude, thr)
% Always drops NaN metric values; additionally drops values at-or-below
% thr when doExclude is true (per-channel-per-metric toggle in USER
% SETTINGS above).
keep = ~isnan(metricVals);
if doExclude
    keep = keep & metricVals > thr;
end
T = T(keep,:);
end

function localAddReplicateSwarm(xBase,yVals,repVals,repColors,pointSize,pointAlpha,jitterAmount)
% Draws one jittered swarmchart per replicate at xBase, colored by
% replicate (so overlapping box/violin plots can still show per-
% replicate spread and identity at a glance).
reps = unique(repVals);
for rr = 1:numel(reps)
    idx = repVals == reps(rr);
    thisX = xBase + jitterAmount*randn(sum(idx),1);
    thisColor = repColors(mod(reps(rr)-1,size(repColors,1))+1,:);
    swarmchart(thisX, yVals(idx), pointSize, thisColor, 'filled', 'MarkerFaceAlpha', pointAlpha, 'MarkerEdgeAlpha', pointAlpha);
end
end

function localMakeSimpleBox(xPos,yVals,colorIn)
% Hand-drawn Tukey box-and-whisker at one x position (median line, IQR
% box, 1.5*IQR whiskers clipped to the actual data range, small end
% caps) -- used instead of MATLAB's built-in boxplot so it can share an
% arbitrary numeric x-axis layout with the swarm points and violins
% plotted alongside it elsewhere in this script.
if isempty(yVals); return; end
q1 = prctile(yVals,25);
q2 = median(yVals);
q3 = prctile(yVals,75);
iqrV = q3-q1;
lowW = max(min(yVals), q1 - 1.5*iqrV);
highW = min(max(yVals), q3 + 1.5*iqrV);
w = 0.25;
plot([xPos xPos],[lowW q1],'-','Color',colorIn,'LineWidth',1.2)
plot([xPos xPos],[q3 highW],'-','Color',colorIn,'LineWidth',1.2)
rectangle('Position',[xPos-w/2 q1 w q3-q1],'EdgeColor',colorIn,'LineWidth',1.2)
plot([xPos-w/2 xPos+w/2],[q2 q2],'-','Color',colorIn,'LineWidth',1.4)
plot([xPos-0.08 xPos+0.08],[lowW lowW],'-','Color',colorIn,'LineWidth',1.2)
plot([xPos-0.08 xPos+0.08],[highW highW],'-','Color',colorIn,'LineWidth',1.2)
end

function ylims = localSafeLogLimits(vals)
% Y-limits for a log-scale axis from strictly-positive finite values
% only, with a fallback range if none exist and padding if min==max
% (a perfectly constant-valued series would otherwise collapse the axis).
vals = vals(:);
vals = vals(~isnan(vals) & vals > 0);
if isempty(vals)
    ylims = [1e-3 1];
else
    ylims = [min(vals)*0.8, max(vals)*1.2];
    if ylims(1) == ylims(2)
        ylims = [ylims(1)*0.8, ylims(2)*1.2];
    end
end
end

function localSaveFig(figHandle,filePath)
% Defensive figure-save wrapper for a long unattended batch run: ensures
% the destination folder exists, exports at 300 DPI, and warns (rather
% than crashing the whole run) if a particular figure fails to save.
[folderPath,~,~] = fileparts(filePath);
if ~exist(folderPath,'dir'); mkdir(folderPath); end
drawnow
try
    exportgraphics(figHandle, filePath, 'Resolution', 300);
catch ME
    warning('Could not save figure: %s\nReason: %s', filePath, ME.message);
end
close(figHandle)
end

function csvPath = localCsvPath(figPath)
% Derives a sibling ".csv" path from a ".png" figure path (same folder,
% same base name) for the many plot types whose CSV export lives right
% next to the PNG it was plotted from.
[folderPath,baseName,~] = fileparts(figPath);
csvPath = fullfile(folderPath, baseName + ".csv");
end

function localSavePlotCSV(T,filePath)
% Writes T to CSV, substituting an empty table if T was empty, and
% ensures the destination folder exists first.
[folderPath,~,~] = fileparts(filePath);
if ~exist(folderPath,'dir'); mkdir(folderPath); end
if isempty(T)
    T = table();
end
writetable(T,filePath);
end

function Tc = localMakeCDFData(xCDF,fCDF,groupName,thisChannel,metricName,analysisName,cellClass,plotGroup,distanceBin)
% Flattens one empirical CDF curve's (x,f) points into a long-format
% table with full bookkeeping columns, for the CDF-by-replicate plots'
% CSV export.
n = numel(xCDF);
Tc = table();
Tc.GroupName = repmat(string(groupName),n,1);
Tc.Channel = repmat(string(thisChannel),n,1);
Tc.Metric = repmat(string(metricName),n,1);
Tc.AnalysisSet = repmat(string(analysisName),n,1);
Tc.CellClass = repmat(string(cellClass),n,1);
Tc.PlotGroup = repmat(string(plotGroup),n,1);
Tc.DistanceBin = repmat(string(distanceBin),n,1);
Tc.CDF_X = xCDF(:);
Tc.CDF_F = fCDF(:);
end

function Tout = localMakePlotData(T,groupName,thisChannel,metricName,analysisName,cellClass,plotGroup,xPos,distanceBin)
% General-purpose "flatten this filtered cell table into a bookkeeping-
% annotated long-format CSV row set" helper, used by nearly every plot
% type in this script. Each data column is only added if it actually
% exists in T, since the tables passed in from different call sites
% (AllCells/Filtered/distance-binned/etc.) don't all carry identical
% column sets.
n = height(T);
Tout = table();
Tout.GroupName = repmat(string(groupName),n,1);
Tout.Channel = repmat(string(thisChannel),n,1);
Tout.Metric = repmat(string(metricName),n,1);
Tout.AnalysisSet = repmat(string(analysisName),n,1);
Tout.CellClass = repmat(string(cellClass),n,1);
Tout.PlotGroup = repmat(string(plotGroup),n,1);
Tout.XPos = repmat(xPos,n,1);
Tout.DistanceBin = repmat(string(distanceBin),n,1);

if ismember('Replicate',T.Properties.VariableNames); Tout.Replicate = T.Replicate; end
if ismember('XYPos',T.Properties.VariableNames); Tout.XYPos = string(T.XYPos); end
if ismember('MetricValue',T.Properties.VariableNames); Tout.MetricValue = T.MetricValue; end
if ismember('MinDistMicro',T.Properties.VariableNames); Tout.MinDistMicro = T.MinDistMicro; end
if ismember('Area',T.Properties.VariableNames); Tout.Area = T.Area; end
if ismember('AspectRatio',T.Properties.VariableNames); Tout.AspectRatio = T.AspectRatio; end
if ismember('WasSplit',T.Properties.VariableNames); Tout.WasSplit = T.WasSplit; end
if ismember('FromBingr',T.Properties.VariableNames); Tout.FromBingr = T.FromBingr; end
if ismember('FromBinbl',T.Properties.VariableNames); Tout.FromBinbl = T.FromBinbl; end
if ismember('IsGreen',T.Properties.VariableNames); Tout.IsGreen = T.IsGreen; end
if ismember('IsNonGreen',T.Properties.VariableNames); Tout.IsNonGreen = T.IsNonGreen; end
if ismember('MeanDistMicro',T.Properties.VariableNames); Tout.MeanDistMicro = T.MeanDistMicro; end
if ismember('NearestMicroInt',T.Properties.VariableNames); Tout.NearestMicroInt = T.NearestMicroInt; end
if ismember('NRed',T.Properties.VariableNames); Tout.NRed = T.NRed; end
if ismember('NFarRed',T.Properties.VariableNames); Tout.NFarRed = T.NFarRed; end
end

function T = localMakeRatioTable(groupName,thisChannel,metricName,analysisName,cellClass,cellClassShort,replicate,thr,nAbove,nBelow,ratioVal)
% One summary row per replicate per cell class per condition, for the
% threshold-ratio plots' CSV export.
T = table();
T.GroupName = string(groupName);
T.Channel = string(thisChannel);
T.Metric = string(metricName);
T.AnalysisSet = string(analysisName);
T.CellClass = string(cellClass);
T.CellClassShort = string(cellClassShort);
T.Replicate = replicate;
T.Threshold = thr;
T.NumAboveThreshold = nAbove;
T.NumBelowOrEqualThreshold = nBelow;
T.RatioAboveVsBelow = ratioVal;
end
