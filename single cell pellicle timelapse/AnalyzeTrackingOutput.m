%% AnalyzeTrackingOutput.m
% Author: Hasreet Gill, Bassler Lab, HHMI/Princeton University
% Developed/revised with Claude (Sonnet 5, Anthropic) - AI-assisted code review and comments.

% Batch driver. Prompts for a PARENT folder containing multiple experiment
% subfolders. For each subfolder found:
%   parent / <subfolder> / Analysis / Tracking_Output
% it deletes any existing Analysis_Output (fresh run every time) and calls
% AnalyzeOneExperiment (a local function below), which does the full
% single-experiment analysis (bead exclusion, Microcolony/Recruited/
% Transient classification, and plots A-Q) exactly as before.
%
% Once every experiment folder has been processed, this driver averages
% two groups of experiments together (group1FolderNames below vs. every
% OTHER processed folder). Before averaging, each experiment's time axis
% is shifted per folderTimeShiftMinutes (so replicates that started
% imaging at different times relative to a shared reference point line up
% correctly). The aligned time axis is NOT cropped to any fixed window --
% it spans the full union of every replicate's own (shifted) time range,
% so no data is dropped from the averaged output. For each of
% plots A/B/C-D/E/F/G-H/L/M/N-O this produces an averaged plot with SEM
% error bars (one error bar per curve for the curve plots, one error bar
% per stacked-bar segment for the stacked-bar plots) plus an Excel
% workbook per group, with one sheet per metric containing each
% replicate's own value, then Mean, then SEM, per category. Plots I/J/K/P/Q
% (microcolony/recruited track-length distribution and volume growth) are
% per-track, not per-frame-per-category, so they aren't averaged here --
% only regenerated per-experiment inside AnalyzeOneExperiment.

clear
clc
close all

%% ===================== PARAMETERS =====================

% Folders in this list are averaged together as "Group 1". Every OTHER
% successfully-processed experiment folder is averaged together as
% "Group 2" -- with exactly 4 experiment folders total, this gives two
% groups of two, as requested.
group1FolderNames = ["20260226", "20260227"];

% Per-folder time alignment: minutes to SUBTRACT from that folder's own
% native time so it lines up with a shared reference clock across
% replicates (e.g. if that experiment's imaging started later relative to
% a shared biological reference point than its partner replicate).
% Folders not listed here default to a shift of 0 (used as-is).
%   20260227's native t=50 -> aligned t=0 (matches 20260226 as-is)
%   20260306's native t=50 -> aligned t=0 (matches 20260320 as-is)
% The aligned time axis for each group is NOT cropped to any fixed window
% -- it automatically spans the full union of every replicate's own
% (shifted) time range, so no data is dropped from the averaged output.
% Aligned time points where only one replicate has data still show that
% replicate's value as the "mean" (SEM is left blank there, since spread
% can't be estimated from a single value).
folderTimeShiftMinutes = containers.Map({'20260227', '20260306'}, [50, 50]);

%% ===================== Discover and process each experiment =====================

parentFolder = uigetdir([], "Select the parent folder containing experiment subfolders");

if parentFolder == 0
    error("No folder selected.")
end

% Only process subfolders that look like a date (e.g. "20260227" --
% exactly 8 digits). This automatically skips "." / ".." and any other
% stray folders (Average_Output from a previous run, notes folders, etc.)
% sitting alongside the experiment folders in the parent.
subDirs = dir(parentFolder);
subDirs = subDirs([subDirs.isdir]);

isDateFolder = false(length(subDirs), 1);
for k = 1:length(subDirs)
    isDateFolder(k) = ~isempty(regexp(subDirs(k).name, '^\d{8}$', 'once'));
end
subDirs = subDirs(isDateFolder);

if isempty(subDirs)
    error("No date-named (8-digit) subfolders found in " + parentFolder)
end

allResults = {};
folderNames = strings(0, 1);

for k = 1:length(subDirs)

    folderName = string(subDirs(k).name);
    experimentAnalysisFolder = fullfile(parentFolder, folderName, "Analysis");

    if ~exist(experimentAnalysisFolder, "dir")
        disp("Skipping " + folderName + ": no Analysis subfolder found.")
        continue
    end

    analysisOutFolder = fullfile(experimentAnalysisFolder, "Analysis_Output");

    if exist(analysisOutFolder, "dir")
        disp("Deleting existing Analysis_Output in " + folderName)
        rmdir(analysisOutFolder, "s")
    end

    disp("===================================================")
    disp("Processing " + folderName)
    disp("===================================================")

    try
        thisResults = AnalyzeOneExperiment(experimentAnalysisFolder);
        allResults{end+1} = thisResults; %#ok<AGROW>
        folderNames(end+1, 1) = folderName; %#ok<AGROW>
    catch ME
        disp("ERROR processing " + folderName + ": " + ME.message)
    end

    close all
end

if isempty(folderNames)
    error("No experiments were successfully processed.")
end

%% ===================== Group averages =====================

averageOutFolder = fullfile(parentFolder, "Average_Output");

if ~exist(averageOutFolder, "dir")
    mkdir(averageOutFolder)
end

isGroup1 = ismember(folderNames, group1FolderNames);

disp("Group 1 (" + strjoin(folderNames(isGroup1), ", ") + ")")
disp("Group 2 (" + strjoin(folderNames(~isGroup1), ", ") + ")")

if nnz(isGroup1) < 2
    disp("WARNING: fewer than 2 folders matched group1FolderNames.")
end
if nnz(~isGroup1) < 2
    disp("WARNING: fewer than 2 folders in group 2 (everything not in group1FolderNames).")
end

averageGroup(allResults(isGroup1), folderNames(isGroup1), "Group1", averageOutFolder, folderTimeShiftMinutes)
averageGroup(allResults(~isGroup1), folderNames(~isGroup1), "Group2", averageOutFolder, folderTimeShiftMinutes)

disp("All experiments processed. Averages written to:")
disp(averageOutFolder)

%% ===================== Helper functions =====================

function results = AnalyzeOneExperiment(analysisFolder)
% Runs the tracking analysis for one experiment: bead exclusion,
% Microcolony/Recruited/Transient classification, and plots A-L (saved
% into analysisFolder/Analysis_Output). Returns the per-frame category
% matrices needed to average multiple experiments together.
%
% analysisFolder must contain a "Tracking_Output" subfolder (the output of
% the tracking script).
%
%   Categories (mutually exclusive and exhaustive -- every non-bead track
%   lands in exactly one):
%     - Microcolony: grows to >= microcolonyGrowthFactor x its initial
%       volume at some point, regardless of initial volume or start frame.
%     - Single recruited cell: not a microcolony, persists >= 2 frames.
%     - Single transient cell: not a microcolony, present for 1 frame only.
%
%   Plots (categories = Microcolony / Single recruited cell / Single
%   transient cell unless noted):
%     A. Volume proportion (stacked bar, normalized to 1)
%     B. Total volume (stacked bar, linear)
%     C/D. Total volume (curve, linear + log)
%     E. Object count proportion (stacked bar, normalized to 1)
%     F. Object count (stacked bar, linear)
%     G/H. Object count (curve, linear + log)
%     I. Microcolony track-length distribution
%     J. Microcolony volume growth, aligned to each track's start frame
%     K. Microcolony volume growth, actual time, colored by duration
%     L. Proportion of objects by cell size (1, 2, ..., 10, >10 cells, stacked bar)
%     M. Proportion of total volume by cell size (same cell-size bins as L,
%        but weighted by volume instead of object count)
%     N/O. New tracks by category (curve, linear + log) -- newly-appearing
%        tracks at each frame, split into all 3 categories
%     P. Recruited volume growth, aligned to each track's start frame
%     Q. Recruited volume growth, actual time, colored by duration
%   Stacked bars are intentionally never shown on a log scale (additive
%   stacking and log axes don't mix -- see curve plots for log versions).

%% ----- PARAMETERS (edit these) -----

% Real-time conversion. Frame 1 is treated as t = 0.
minutesPerFrame = 25;

% Single-cell volume, in voxels (from regionprops3 on a representative
% single-cell object). Used for the bead exclusion below.
singleCellVolume_voxels = 500;

% --- Bead exclusion (applied globally -- excluded objects are removed
% from EVERY plot and every denominator below) ---
% A track is treated as a non-growing "bead" (and excluded) if its first
% surviving frame's volume is below singleCellVolume_voxels, AND EITHER:
%   - it only appears in 1 frame, OR
%   - it appears in >1 frame but never reaches beadGrowthFactor x its
%     first-frame volume.
% Tracks whose first-frame volume is already >= singleCellVolume_voxels
% are not evaluated by this rule at all (kept regardless).
beadGrowthFactor = 2;

% --- Microcolony ---
% Grows at its maximum to >= microcolonyGrowthFactor x its first-frame
% volume, regardless of initial volume or start frame.
microcolonyGrowthFactor = 3;

% --- Single recruited cell ---
% Not a microcolony, and persists >= recruitedMinDurationFrames.
recruitedMinDurationFrames = 2;

% Single transient cell: not a microcolony, and present for exactly 1 frame.
% No extra parameter needed.

% Plots below only show this frame onwards.
plotStartFrame = 1;

disp("Single-cell volume: " + singleCellVolume_voxels + " voxels")

%% ----- Locate and load data -----

outputRoot = fullfile(analysisFolder, "Tracking_Output");

if ~exist(outputRoot, "dir")
    error("Could not find Tracking_Output subfolder in " + analysisFolder)
end

analysisOutFolder = fullfile(analysisFolder, "Analysis_Output");

if ~exist(analysisOutFolder, "dir")
    mkdir(analysisOutFolder)
end

tileDirs = dir(fullfile(outputRoot, "Tile*"));
tileDirs = tileDirs([tileDirs.isdir]);

if isempty(tileDirs)
    error("No Tile* subfolders found in " + outputRoot)
end

objectsAll = table();
lineageAll = table();

for k = 1:length(tileDirs)

    tileFolder = fullfile(tileDirs(k).folder, tileDirs(k).name);

    xlsxFiles = dir(fullfile(tileFolder, "*_TrackedObjects_2D3D.xlsx"));

    if isempty(xlsxFiles)
        disp("WARNING: no TrackedObjects_2D3D.xlsx found in " + tileFolder + ", skipping.")
        continue
    end

    xlsxFile = fullfile(tileFolder, xlsxFiles(1).name);

    disp("Loading " + xlsxFile)

    thisObjects = readtable(xlsxFile, "Sheet", "Objects", "TextType", "string");
    thisLineage = readtable(xlsxFile, "Sheet", "TrackLineage", "TextType", "string");

    objectsAll = [objectsAll; thisObjects]; %#ok<AGROW>
    lineageAll = [lineageAll; thisLineage]; %#ok<AGROW>
end

if height(objectsAll) == 0
    error("No object data loaded from any tile.")
end

% Global track ID, unique across tiles (TrackID is only unique within a tile).
objectsAll.GlobalTrackID = "T" + string(objectsAll.Tile) + "_" + string(objectsAll.TrackID);
lineageAll.GlobalTrackID = "T" + string(lineageAll.LineageTile) + "_" + string(lineageAll.LineageTrackID);

objectsAll.Time_min = (objectsAll.Frame - 1) * minutesPerFrame;

maxFrame = max(objectsAll.Frame);
frames = (1:maxFrame)';
timeVec = (frames - 1) * minutesPerFrame;
plotFrameMask = frames >= plotStartFrame;

% Per-tile last frame, computed from the raw (unfiltered) data -- tiles can
% have different movie lengths, so "did this track get cut off by the
% movie ending" has to be checked against ITS OWN tile's last frame, not
% the global maxFrame across all tiles.
[tileIDsForMaxFrame, ~, tileRowIdxForMaxFrame] = unique(objectsAll.Tile);
maxFrameByTile = accumarray(tileRowIdxForMaxFrame, objectsAll.Frame, [], @max);

disp("Loaded " + height(objectsAll) + " object instances across " + ...
     length(unique(objectsAll.GlobalTrackID)) + " tracks and " + ...
     length(tileDirs) + " tiles.")

%% ----- Bead exclusion -----

% Drop instances with no matched 3D volume -- can't evaluate the bead
% rule, or any category, without a volume.
objectsAll = objectsAll(~isnan(objectsAll.Volume_voxels), :);

[uniqueTrackIDs, ~, trackRowIdx] = unique(objectsAll.GlobalTrackID);
nUniqueTracks = length(uniqueTrackIDs);

minFramePerTrack = accumarray(trackRowIdx, objectsAll.Frame, [nUniqueTracks, 1], @min);
objectsAll.MinSurvivingFrame = minFramePerTrack(trackRowIdx);
isFirstSurvivingRow = objectsAll.Frame == objectsAll.MinSurvivingFrame;

firstFrameVolumeVoxels = nan(nUniqueTracks, 1);
firstFrameVolumeVoxels(trackRowIdx(isFirstSurvivingRow)) = objectsAll.Volume_voxels(isFirstSurvivingRow);

maxVolumeVoxelsPerTrack = accumarray(trackRowIdx, objectsAll.Volume_voxels, [nUniqueTracks, 1], @max, NaN);
maxVolumeUm3PerTrack = accumarray(trackRowIdx, objectsAll.Volume_um3, [nUniqueTracks, 1], @max, NaN);
instanceCountPerTrack = accumarray(trackRowIdx, 1, [nUniqueTracks, 1]);

isBead = (firstFrameVolumeVoxels < singleCellVolume_voxels) & ...
    ((instanceCountPerTrack == 1) | (maxVolumeVoxelsPerTrack < beadGrowthFactor * firstFrameVolumeVoxels));

disp("Excluding " + nnz(isBead) + " / " + nUniqueTracks + " tracks as non-growing beads.")

objectsAll = objectsAll(ismember(objectsAll.GlobalTrackID, uniqueTrackIDs(~isBead)), :);

%% ----- Track classification -----

[~, locLineage] = ismember(uniqueTrackIDs, lineageAll.GlobalTrackID);
startFramePerTrack = lineageAll.StartFrame(locLineage);
durationPerTrack = lineageAll.Duration_frames(locLineage);
endFramePerTrack = lineageAll.EndFrame(locLineage);

tilePerTrack = lineageAll.LineageTile(locLineage);
[~, tileLookupIdx] = ismember(tilePerTrack, tileIDsForMaxFrame);
maxFramePerTrackTile = maxFrameByTile(tileLookupIdx);

isMicrocolony = ~isBead & (maxVolumeVoxelsPerTrack >= microcolonyGrowthFactor * firstFrameVolumeVoxels);

% Reclassify microcolonies that never actually reach a meaningful absolute
% size as recruited instead -- tripling in volume doesn't mean much if the
% peak volume is still tiny.
microcolonyMinPeakVolume_um3 = 150;
isMicrocolonyTooSmall = isMicrocolony & (maxVolumeUm3PerTrack < microcolonyMinPeakVolume_um3);
disp("Reclassifying " + nnz(isMicrocolonyTooSmall) + " microcolonies that never exceed " + ...
     microcolonyMinPeakVolume_um3 + " um^3 as recruited.")

isMicrocolony = isMicrocolony & ~isMicrocolonyTooSmall;

% Reclassify microcolonies that are still present at the last frame of
% THEIR OWN TILE'S movie (tiles can have different lengths) but are
% short-lived (< 5 frames) as recruited instead -- too little track
% history to trust as a genuine microcolony call.
isMicrocolonyShortAtEnd = isMicrocolony & (endFramePerTrack == maxFramePerTrackTile) & (durationPerTrack < 5);
disp("Reclassifying " + nnz(isMicrocolonyShortAtEnd) + " microcolonies present at their tile's last " + ...
     "frame with duration < 5 frames as recruited.")

isMicrocolony = isMicrocolony & ~isMicrocolonyShortAtEnd;

isRecruited = ~isBead & ~isMicrocolony & (durationPerTrack >= recruitedMinDurationFrames);

isTransient = ~isBead & ~isMicrocolony & (durationPerTrack == 1);

% These three are mutually exclusive and exhaustive over the non-bead
% population by construction: a track either grows >=3x (Microcolony), or
% it doesn't -- in which case duration is either == 1 (Transient) or >= 2
% (Recruited). Sanity-check that instead of relying on it silently.
allCatMasks = [isMicrocolony, isRecruited, isTransient];
nOverlap = nnz(sum(allCatMasks, 2) > 1);
nUncategorized = nnz(~isBead & sum(allCatMasks, 2) == 0);

if nOverlap > 0
    disp("WARNING: " + nOverlap + " tracks matched more than one category -- check definitions.")
end
if nUncategorized > 0
    disp("WARNING: " + nUncategorized + " non-bead tracks matched no category -- check definitions.")
end

catLabels = ["Microcolony", "Single recruited cell", "Single transient cell"];
nCats = length(catLabels);

catIdxPerTrack = zeros(nUniqueTracks, 1);   % 0 = uncategorized (should not occur -- see checks above)
catIdxPerTrack(isMicrocolony) = 1;
catIdxPerTrack(isRecruited) = 2;
catIdxPerTrack(isTransient) = 3;

disp("Microcolony tracks: " + nnz(isMicrocolony))
disp("Single recruited cell tracks: " + nnz(isRecruited))
disp("Single transient cell tracks: " + nnz(isTransient))

classificationLabel = repmat("Uncategorized", nUniqueTracks, 1);
classificationLabel(isBead) = "Bead";
classificationLabel(catIdxPerTrack == 1) = catLabels(1);
classificationLabel(catIdxPerTrack == 2) = catLabels(2);
classificationLabel(catIdxPerTrack == 3) = catLabels(3);

writeMixedCSV(fullfile(analysisOutFolder, "TrackClassification.csv"), ...
    ["GlobalTrackID", "StartFrame", "Duration_frames", "FirstFrameVolume_voxels", "MaxVolume_voxels", "Classification"], ...
    {uniqueTrackIDs, startFramePerTrack, durationPerTrack, firstFrameVolumeVoxels, maxVolumeVoxelsPerTrack, classificationLabel});

% Restrict to categorized tracks only for all plots below.
[~, locForRows] = ismember(objectsAll.GlobalTrackID, uniqueTrackIDs);
objectsAll.CatIdx = catIdxPerTrack(locForRows);
objectsAll = objectsAll(objectsAll.CatIdx > 0, :);

[~, frameIdxA] = ismember(objectsAll.Frame, frames);
volMatrix = accumarray([frameIdxA, objectsAll.CatIdx], objectsAll.Volume_um3, [maxFrame, nCats]);
countMatrix = accumarray([frameIdxA, objectsAll.CatIdx], 1, [maxFrame, nCats]);

totalVol = sum(volMatrix, 2);
propVol = volMatrix ./ totalVol;

totalCount = sum(countMatrix, 2);
propCount = countMatrix ./ totalCount;

%% ----- PLOT A: volume proportion by category (stacked bar, normalized) -----

headerA = ["Frame", "Time_min", catLabels];
writeWideCSV(fullfile(analysisOutFolder, "PlotA_VolumeProportion.csv"), ...
    headerA, [frames(plotFrameMask), timeVec(plotFrameMask), propVol(plotFrameMask, :)]);

figure("Name", "Volume proportion by category")
bar(timeVec(plotFrameMask), propVol(plotFrameMask, :), "stacked")
xlabel("Time (min)")
ylabel("Proportion of total volume")
title("Volume proportion by category")
legend(catLabels, "Location", "eastoutside")
savePng(analysisOutFolder, "PlotA_VolumeProportion")

%% ----- PLOT B: total volume by category (stacked bar, linear) -----

headerB = ["Frame", "Time_min", catLabels];
writeWideCSV(fullfile(analysisOutFolder, "PlotB_TotalVolume_StackedBar.csv"), ...
    headerB, [frames(plotFrameMask), timeVec(plotFrameMask), volMatrix(plotFrameMask, :)]);

figure("Name", "Total volume by category (stacked bar)")
bar(timeVec(plotFrameMask), volMatrix(plotFrameMask, :), "stacked")
xlabel("Time (min)")
ylabel("Total volume (um^3)")
title("Total volume by category (stacked bar)")
legend(catLabels, "Location", "eastoutside")
savePng(analysisOutFolder, "PlotB_TotalVolume_StackedBar")

%% ----- PLOTS C/D: total volume by category (curves) -----

headerCD = ["Frame", "Time_min", catLabels];
writeWideCSV(fullfile(analysisOutFolder, "PlotCD_TotalVolume_Curve.csv"), ...
    headerCD, [frames(plotFrameMask), timeVec(plotFrameMask), volMatrix(plotFrameMask, :)]);

figure("Name", "Total volume by category (curve, linear)")
plot(timeVec(plotFrameMask), volMatrix(plotFrameMask, :), "LineWidth", 2)
xlabel("Time (min)")
ylabel("Total volume (um^3)")
title("Total volume by category (curve, linear)")
legend(catLabels, "Location", "best")
savePng(analysisOutFolder, "PlotC_TotalVolume_Curve_Linear")

figure("Name", "Total volume by category (curve, log)")
semilogy(timeVec(plotFrameMask), volMatrix(plotFrameMask, :), "LineWidth", 2)
xlabel("Time (min)")
ylabel("Total volume (um^3, log scale)")
title("Total volume by category (curve, log)")
legend(catLabels, "Location", "best")
savePng(analysisOutFolder, "PlotD_TotalVolume_Curve_Log")

%% ----- PLOT E: object count by category (stacked bar, normalized) -----

headerE = ["Frame", "Time_min", catLabels];
writeWideCSV(fullfile(analysisOutFolder, "PlotE_ObjectCountProportion.csv"), ...
    headerE, [frames(plotFrameMask), timeVec(plotFrameMask), propCount(plotFrameMask, :)]);

figure("Name", "Object count proportion by category")
bar(timeVec(plotFrameMask), propCount(plotFrameMask, :), "stacked")
xlabel("Time (min)")
ylabel("Proportion of object count")
title("Object count proportion by category")
legend(catLabels, "Location", "eastoutside")
savePng(analysisOutFolder, "PlotE_ObjectCountProportion")

%% ----- PLOT F: object count by category (stacked bar, linear) -----

headerF = ["Frame", "Time_min", catLabels];
writeWideCSV(fullfile(analysisOutFolder, "PlotF_ObjectCount_StackedBar.csv"), ...
    headerF, [frames(plotFrameMask), timeVec(plotFrameMask), countMatrix(plotFrameMask, :)]);

figure("Name", "Object count by category (stacked bar)")
bar(timeVec(plotFrameMask), countMatrix(plotFrameMask, :), "stacked")
xlabel("Time (min)")
ylabel("Number of objects")
title("Object count by category (stacked bar)")
legend(catLabels, "Location", "eastoutside")
savePng(analysisOutFolder, "PlotF_ObjectCount_StackedBar")

%% ----- PLOTS G/H: object count by category (curves) -----

headerGH = ["Frame", "Time_min", catLabels];
writeWideCSV(fullfile(analysisOutFolder, "PlotGH_ObjectCount_Curve.csv"), ...
    headerGH, [frames(plotFrameMask), timeVec(plotFrameMask), countMatrix(plotFrameMask, :)]);

figure("Name", "Object count by category (curve, linear)")
plot(timeVec(plotFrameMask), countMatrix(plotFrameMask, :), "LineWidth", 2)
xlabel("Time (min)")
ylabel("Number of objects")
title("Object count by category (curve, linear)")
legend(catLabels, "Location", "best")
savePng(analysisOutFolder, "PlotG_ObjectCount_Curve_Linear")

figure("Name", "Object count by category (curve, log)")
semilogy(timeVec(plotFrameMask), countMatrix(plotFrameMask, :), "LineWidth", 2)
xlabel("Time (min)")
ylabel("Number of objects (log scale)")
title("Object count by category (curve, log)")
legend(catLabels, "Location", "best")
savePng(analysisOutFolder, "PlotH_ObjectCount_Curve_Log")

%% ----- PLOT I: microcolony track-length distribution -----

if nnz(isMicrocolony) == 0

    disp("Skipping Plot I: no microcolony tracks found.")

else

    microcolonyDurations = durationPerTrack(isMicrocolony);
    maxMicrocolonyDuration = max(microcolonyDurations);

    durationBins = (1:maxMicrocolonyDuration)';
    countByDuration = accumarray(microcolonyDurations, 1, [maxMicrocolonyDuration, 1]);
    proportionByDuration = countByDuration / nnz(isMicrocolony);

    writeWideCSV(fullfile(analysisOutFolder, "PlotI_MicrocolonyTrackLength.csv"), ...
        ["Duration_frames", "Count", "Proportion"], [durationBins, countByDuration, proportionByDuration]);

    figure("Name", "Microcolony track length distribution")
    bar(durationBins, proportionByDuration)
    xlabel("Track duration (frames)")
    ylabel("Proportion of microcolonies")
    title("Microcolony track length distribution")
    savePng(analysisOutFolder, "PlotI_MicrocolonyTrackLength")

end

%% ----- PLOT J: microcolony volume growth, aligned to start -----

microcolonyObjects = objectsAll(objectsAll.CatIdx == 1, :);
microcolonyTrackIDs = unique(microcolonyObjects.GlobalTrackID);
nMicrocolonies = length(microcolonyTrackIDs);

disp(nMicrocolonies + " microcolony tracks plotted in Plot J.")

if nMicrocolonies == 0

    disp("Skipping Plot J: no microcolony tracks found.")

else

    [~, trackCol] = ismember(microcolonyObjects.GlobalTrackID, microcolonyTrackIDs);

    startFrameByTrack = nan(nMicrocolonies, 1);

    for t = 1:nMicrocolonies
        lineageRow = find(lineageAll.GlobalTrackID == microcolonyTrackIDs(t), 1);
        startFrameByTrack(t) = lineageAll.StartFrame(lineageRow);
    end

    frameSinceStart = microcolonyObjects.Frame - startFrameByTrack(trackCol);
    maxFrameSinceStart = max(frameSinceStart);
    alignedFrames = (0:maxFrameSinceStart)';
    alignedTime = alignedFrames * minutesPerFrame;

    matrixAligned = nan(length(alignedFrames), nMicrocolonies);
    linIdxAligned = sub2ind(size(matrixAligned), frameSinceStart + 1, trackCol);
    matrixAligned(linIdxAligned) = microcolonyObjects.Volume_um3;

    headerAligned = ["FramesSinceStart", "TimeSinceStart_min", microcolonyTrackIDs'];
    writeWideCSV(fullfile(analysisOutFolder, "PlotJ_MicrocolonyVolumeGrowth.csv"), ...
        headerAligned, [alignedFrames, alignedTime, matrixAligned]);

    cmap = parula(256);
    startRange = max(startFrameByTrack) - min(startFrameByTrack);

    if startRange == 0
        colorIdx = ones(nMicrocolonies, 1);
    else
        normStart = (startFrameByTrack - min(startFrameByTrack)) / startRange;
        colorIdx = round(normStart * (size(cmap, 1) - 1)) + 1;
    end

    figure("Name", "Microcolony volume growth")
    hold on

    for t = 1:nMicrocolonies
        plot(alignedTime, matrixAligned(:, t), "Color", cmap(colorIdx(t), :), "LineWidth", 1.2)
    end

    colormap(cmap)
    caxis([min(startFrameByTrack), max(startFrameByTrack)])
    cb = colorbar;
    cb.Label.String = "Track start frame";

    xlabel("Time since track start (min)")
    ylabel("Volume (um^3)")
    title("Microcolony volume growth")
    savePng(analysisOutFolder, "PlotJ_MicrocolonyVolumeGrowth")

end

%% ----- PLOT K: microcolony volume growth, actual time, colored by duration -----

microcolonyObjectsH = objectsAll(objectsAll.CatIdx == 1, :);
microcolonyTrackIDsH = unique(microcolonyObjectsH.GlobalTrackID);
nMicrocoloniesH = length(microcolonyTrackIDsH);

if nMicrocoloniesH == 0

    disp("Skipping Plot K: no microcolony tracks found.")

else

    [~, trackColH] = ismember(microcolonyObjectsH.GlobalTrackID, microcolonyTrackIDsH);

    durationByTrackH = nan(nMicrocoloniesH, 1);

    for t = 1:nMicrocoloniesH
        lineageRow = find(lineageAll.GlobalTrackID == microcolonyTrackIDsH(t), 1);
        durationByTrackH(t) = lineageAll.Duration_frames(lineageRow);
    end

    [~, frameRowH] = ismember(microcolonyObjectsH.Frame, frames);

    matrixActualH = nan(maxFrame, nMicrocoloniesH);
    linIdxActualH = sub2ind(size(matrixActualH), frameRowH, trackColH);
    matrixActualH(linIdxActualH) = microcolonyObjectsH.Volume_um3;

    headerActualH = ["Frame", "Time_min", microcolonyTrackIDsH'];
    writeWideCSV(fullfile(analysisOutFolder, "PlotK_MicrocolonyVolumeGrowth_ActualTime.csv"), ...
        headerActualH, [frames, timeVec, matrixActualH]);

    cmap = parula(256);
    durRange = max(durationByTrackH) - min(durationByTrackH);

    if durRange == 0
        colorIdxH = ones(nMicrocoloniesH, 1);
    else
        normDur = (durationByTrackH - min(durationByTrackH)) / durRange;
        colorIdxH = round(normDur * (size(cmap, 1) - 1)) + 1;
    end

    figure("Name", "Microcolony volume growth (actual time)")
    hold on

    for t = 1:nMicrocoloniesH
        plot(timeVec, matrixActualH(:, t), "Color", cmap(colorIdxH(t), :), "LineWidth", 1.2)
    end

    colormap(cmap)
    caxis([min(durationByTrackH), max(durationByTrackH)])
    cb = colorbar;
    cb.Label.String = "Track duration (frames)";

    xlabel("Time (min)")
    ylabel("Volume (um^3)")
    title("Microcolony volume growth (actual time)")
    savePng(analysisOutFolder, "PlotK_MicrocolonyVolumeGrowth_ActualTime")

end

%% ----- PLOT L: proportion of objects by cell size (1, 2, ..., 10, >10) -----

nCellsPerInstance = round(objectsAll.Volume_voxels / singleCellVolume_voxels);

cellSizeCatLabels = ["1", "2", "3", "4", "5", "6", "7", "8", "9", "10", ">10"];
nCellSizeCats = length(cellSizeCatLabels);

cellSizeCatIdx = nCellsPerInstance;
cellSizeCatIdx(cellSizeCatIdx < 1) = 1;              % clamp any sub-single-cell instance into bin "1"
cellSizeCatIdx(cellSizeCatIdx > 10) = nCellSizeCats;  % anything above 10 -> ">10" bin

cellSizeCountMatrix = accumarray([frameIdxA, cellSizeCatIdx], 1, [maxFrame, nCellSizeCats]);
totalCellSizeCount = sum(cellSizeCountMatrix, 2);
propCellSize = cellSizeCountMatrix ./ totalCellSizeCount;

headerL = ["Frame", "Time_min", cellSizeCatLabels];
writeWideCSV(fullfile(analysisOutFolder, "PlotL_CellSizeProportion.csv"), ...
    headerL, [frames(plotFrameMask), timeVec(plotFrameMask), propCellSize(plotFrameMask, :)]);

figure("Name", "Proportion of objects by cell size")
bar(timeVec(plotFrameMask), propCellSize(plotFrameMask, :), "stacked")
xlabel("Time (min)")
ylabel("Proportion of objects")
title("Proportion of objects by cell size")
legend(cellSizeCatLabels, "Location", "eastoutside")
savePng(analysisOutFolder, "PlotL_CellSizeProportion")

%% ----- PLOT M: proportion of total volume by cell size (1, 2, ..., 10, >10) -----

cellSizeVolMatrix = accumarray([frameIdxA, cellSizeCatIdx], objectsAll.Volume_um3, [maxFrame, nCellSizeCats]);
totalCellSizeVol = sum(cellSizeVolMatrix, 2);
propCellSizeVol = cellSizeVolMatrix ./ totalCellSizeVol;

headerM = ["Frame", "Time_min", cellSizeCatLabels];
writeWideCSV(fullfile(analysisOutFolder, "PlotM_CellSizeVolumeProportion.csv"), ...
    headerM, [frames(plotFrameMask), timeVec(plotFrameMask), propCellSizeVol(plotFrameMask, :)]);

figure("Name", "Proportion of total volume by cell size")
bar(timeVec(plotFrameMask), propCellSizeVol(plotFrameMask, :), "stacked")
xlabel("Time (min)")
ylabel("Proportion of total volume")
title("Proportion of total volume by cell size")
legend(cellSizeCatLabels, "Location", "eastoutside")
savePng(analysisOutFolder, "PlotM_CellSizeVolumeProportion")

%% ----- PLOTS N/O: newly appearing tracks by category (curves) -----
% New tracks starting at each time point, split into all 3 categories
% (Microcolony start events included, alongside newly-recruited and
% newly-transient tracks).

newTrackCatLabels = catLabels;
nNewTrackCats = nCats;

isNewTrackCounted = catIdxPerTrack > 0;

[~, frameIdxNew] = ismember(startFramePerTrack(isNewTrackCounted), frames);
newTrackMatrix = accumarray([frameIdxNew, catIdxPerTrack(isNewTrackCounted)], 1, [maxFrame, nNewTrackCats]);

headerNO = ["Frame", "Time_min", newTrackCatLabels];
writeWideCSV(fullfile(analysisOutFolder, "PlotNO_NewTracksByCategory_Curve.csv"), ...
    headerNO, [frames(plotFrameMask), timeVec(plotFrameMask), newTrackMatrix(plotFrameMask, :)]);

figure("Name", "New tracks by category (curve, linear)")
plot(timeVec(plotFrameMask), newTrackMatrix(plotFrameMask, :), "LineWidth", 2)
xlabel("Time (min)")
ylabel("Number of new tracks")
title("New tracks by category (curve, linear)")
legend(newTrackCatLabels, "Location", "best")
savePng(analysisOutFolder, "PlotN_NewTracksByCategory_Curve_Linear")

figure("Name", "New tracks by category (curve, log)")
semilogy(timeVec(plotFrameMask), newTrackMatrix(plotFrameMask, :), "LineWidth", 2)
xlabel("Time (min)")
ylabel("Number of new tracks (log scale)")
title("New tracks by category (curve, log)")
legend(newTrackCatLabels, "Location", "best")
savePng(analysisOutFolder, "PlotO_NewTracksByCategory_Curve_Log")

%% ----- PLOT P: recruited volume growth, aligned to start -----

recruitedObjects = objectsAll(objectsAll.CatIdx == 2, :);
recruitedTrackIDs = unique(recruitedObjects.GlobalTrackID);
nRecruited = length(recruitedTrackIDs);

disp(nRecruited + " recruited tracks plotted in Plot P.")

if nRecruited == 0

    disp("Skipping Plot P: no recruited tracks found.")

else

    [~, trackColP] = ismember(recruitedObjects.GlobalTrackID, recruitedTrackIDs);

    startFrameByTrackP = nan(nRecruited, 1);

    for t = 1:nRecruited
        lineageRow = find(lineageAll.GlobalTrackID == recruitedTrackIDs(t), 1);
        startFrameByTrackP(t) = lineageAll.StartFrame(lineageRow);
    end

    frameSinceStartP = recruitedObjects.Frame - startFrameByTrackP(trackColP);
    maxFrameSinceStartP = max(frameSinceStartP);
    alignedFramesP = (0:maxFrameSinceStartP)';
    alignedTimeP = alignedFramesP * minutesPerFrame;

    matrixAlignedP = nan(length(alignedFramesP), nRecruited);
    linIdxAlignedP = sub2ind(size(matrixAlignedP), frameSinceStartP + 1, trackColP);
    matrixAlignedP(linIdxAlignedP) = recruitedObjects.Volume_um3;

    headerAlignedP = ["FramesSinceStart", "TimeSinceStart_min", recruitedTrackIDs'];
    writeWideCSV(fullfile(analysisOutFolder, "PlotP_RecruitedVolumeGrowth.csv"), ...
        headerAlignedP, [alignedFramesP, alignedTimeP, matrixAlignedP]);

    cmap = parula(256);
    startRangeP = max(startFrameByTrackP) - min(startFrameByTrackP);

    if startRangeP == 0
        colorIdxP = ones(nRecruited, 1);
    else
        normStartP = (startFrameByTrackP - min(startFrameByTrackP)) / startRangeP;
        colorIdxP = round(normStartP * (size(cmap, 1) - 1)) + 1;
    end

    figure("Name", "Recruited volume growth")
    hold on

    for t = 1:nRecruited
        plot(alignedTimeP, matrixAlignedP(:, t), "Color", cmap(colorIdxP(t), :), "LineWidth", 1.2)
    end

    colormap(cmap)
    caxis([min(startFrameByTrackP), max(startFrameByTrackP)])
    cb = colorbar;
    cb.Label.String = "Track start frame";

    xlabel("Time since track start (min)")
    ylabel("Volume (um^3)")
    title("Recruited volume growth")
    savePng(analysisOutFolder, "PlotP_RecruitedVolumeGrowth")

end

%% ----- PLOT Q: recruited volume growth, actual time, colored by duration -----

recruitedObjectsQ = objectsAll(objectsAll.CatIdx == 2, :);
recruitedTrackIDsQ = unique(recruitedObjectsQ.GlobalTrackID);
nRecruitedQ = length(recruitedTrackIDsQ);

if nRecruitedQ == 0

    disp("Skipping Plot Q: no recruited tracks found.")

else

    [~, trackColQ] = ismember(recruitedObjectsQ.GlobalTrackID, recruitedTrackIDsQ);

    durationByTrackQ = nan(nRecruitedQ, 1);

    for t = 1:nRecruitedQ
        lineageRow = find(lineageAll.GlobalTrackID == recruitedTrackIDsQ(t), 1);
        durationByTrackQ(t) = lineageAll.Duration_frames(lineageRow);
    end

    [~, frameRowQ] = ismember(recruitedObjectsQ.Frame, frames);

    matrixActualQ = nan(maxFrame, nRecruitedQ);
    linIdxActualQ = sub2ind(size(matrixActualQ), frameRowQ, trackColQ);
    matrixActualQ(linIdxActualQ) = recruitedObjectsQ.Volume_um3;

    headerActualQ = ["Frame", "Time_min", recruitedTrackIDsQ'];
    writeWideCSV(fullfile(analysisOutFolder, "PlotQ_RecruitedVolumeGrowth_ActualTime.csv"), ...
        headerActualQ, [frames, timeVec, matrixActualQ]);

    cmap = parula(256);
    durRangeQ = max(durationByTrackQ) - min(durationByTrackQ);

    if durRangeQ == 0
        colorIdxQ = ones(nRecruitedQ, 1);
    else
        normDurQ = (durationByTrackQ - min(durationByTrackQ)) / durRangeQ;
        colorIdxQ = round(normDurQ * (size(cmap, 1) - 1)) + 1;
    end

    figure("Name", "Recruited volume growth (actual time)")
    hold on

    for t = 1:nRecruitedQ
        plot(timeVec, matrixActualQ(:, t), "Color", cmap(colorIdxQ(t), :), "LineWidth", 1.2)
    end

    colormap(cmap)
    caxis([min(durationByTrackQ), max(durationByTrackQ)])
    cb = colorbar;
    cb.Label.String = "Track duration (frames)";

    xlabel("Time (min)")
    ylabel("Volume (um^3)")
    title("Recruited volume growth (actual time)")
    savePng(analysisOutFolder, "PlotQ_RecruitedVolumeGrowth_ActualTime")

end

disp("All plots and CSVs written to:")
disp(analysisOutFolder)

%% ----- Return results for cross-experiment averaging -----

results.catLabels = catLabels;
results.cellSizeCatLabels = cellSizeCatLabels;
results.newTrackCatLabels = newTrackCatLabels;
results.minutesPerFrame = minutesPerFrame;
results.maxFrame = maxFrame;
results.frames = frames;
results.timeVec = timeVec;
results.volMatrix = volMatrix;
results.countMatrix = countMatrix;
results.propVol = propVol;
results.propCount = propCount;
results.propCellSize = propCellSize;
results.propCellSizeVol = propCellSizeVol;
results.newTrackMatrix = newTrackMatrix;

end

function savePng(folder, name)
% Saves the current figure as a PNG under FOLDER/NAME.png.
    saveas(gcf, fullfile(folder, name + ".png"));
end

function writeWideCSV(filename, headerCells, numericData)
% Writes a header row (as given, not sanitized into MATLAB identifiers)
% followed by a numeric data matrix, using writecell so column headers
% like "Single recruited cell" or track IDs like "T3_15" are preserved exactly.

    outCell = [cellstr(headerCells); num2cell(numericData)];
    writecell(outCell, filename);

end

function writeMixedCSV(filename, header, columns)
% Writes a header row followed by data columns that may be a mix of
% numeric and string column vectors. Avoids the table() constructor
% entirely, since its "VariableNames" name-value syntax is not accepted
% the same way across MATLAB versions.

    nCols = numel(columns);
    nRows = length(columns{1});

    outCell = cell(nRows + 1, nCols);
    outCell(1, :) = cellstr(header);

    for c = 1:nCols

        col = columns{c};

        if isstring(col)
            outCell(2:end, c) = cellstr(col);
        else
            outCell(2:end, c) = num2cell(col);
        end
    end

    writecell(outCell, filename);

end

function averageGroup(resultsList, names, groupLabel, outFolder, folderTimeShiftMinutes)
% Computes mean +/- SEM across the experiments in resultsList, after
% shifting each experiment's time axis per folderTimeShiftMinutes. The
% aligned time axis spans the full UNION of every replicate's own shifted
% time range -- no data is cropped out of the final output. Aligned time
% points covered by only one replicate simply have that replicate's own
% value as the mean, with SEM left blank (can't estimate spread from one
% value). Writes the averaged plots + one Excel workbook for this group.

    if isempty(resultsList)
        disp("Skipping " + groupLabel + ": no experiments in this group.")
        return
    end

    minutesPerFrame = resultsList{1}.minutesPerFrame;
    catLabels = resultsList{1}.catLabels;
    cellSizeCatLabels = resultsList{1}.cellSizeCatLabels;
    newTrackCatLabels = resultsList{1}.newTrackCatLabels;

    nExp = length(resultsList);
    alignedStartMin = zeros(nExp, 1);
    alignedEndMin = zeros(nExp, 1);

    for e = 1:nExp
        name = char(names(e));
        if isKey(folderTimeShiftMinutes, name)
            shiftMin = folderTimeShiftMinutes(name);
        else
            shiftMin = 0;
        end
        alignedStartMin(e) = 0 - shiftMin;
        alignedEndMin(e) = (resultsList{e}.maxFrame - 1) * minutesPerFrame - shiftMin;
    end

    alignedWindowMinutes = [min(alignedStartMin), max(alignedEndMin)];

    disp(groupLabel + ": averaging " + strjoin(names, ", ") + ...
         " (aligned window " + alignedWindowMinutes(1) + "-" + alignedWindowMinutes(2) + ...
         " min, full union -- no data cropped)")

    alignedTimeVec = (alignedWindowMinutes(1):minutesPerFrame:alignedWindowMinutes(2))';

    volMatrix3D = stackAndPadAligned(resultsList, names, "volMatrix", alignedTimeVec, folderTimeShiftMinutes, minutesPerFrame);
    countMatrix3D = stackAndPadAligned(resultsList, names, "countMatrix", alignedTimeVec, folderTimeShiftMinutes, minutesPerFrame);
    propVol3D = stackAndPadAligned(resultsList, names, "propVol", alignedTimeVec, folderTimeShiftMinutes, minutesPerFrame);
    propCount3D = stackAndPadAligned(resultsList, names, "propCount", alignedTimeVec, folderTimeShiftMinutes, minutesPerFrame);
    propCellSize3D = stackAndPadAligned(resultsList, names, "propCellSize", alignedTimeVec, folderTimeShiftMinutes, minutesPerFrame);
    propCellSizeVol3D = stackAndPadAligned(resultsList, names, "propCellSizeVol", alignedTimeVec, folderTimeShiftMinutes, minutesPerFrame);
    newTrackMatrix3D = stackAndPadAligned(resultsList, names, "newTrackMatrix", alignedTimeVec, folderTimeShiftMinutes, minutesPerFrame);

    excelFile = fullfile(outFolder, groupLabel + "_MeanSEM.xlsx");

    if exist(excelFile, "file")
        delete(excelFile)
    end

    plotAndExportStackedBar(propVol3D, catLabels, names, alignedTimeVec, ...
        groupLabel + ": volume proportion", "Proportion of total volume", ...
        outFolder, groupLabel + "_A_VolumeProportion", excelFile, "VolumeProportion")

    plotAndExportStackedBar(volMatrix3D, catLabels, names, alignedTimeVec, ...
        groupLabel + ": total volume (stacked bar)", "Total volume (um^3)", ...
        outFolder, groupLabel + "_B_TotalVolume_StackedBar", excelFile, "TotalVolume")

    plotAndExportCurve(volMatrix3D, catLabels, names, alignedTimeVec, ...
        groupLabel + ": total volume (curve)", "Total volume (um^3)", ...
        outFolder, groupLabel + "_CD_TotalVolume_Curve", excelFile, "TotalVolumeCurve")

    plotAndExportStackedBar(propCount3D, catLabels, names, alignedTimeVec, ...
        groupLabel + ": object count proportion", "Proportion of object count", ...
        outFolder, groupLabel + "_E_ObjectCountProportion", excelFile, "ObjectCountProportion")

    plotAndExportStackedBar(countMatrix3D, catLabels, names, alignedTimeVec, ...
        groupLabel + ": object count (stacked bar)", "Number of objects", ...
        outFolder, groupLabel + "_F_ObjectCount_StackedBar", excelFile, "ObjectCount")

    plotAndExportCurve(countMatrix3D, catLabels, names, alignedTimeVec, ...
        groupLabel + ": object count (curve)", "Number of objects", ...
        outFolder, groupLabel + "_GH_ObjectCount_Curve", excelFile, "ObjectCountCurve")

    plotAndExportStackedBar(propCellSize3D, cellSizeCatLabels, names, alignedTimeVec, ...
        groupLabel + ": cell size proportion", "Proportion of objects", ...
        outFolder, groupLabel + "_L_CellSizeProportion", excelFile, "CellSizeProportion")

    plotAndExportStackedBar(propCellSizeVol3D, cellSizeCatLabels, names, alignedTimeVec, ...
        groupLabel + ": cell size volume proportion", "Proportion of total volume", ...
        outFolder, groupLabel + "_M_CellSizeVolumeProportion", excelFile, "CellSizeVolumeProportion")

    plotAndExportCurve(newTrackMatrix3D, newTrackCatLabels, names, alignedTimeVec, ...
        groupLabel + ": new tracks by category", "Number of new tracks", ...
        outFolder, groupLabel + "_NO_NewTracksByCategory_Curve", excelFile, "NewTracksByCategory")

    % --- Separate workbook: every replicate's own RAW data, indexed by its
    % own native Frame/Time_min, with NO shifting and NO cropping. Purely
    % for reference/cross-checking against the aligned+averaged workbook
    % above -- no averaging is done here. ---
    rawExcelFile = fullfile(outFolder, groupLabel + "_RawReplicates.xlsx");

    if exist(rawExcelFile, "file")
        delete(rawExcelFile)
    end

    writeRawReplicateSheet(resultsList, names, "propVol", catLabels, "VolumeProportion", rawExcelFile, minutesPerFrame)
    writeRawReplicateSheet(resultsList, names, "volMatrix", catLabels, "TotalVolume", rawExcelFile, minutesPerFrame)
    writeRawReplicateSheet(resultsList, names, "volMatrix", catLabels, "TotalVolumeCurve", rawExcelFile, minutesPerFrame)
    writeRawReplicateSheet(resultsList, names, "propCount", catLabels, "ObjectCountProportion", rawExcelFile, minutesPerFrame)
    writeRawReplicateSheet(resultsList, names, "countMatrix", catLabels, "ObjectCount", rawExcelFile, minutesPerFrame)
    writeRawReplicateSheet(resultsList, names, "propCellSize", cellSizeCatLabels, "CellSizeProportion", rawExcelFile, minutesPerFrame)
    writeRawReplicateSheet(resultsList, names, "propCellSizeVol", cellSizeCatLabels, "CellSizeVolumeProportion", rawExcelFile, minutesPerFrame)
    writeRawReplicateSheet(resultsList, names, "newTrackMatrix", newTrackCatLabels, "NewTracksByCategory", rawExcelFile, minutesPerFrame)

end

function writeRawReplicateSheet(resultsList, names, fieldName, catLabels, sheetName, rawExcelFile, minutesPerFrame)
% Writes ONE sheet: every replicate's own value for fieldName, indexed by
% that replicate's own native Frame/Time_min (no shift, no alignment).
% Shorter replicates are padded with NaN only where their own data has
% genuinely ended -- nothing is cropped off the longer replicate(s).

    nExp = length(resultsList);
    nCats = length(catLabels);

    nativeMaxFrame = max(cellfun(@(r) size(r.(fieldName), 1), resultsList));
    nativeFrames = (1:nativeMaxFrame)';
    nativeTimeVec = (nativeFrames - 1) * minutesPerFrame;

    header = ["Frame", "Time_min"];
    dataCols = {nativeFrames, nativeTimeVec};

    for c = 1:nCats
        for e = 1:nExp
            header(end+1) = catLabels(c) + "_" + names(e); %#ok<AGROW>

            thisData = resultsList{e}.(fieldName)(:, c);
            paddedData = nan(nativeMaxFrame, 1);
            paddedData(1:length(thisData)) = thisData;

            dataCols{end+1} = paddedData; %#ok<AGROW>
        end
    end

    nRows = nativeMaxFrame;
    nCols = length(dataCols);

    outCell = cell(nRows + 1, nCols);
    outCell(1, :) = cellstr(header);

    for col = 1:nCols
        outCell(2:end, col) = num2cell(dataCols{col});
    end

    writecell(outCell, rawExcelFile, "Sheet", sheetName)

end

function data3D = stackAndPadAligned(resultsList, names, fieldName, alignedTimeVec, folderTimeShiftMinutes, minutesPerFrame)
% Stacks a (frame x category) matrix from each experiment into a
% (alignedTime x category x experiment) array. Each experiment's own
% native time axis is shifted by folderTimeShiftMinutes(name) (0 if that
% name isn't listed) before being matched up to the shared alignedTimeVec,
% so experiments that started imaging at different times relative to a
% shared reference point line up correctly. Aligned time points that fall
% outside an experiment's own native data range are left as NaN.

    nExp = length(resultsList);
    nCats = size(resultsList{1}.(fieldName), 2);
    nAligned = length(alignedTimeVec);

    data3D = nan(nAligned, nCats, nExp);

    for e = 1:nExp

        name = char(names(e));

        if isKey(folderTimeShiftMinutes, name)
            shiftMin = folderTimeShiftMinutes(name);
        else
            shiftMin = 0;
        end

        nativeTime = alignedTimeVec + shiftMin;
        nativeFrame = round(nativeTime / minutesPerFrame) + 1;

        thisData = resultsList{e}.(fieldName);
        validRows = nativeFrame >= 1 & nativeFrame <= size(thisData, 1);

        data3D(validRows, :, e) = thisData(nativeFrame(validRows), :);
    end

end

function [meanMat, semMat] = meanSEM(data3D)
% Mean and SEM across the 3rd (experiment) dimension, ignoring NaN so
% experiments missing data at a given aligned time point don't drag it down.

    meanMat = mean(data3D, 3, "omitnan");
    nValid = sum(~isnan(data3D), 3);
    semMat = std(data3D, 0, 3, "omitnan") ./ sqrt(nValid);
    semMat(nValid < 2) = NaN;   % can't estimate spread from < 2 experiments

end

function plotAndExportStackedBar(data3D, catLabels, names, alignedTimeVec, ttl, yLab, outFolder, pngName, excelFile, sheetName)

    [meanMat, semMat] = meanSEM(data3D);

    figure("Name", ttl)
    bar(alignedTimeVec, meanMat, "stacked")
    hold on

    cumBase = [zeros(size(meanMat, 1), 1), cumsum(meanMat(:, 1:end-1), 2)];

    for c = 1:size(meanMat, 2)
        errorbar(alignedTimeVec, cumBase(:, c) + meanMat(:, c), semMat(:, c), ...
            "LineStyle", "none", "Marker", ".", "Color", "k", "CapSize", 3)
    end

    xlabel("Aligned time (min)")
    ylabel(yLab)
    title(ttl)
    legend(catLabels, "Location", "eastoutside")
    saveas(gcf, fullfile(outFolder, pngName + ".png"))

    writeReplicateMeanSEMSheet(excelFile, sheetName, alignedTimeVec, catLabels, names, data3D, meanMat, semMat)

end

function plotAndExportCurve(data3D, catLabels, names, alignedTimeVec, ttl, yLab, outFolder, pngNameBase, excelFile, sheetName)

    [meanMat, semMat] = meanSEM(data3D);

    figure("Name", ttl + " (linear)")
    hold on
    for c = 1:size(meanMat, 2)
        errorbar(alignedTimeVec, meanMat(:, c), semMat(:, c), "LineWidth", 1.5)
    end
    xlabel("Aligned time (min)")
    ylabel(yLab)
    title(ttl + " (linear)")
    legend(catLabels, "Location", "best")
    saveas(gcf, fullfile(outFolder, pngNameBase + "_Linear.png"))

    figure("Name", ttl + " (log)")
    hold on
    for c = 1:size(meanMat, 2)
        errorbar(alignedTimeVec, meanMat(:, c), semMat(:, c), "LineWidth", 1.5)
    end
    set(gca, "YScale", "log")
    xlabel("Aligned time (min)")
    ylabel(yLab + " (log scale)")
    title(ttl + " (log)")
    legend(catLabels, "Location", "best")
    saveas(gcf, fullfile(outFolder, pngNameBase + "_Log.png"))

    writeReplicateMeanSEMSheet(excelFile, sheetName, alignedTimeVec, catLabels, names, data3D, meanMat, semMat)

end

function writeReplicateMeanSEMSheet(excelFile, sheetName, alignedTimeVec, catLabels, names, data3D, meanMat, semMat)
% Writes ONE sheet per metric. For each category, columns are: every
% replicate's own value (labeled by folder name), then Mean, then SEM --
% e.g. "Microcolony_20260226", "Microcolony_20260227", "Microcolony_Mean",
% "Microcolony_SEM".

    nCats = length(catLabels);
    nExp = length(names);

    header = "Time_min";
    dataCols = {alignedTimeVec};

    for c = 1:nCats
        for e = 1:nExp
            header(end+1) = catLabels(c) + "_" + names(e); %#ok<AGROW>
            dataCols{end+1} = data3D(:, c, e); %#ok<AGROW>
        end
        header(end+1) = catLabels(c) + "_Mean"; %#ok<AGROW>
        dataCols{end+1} = meanMat(:, c); %#ok<AGROW>
        header(end+1) = catLabels(c) + "_SEM"; %#ok<AGROW>
        dataCols{end+1} = semMat(:, c); %#ok<AGROW>
    end

    nRows = length(alignedTimeVec);
    nCols = length(dataCols);

    outCell = cell(nRows + 1, nCols);
    outCell(1, :) = cellstr(header);

    for col = 1:nCols
        outCell(2:end, col) = num2cell(dataCols{col});
    end

    writecell(outCell, excelFile, "Sheet", sheetName)

end
