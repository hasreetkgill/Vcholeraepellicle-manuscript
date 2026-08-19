%% ClassifyObjectsDistanceDensity_Coculture.m
% Author: Hasreet Gill, Bassler Lab, HHMI/Princeton University
% Revised with Claude (Sonnet 5, Anthropic)

% For red/green COCULTURE imaging: processes every matched 488(green)/
% 561(red) "*_Denoised*.tif" Z-stack pair in a folder (pairs are
% mandatory here -- a 488 file with no matching 561 file is skipped
% entirely). For each channel independently:
%   1. Thresholds the volume (image-adaptive percentile, with several
%      hand-tuned per-experiment alternatives left as comments below --
%      edit the active branch to match whatever dataset this is run on).
%   2. Picks, per 3D connected object, its single best Z-slice (largest
%      area) to build one flat "selected slice" image -- discarding
%      objects that span fewer than 3 slices or whose max-projection
%      footprint is much larger than their best single slice (both
%      signs of a transient/motion artifact rather than a real static
%      object, same idea as the "swimmer" exclusion elsewhere in this
%      pipeline family).
%   3. Classifies the resulting 2D objects as SingleCell vs.
%      GreaterThanSingleCell (by area/aspect-ratio/watershed-split
%      status), and separately flags which ones are "Large" (large
%      objects are inherently a subset of GreaterThanSingleCell, since
%      largeAreaThresh > singleMaxArea).
% Then, for the PAIR together: builds a green+red composite/boundary
% overlay, and computes -- separately using EACH color's own large
% objects as the distance-reference surface -- how the mix of nearby
% green vs. red small objects changes with distance from that surface,
% binned and density-normalized, then compiled to mean +/- SEM plots
% across replicates.

clear
clc
close all

%% PARAMETERS

binWidth = 45.5;                % distance bin width in pixels
densityScale = 45455;           % report density as objects per this many available pixels

minObjectArea = 50;             % remove 3D objects if largest single-slice area is below this
% projRatioThresh = 5;          % WT delmshA % put object in less-than-3 group if max projection area is too spread out
projRatioThresh = 1.5;          % Protein addition % put object in less-than-3 group if max projection area is too spread out

singleMinArea = 50;             % exclude objects below this from object output
singleMaxArea = 500;            % single-cell upper area cutoff
singleMaxAspectRatio = 4;       % single-cell aspect ratio cutoff
largeAreaThresh = 1500;         % large objects must be > this area

objectSetNames = ["AllObjects"; ...
                  "SingleCellObjects"; ...
                  "GreaterThanSingleCellObjects"];

%% OUTPUT CONTROL FLAGS
% Turn these on when tuning thresholds on new data and you need to
% visually inspect intermediate steps. Leave off for routine runs -
% the compiled tables/plots below are the outputs actually used for
% quantification, and everything gated by these flags is either
% redundant with those compiled outputs or derivable from
% *_IncludedObjectLabels.tif plus the classified object table.

saveDebugMasks = false;         % per-channel intermediate tif outputs: masked
                                 % <3-slice/>=3-slice stacks, selected-slice
                                 % binary/gray, per-class binary masks,
                                 % boundary mask, classification overlay

savePerFileTables = false;      % per-channel and per-pair CSVs that duplicate
                                 % rows already present in the compiled CSVs

savePerReplicatePlots = false;  % per-pair distance-density bar plots (one per
                                 % object set per replicate); the compiled
                                 % per-condition mean+SEM plots are kept regardless

saveFigFiles = false;           % also save MATLAB .fig alongside .png for every plot

%% CHOOSE FOLDER AND FIND INPUT FILES

path = uigetdir();
cd(path)

allTifs = dir("*Denoised*.tif");
allNames = string({allTifs.name});

% Excludes every kind of output this script (or a re-run of it) writes,
% so re-running on the same folder doesn't try to re-ingest its own
% results as if they were fresh input images.
isInput = ~contains(allNames, "_Masked_") & ...
          ~contains(allNames, "_SelectedSlice_") & ...
          ~contains(allNames, "_SmallObjects_") & ...
          ~contains(allNames, "_Objects_") & ...
          ~contains(allNames, "_Compiled_") & ...
          ~contains(allNames, "_Intermixing_") & ...
          ~contains(allNames, "_SingleCells_") & ...
          ~contains(allNames, "_GreaterThanSingleCell_") & ...
          ~contains(allNames, "_LargeObjects_") & ...
          ~contains(allNames, "_IncludedObject_") & ...
          ~contains(allNames, "_DistanceToLarge");

inputNames = allNames(isInput);

greenFiles = inputNames(contains(inputNames, "488"));

hasRedPair = false(size(greenFiles));

for i = 1:length(greenFiles)

    testPair = replace(greenFiles(i), "488", "561");

    if any(inputNames == testPair)
        hasRedPair(i) = true;
    end

end

% Only 488 files with a genuine 561 partner drive the main loop below --
% unpaired files are silently dropped (a 561-only file never gets its
% own turn since the loop is keyed on the 488 side).
startFiles = greenFiles(hasRedPair);

if isempty(startFiles)
    error("No paired 488/561 Denoised files were found.")
end

compiledObjectTable = table();
compiledSummaryTable = table();
compiledDistanceBinTableAll = table();

%% PROCESS EACH 488/561 PAIR

for s = 1:length(startFiles)

    fname488 = startFiles(s);
    fname561 = replace(fname488, "488", "561");

    files = [fname488 fname561];
    channelNames = ["488" "561"];
    colorNames = ["Green" "Red"];

    tok = regexp(char(fname488), '^(.*)_(\d+)_Denoised', 'tokens', 'once');

    if ~isempty(tok)
        conditionName = string(tok{1});
        replicateNum = str2double(tok{2});
    else
        conditionName = erase(fname488, ".tif");
        replicateNum = s;
    end

    nFiles = length(files);

    selectedBinCell = cell(nFiles,1);
    selectedGrayCell = cell(nFiles,1);
    includedMaskCell = cell(nFiles,1);
    singleMaskCell = cell(nFiles,1);
    greaterMaskCell = cell(nFiles,1);
    largeMaskCell = cell(nFiles,1);
    boundaryMaskCell = cell(nFiles,1);
    labelImageCell = cell(nFiles,1);
    rawTableCell = cell(nFiles,1);

%% PROCESS EACH CHANNEL

    for f = 1:nFiles

        im = tiffreadVolume(files(f));
        imDouble = double(im);

        [H,W,Z] = size(im);

        % thresh = prctile(imDouble(:), 99); % Protein addition = 98, WTdelmshA = 99

        % Active threshold rule for this dataset: strain code "340"
        % (a spike-in WT code, per the SpikeIn labeling convention used
        % elsewhere in this pipeline) gets its own looser cutoff;
        % otherwise green (f==1) and red (f==2) get different percentile
        % cutoffs, tuned for this "spike-in QS mutants" condition.
        if contains(fname488, "340")
            thresh = prctile(imDouble(:), 97);
        else
            if f == 1
                thresh = prctile(imDouble(:), 99.5); % Spike in QS mutants = 99.8
            else
                thresh = prctile(imDouble(:), 98.5); % Spike in QS mutants = 98.5
            end
        end

        % WT delmshA alternative
        % if contains(fname1, "del")
        %     thresh = prctile(imDouble(:), 99.5);
        % else
        %     thresh = prctile(imDouble(:), 99);
        % end

        bin = imDouble >= thresh;

        CC = bwconncomp(bin, 26);

        binLessThan3 = false(size(bin));
        bin3OrMore = false(size(bin));

        selectedBin = false(H,W);

        bgVals = imDouble(~bin);

        if isempty(bgVals)
            bg = median(imDouble(:));
        else
            bg = median(bgVals);
        end

        selectedGray = zeros(H,W, "like", im);
        selectedGray(:) = bg;

        % For every 3D thresholded object: keep it only if its largest
        % single-slice area clears minObjectArea, then decide whether
        % it's a real static object (>=3 slices AND not badly spread out
        % in max-projection relative to its best slice) or a likely
        % transient/motion artifact (binLessThan3, excluded from
        % everything downstream). Real objects contribute only their
        % single best (largest-area) slice into the flat 2D
        % selectedBin/selectedGray images used for all classification
        % below -- not their full 3D extent.
        for j = 1:CC.NumObjects

            pix = CC.PixelIdxList{j};
            [y,x,z] = ind2sub(size(bin), pix);

            zList = unique(z);
            areas = zeros(length(zList),1);

            for k = 1:length(zList)
                areas(k) = sum(z == zList(k));
            end

            largestSingleSliceArea = max(areas);

            xyPix = unique(sub2ind([H W], y, x));
            maxProjArea = length(xyPix);

            if largestSingleSliceArea >= minObjectArea

                if length(zList) < 3 || maxProjArea > projRatioThresh * largestSingleSliceArea

                    binLessThan3(pix) = true;

                else

                    bin3OrMore(pix) = true;

                    [~,idxMax] = max(areas);
                    bestZ = zList(idxMax);

                    keep = z == bestZ;
                    idx2D = sub2ind([H W], y(keep), x(keep));

                    selectedBin(idx2D) = true;

                    oneSlice = im(:,:,bestZ);
                    selectedGray(idx2D) = oneSlice(idx2D);

                end

            end

        end

        selectedBinCell{f} = selectedBin;
        selectedGrayCell{f} = selectedGray;

        [~,baseName,~] = fileparts(files(f));

        if saveDebugMasks

            maskedLessThan3 = im;
            maskedLessThan3(~binLessThan3) = 0;

            masked3OrMore = im;
            masked3OrMore(~bin3OrMore) = 0;

            outLess = baseName + "_Masked_LessThan3Slices.tif";
            outMore = baseName + "_Masked_3OrMoreSlices.tif";
            outBin = baseName + "_SelectedSlice_Binary.tif";
            outGray = baseName + "_SelectedSlice_Gray.tif";

            if exist(outLess, "file")
                delete(outLess)
            end

            if exist(outMore, "file")
                delete(outMore)
            end

            for z = 1:Z

                if z == 1
                    imwrite(maskedLessThan3(:,:,z), outLess, "tif")
                    imwrite(masked3OrMore(:,:,z), outMore, "tif")
                else
                    imwrite(maskedLessThan3(:,:,z), outLess, "tif", "WriteMode", "append")
                    imwrite(masked3OrMore(:,:,z), outMore, "tif", "WriteMode", "append")
                end

            end

            imwrite(uint8(selectedBin)*255, outBin, "tif")
            imwrite(selectedGray, outGray, "tif")

        end

%% WATERSHED SPLIT TEST FOR THIS CHANNEL
% Same distance-transform + h-maxima + watershed recipe used elsewhere in
% this pipeline family: ridge marks exactly which pixels are the
% dividing lines between two objects that were touching before this
% split -- used below to flag which final classified objects had to be
% separated from a neighbor (WasSplit), not to actually change selectedBin.

        D = bwdist(~selectedBin);
        D = imhmax(D, 1);
        Lw = watershed(-D);
        ridge = Lw == 0;
        ridge = ridge & selectedBin;

%% CLASSIFY OBJECTS FOR THIS CHANNEL

        CCobjects = bwconncomp(selectedBin, 8);

        propsObjects = regionprops(CCobjects, "Area", "Centroid", "PixelIdxList", ...
            "MajorAxisLength", "MinorAxisLength", "BoundingBox");

        objectID = [];
        includedObjectID = [];
        objectArea = [];
        centroidX = [];
        centroidY = [];
        aspectRatio = [];
        wasSplit = [];
        watershedResult = strings(0,1);
        objectClass = strings(0,1);
        isLargeObject = [];
        boundingBoxX = [];
        boundingBoxY = [];
        boundingBoxWidth = [];
        boundingBoxHeight = [];

        includedMask = false(size(selectedBin));
        singleMask = false(size(selectedBin));
        greaterMask = false(size(selectedBin));
        largeMask = false(size(selectedBin));
        includedLabelImage = zeros(size(selectedBin), "uint16");

        includedCounter = 0;

        for j = 1:length(propsObjects)

            includeThisObject = propsObjects(j).Area >= singleMinArea;

            if includeThisObject

                includedCounter = includedCounter + 1;

                c = propsObjects(j).Centroid;

                if propsObjects(j).MinorAxisLength == 0
                    ar = inf;
                else
                    ar = propsObjects(j).MajorAxisLength / propsObjects(j).MinorAxisLength;
                end

                splitNow = any(ridge(propsObjects(j).PixelIdxList));

                isSingleCell = propsObjects(j).Area > singleMinArea && ...
                               propsObjects(j).Area < singleMaxArea && ...
                               ar < singleMaxAspectRatio && ...
                               splitNow == false;

                largeNow = propsObjects(j).Area > largeAreaThresh;

                % Every object is exactly one of SingleCell or
                % GreaterThanSingleCell; IsLargeObject is a further,
                % independent size flag layered on top (large objects
                % are always GreaterThanSingleCell too, since
                % largeAreaThresh > singleMaxArea).
                if isSingleCell
                    classNow = "SingleCell";
                    singleMask(propsObjects(j).PixelIdxList) = true;
                else
                    classNow = "GreaterThanSingleCell";
                    greaterMask(propsObjects(j).PixelIdxList) = true;
                end

                if largeNow
                    largeMask(propsObjects(j).PixelIdxList) = true;
                end

                if splitNow
                    watershedNow = "Split";
                else
                    watershedNow = "NotSplit";
                end

                includedMask(propsObjects(j).PixelIdxList) = true;
                includedLabelImage(propsObjects(j).PixelIdxList) = uint16(includedCounter);

                bb = propsObjects(j).BoundingBox;

                objectID = [objectID; j];
                includedObjectID = [includedObjectID; includedCounter];
                objectArea = [objectArea; propsObjects(j).Area];
                centroidX = [centroidX; c(1)];
                centroidY = [centroidY; c(2)];
                aspectRatio = [aspectRatio; ar];
                wasSplit = [wasSplit; splitNow];
                watershedResult = [watershedResult; watershedNow];
                objectClass = [objectClass; classNow];
                isLargeObject = [isLargeObject; largeNow];
                boundingBoxX = [boundingBoxX; bb(1)];
                boundingBoxY = [boundingBoxY; bb(2)];
                boundingBoxWidth = [boundingBoxWidth; bb(3)];
                boundingBoxHeight = [boundingBoxHeight; bb(4)];

            end

        end

        wasSplit = logical(wasSplit);
        isLargeObject = logical(isLargeObject);

        includedBoundary = boundarymask(includedMask);
        singleBoundary = boundarymask(singleMask);
        greaterBoundary = boundarymask(greaterMask);
        largeBoundary = boundarymask(largeMask);

        includedMaskCell{f} = includedMask;
        singleMaskCell{f} = singleMask;
        greaterMaskCell{f} = greaterMask;
        largeMaskCell{f} = largeMask;
        boundaryMaskCell{f} = includedBoundary;
        labelImageCell{f} = includedLabelImage;

%% SAVE OBJECT TABLE FOR THIS CHANNEL

        channelCol = repmat(channelNames(f), length(objectID), 1);
        colorCol = repmat(colorNames(f), length(objectID), 1);

        rawTable = table(objectID, includedObjectID, channelCol, colorCol, ...
            objectArea, centroidX, centroidY, aspectRatio, ...
            wasSplit, watershedResult, objectClass, isLargeObject, ...
            boundingBoxX, boundingBoxY, boundingBoxWidth, boundingBoxHeight, ...
            'VariableNames', {'ObjectID','IncludedObjectID','Channel','Color', ...
            'ObjectArea_pixels','CentroidX_pixels','CentroidY_pixels', ...
            'AspectRatio','WasSplit','WatershedResult','ObjectClass', ...
            'IsLargeObject','BoundingBoxX_pixels','BoundingBoxY_pixels', ...
            'BoundingBoxWidth_pixels','BoundingBoxHeight_pixels'});

        rawTable.Condition = repmat(conditionName, height(rawTable), 1);
        rawTable.Replicate = repmat(replicateNum, height(rawTable), 1);
        rawTable.SourceFile = repmat(files(f), height(rawTable), 1);
        rawTable.ThresholdUsed = repmat(thresh, height(rawTable), 1);
        rawTable = movevars(rawTable, ...
            {'Condition','Replicate','SourceFile','ThresholdUsed'}, ...
            'Before', 1);

        if savePerFileTables
            writetable(rawTable, baseName + "_Objects_Classified.csv")
        end

        rawTableCell{f} = rawTable;

        compiledObjectTable = [compiledObjectTable; rawTable];

%% SAVE OBJECT MASKS FOR THIS CHANNEL

        imwrite(uint16(includedLabelImage), baseName + "_IncludedObjectLabels.tif", "tif")

        if saveDebugMasks
            imwrite(uint8(includedMask)*255, baseName + "_IncludedObjects_Binary.tif", "tif")
            imwrite(uint8(singleMask)*255, baseName + "_SingleCells_Binary.tif", "tif")
            imwrite(uint8(greaterMask)*255, baseName + "_GreaterThanSingleCell_Binary.tif", "tif")
            imwrite(uint8(largeMask)*255, baseName + "_LargeObjects_Binary.tif", "tif")
            imwrite(uint8(includedBoundary)*255, baseName + "_IncludedObjectBoundaries_Binary.tif", "tif")
        end

%% SAVE CLASSIFICATION OVERLAY FOR THIS CHANNEL
% Green outline = SingleCell, blue outline = GreaterThanSingleCell, red
% outline = Large (drawn last, so a large object's outline always shows
% red even though it's also a GreaterThanSingleCell object underneath).

        if saveDebugMasks

            baseGray = mat2gray(selectedGray);
            rgb = repmat(baseGray, [1 1 3]);

            red = rgb(:,:,1);
            green = rgb(:,:,2);
            blue = rgb(:,:,3);

            red(singleBoundary) = 0;
            green(singleBoundary) = 1;
            blue(singleBoundary) = 0;

            red(greaterBoundary) = 0;
            green(greaterBoundary) = 0;
            blue(greaterBoundary) = 1;

            red(largeBoundary) = 1;
            green(largeBoundary) = 0;
            blue(largeBoundary) = 0;

            rgb(:,:,1) = red;
            rgb(:,:,2) = green;
            rgb(:,:,3) = blue;

            imwrite(im2uint8(rgb), baseName + "_SelectedSlice_ObjectClassOverlay.tif", "tif")

        end

%% SAVE SUMMARY ROW FOR THIS CHANNEL

        totalIncludedObjects = length(objectID);
        totalSingleCells = sum(objectClass == "SingleCell");
        totalGreaterThanSingleCell = sum(objectClass == "GreaterThanSingleCell");
        totalLargeObjects = sum(isLargeObject);

        if totalIncludedObjects > 0
            proportionSingleCells = totalSingleCells / totalIncludedObjects;
            proportionGreaterThanSingleCell = totalGreaterThanSingleCell / totalIncludedObjects;
            proportionLargeObjects = totalLargeObjects / totalIncludedObjects;
        else
            proportionSingleCells = NaN;
            proportionGreaterThanSingleCell = NaN;
            proportionLargeObjects = NaN;
        end

        summaryTable = table(conditionName, replicateNum, files(f), ...
            channelNames(f), colorNames(f), thresh, ...
            totalIncludedObjects, totalSingleCells, ...
            totalGreaterThanSingleCell, totalLargeObjects, ...
            proportionSingleCells, proportionGreaterThanSingleCell, ...
            proportionLargeObjects, ...
            'VariableNames', {'Condition','Replicate','SourceFile', ...
            'Channel','Color','ThresholdUsed', ...
            'TotalIncludedObjects','TotalSingleCells', ...
            'TotalGreaterThanSingleCellObjects','TotalLargeObjects', ...
            'ProportionSingleCells', ...
            'ProportionGreaterThanSingleCellObjects', ...
            'ProportionLargeObjects'});

        compiledSummaryTable = [compiledSummaryTable; summaryTable];

    end

%% SAVE PAIRED RED/GREEN SELECTED-SLICE COMPOSITE AND BOUNDARY OVERLAY
% rgbPair: red = 561 selected-slice, green = 488 selected-slice, no blue
% -- a direct two-color merge of the two channels' chosen slices.
% rgbBoundary overlays each channel's own object boundaries in its own
% color, with any pixel that is a boundary in BOTH channels drawn yellow
% (red+green) so overlapping outlines are visually distinguishable from
% either channel's outline alone.

    pairBase = erase(fname488, ".tif");
    pairBase = replace(pairBase, "488", "488_561");

    greenGray = mat2gray(selectedGrayCell{1});
    redGray = mat2gray(selectedGrayCell{2});

    rgbPair = zeros(size(greenGray,1), size(greenGray,2), 3);

    rgbPair(:,:,1) = redGray;
    rgbPair(:,:,2) = greenGray;
    rgbPair(:,:,3) = 0;

    imwrite(im2uint8(rgbPair), pairBase + "_Paired_SelectedSlice_Composite.tif", "tif")

    rgbBoundary = rgbPair;

    redBoundary = boundaryMaskCell{2};
    greenBoundary = boundaryMaskCell{1};

    redPlane = rgbBoundary(:,:,1);
    greenPlane = rgbBoundary(:,:,2);
    bluePlane = rgbBoundary(:,:,3);

    redPlane(redBoundary) = 1;
    greenPlane(redBoundary) = 0;
    bluePlane(redBoundary) = 0;

    redPlane(greenBoundary) = 0;
    greenPlane(greenBoundary) = 1;
    bluePlane(greenBoundary) = 0;

    bothBoundary = redBoundary & greenBoundary;

    redPlane(bothBoundary) = 1;
    greenPlane(bothBoundary) = 1;
    bluePlane(bothBoundary) = 0;

    rgbBoundary(:,:,1) = redPlane;
    rgbBoundary(:,:,2) = greenPlane;
    rgbBoundary(:,:,3) = bluePlane;

    imwrite(im2uint8(rgbBoundary), pairBase + "_Paired_SelectedSlice_BoundaryOverlay.tif", "tif")

%% DISTANCE-BINNED GREEN/RED OBJECT DENSITY AROUND LARGE OBJECTS (SPLIT BY LARGE-OBJECT COLOR)
% Core coculture-specific analysis: run TWICE per pair, once treating
% green's own large objects as the distance-reference surface, once
% treating red's -- for each, bins every non-large small object (all
% objects / SingleCell only / GreaterThanSingleCell only) by distance to
% that surface, and counts/densities/proportions how many are green vs.
% red in each bin. Large objects of EITHER color are always excluded
% from the counted population and from the available-pixel pool, no
% matter which color is currently the reference surface, so the two
% passes are directly comparable.

    pairObjectTable = [rawTableCell{1}; rawTableCell{2}];

    if height(pairObjectTable) > 0
        pairObjectTable.PairBase = repmat(pairBase, height(pairObjectTable), 1);
        pairObjectTable = movevars(pairObjectTable, "PairBase", "After", "Replicate");
        pairObjectTable.IsLargeObject = logical(pairObjectTable.IsLargeObject);
        pairObjectTable.WasSplit = logical(pairObjectTable.WasSplit);
    end

    largeColorNames = ["Green" "Red"];

    % Large objects of either color are excluded from every counted
    % population and from the available-pixel pool, regardless of which
    % color's large objects are currently the distance reference surface.
    excludeLargeMask = largeMaskCell{1} | largeMaskCell{2};

    pairDistanceBinTable = table();

    for largeColorNum = 1:length(largeColorNames)

        largeColorNow = largeColorNames(largeColorNum);
        largeMaskNow = largeMaskCell{largeColorNum};

        if any(largeMaskNow(:))
            distToLarge = bwdist(largeMaskNow);
            availableDistVals = distToLarge(~excludeLargeMask);
            availableDistVals(isinf(availableDistVals)) = nan;
            maxDist = max(availableDistVals, [], "omitnan");
        else
            distToLarge = nan(size(largeMaskNow));
            maxDist = 0;
        end

        if isempty(maxDist) || isnan(maxDist)
            maxDist = 0;
        end

        edges = 0:binWidth:(ceil(maxDist/binWidth)*binWidth + binWidth);

        if length(edges) < 2
            edges = [0 binWidth];
        end

        binStart = edges(1:end-1)';
        binEnd = edges(2:end)';
        binLabel = strings(length(binStart),1);

        for j = 1:length(binStart)
            binLabel(j) = string(binStart(j)) + "-" + string(binEnd(j));
        end

        distanceToLarge_pixels = nan(height(pairObjectTable),1);

        for j = 1:height(pairObjectTable)

            cx = pairObjectTable.CentroidX_pixels(j);
            cy = pairObjectTable.CentroidY_pixels(j);

            if any(largeMaskNow(:))
                distanceToLarge_pixels(j) = interp2(distToLarge, cx, cy, "linear");
            else
                distanceToLarge_pixels(j) = NaN;
            end

        end

        pairObjectTableNow = pairObjectTable;
        pairObjectTableNow.DistanceToNearestLargeObject_pixels = distanceToLarge_pixels;

        for setNum = 1:length(objectSetNames)

            objectSetNow = objectSetNames(setNum);

            isLargeObjectNow = logical(pairObjectTableNow.IsLargeObject);

            if objectSetNow == "AllObjects"
                objectRows = ~isLargeObjectNow;
            elseif objectSetNow == "SingleCellObjects"
                objectRows = ~isLargeObjectNow & pairObjectTableNow.ObjectClass == "SingleCell";
            elseif objectSetNow == "GreaterThanSingleCellObjects"
                objectRows = ~isLargeObjectNow & pairObjectTableNow.ObjectClass == "GreaterThanSingleCell";
            else
                objectRows = false(height(pairObjectTableNow),1);
            end

            numGreenObjects = zeros(length(binStart),1);
            numRedObjects = zeros(length(binStart),1);
            totalObjects = zeros(length(binStart),1);

            pixelsInBin = zeros(length(binStart),1);

            densityGreenObjects = nan(length(binStart),1);
            densityRedObjects = nan(length(binStart),1);
            densityTotalObjects = nan(length(binStart),1);

            proportionGreen = nan(length(binStart),1);
            proportionRed = nan(length(binStart),1);

            for bnum = 1:length(binStart)

                inBinObjects = objectRows & ...
                               pairObjectTableNow.DistanceToNearestLargeObject_pixels >= binStart(bnum) & ...
                               pairObjectTableNow.DistanceToNearestLargeObject_pixels < binEnd(bnum);

                numGreenObjects(bnum) = sum(inBinObjects & pairObjectTableNow.Color == "Green");
                numRedObjects(bnum) = sum(inBinObjects & pairObjectTableNow.Color == "Red");
                totalObjects(bnum) = numGreenObjects(bnum) + numRedObjects(bnum);

                if any(largeMaskNow(:))
                    availablePixels = distToLarge >= binStart(bnum) & ...
                                      distToLarge < binEnd(bnum) & ...
                                      ~excludeLargeMask;

                    pixelsInBin(bnum) = sum(availablePixels(:));
                else
                    pixelsInBin(bnum) = 0;
                end

                if pixelsInBin(bnum) > 0
                    densityGreenObjects(bnum) = numGreenObjects(bnum) / pixelsInBin(bnum) * densityScale;
                    densityRedObjects(bnum) = numRedObjects(bnum) / pixelsInBin(bnum) * densityScale;
                    densityTotalObjects(bnum) = totalObjects(bnum) / pixelsInBin(bnum) * densityScale;
                else
                    densityGreenObjects(bnum) = NaN;
                    densityRedObjects(bnum) = NaN;
                    densityTotalObjects(bnum) = NaN;
                end

                if totalObjects(bnum) > 0
                    proportionGreen(bnum) = numGreenObjects(bnum) / totalObjects(bnum);
                    proportionRed(bnum) = numRedObjects(bnum) / totalObjects(bnum);
                else
                    proportionGreen(bnum) = NaN;
                    proportionRed(bnum) = NaN;
                end

            end

            oneSetTable = table(repmat(conditionName, length(binStart), 1), ...
                repmat(replicateNum, length(binStart), 1), ...
                repmat(pairBase, length(binStart), 1), ...
                repmat(largeColorNow, length(binStart), 1), ...
                repmat(objectSetNow, length(binStart), 1), ...
                binLabel, binStart, binEnd, ...
                repmat(densityScale, length(binStart), 1), ...
                pixelsInBin, ...
                numGreenObjects, numRedObjects, totalObjects, ...
                densityGreenObjects, densityRedObjects, densityTotalObjects, ...
                proportionGreen, proportionRed, ...
                'VariableNames', {'Condition','Replicate','PairBase','LargeObjectColor','ObjectSet', ...
                'DistanceBin','BinStart_pixels','BinEnd_pixels', ...
                'DensityScale_pixels','PixelsInBin', ...
                'NumGreenObjects','NumRedObjects','TotalObjects', ...
                'DensityGreenObjects_perDensityScalePixels', ...
                'DensityRedObjects_perDensityScalePixels', ...
                'DensityTotalObjects_perDensityScalePixels', ...
                'ProportionGreen','ProportionRed'});

            pairDistanceBinTable = [pairDistanceBinTable; oneSetTable];

%% SAVE PER-PAIR DISTANCE-BIN PLOT FOR THIS LARGE-OBJECT COLOR AND OBJECT SET

            if savePerReplicatePlots

                figure

                binCats = categorical(binLabel, binLabel, "Ordinal", true);

                Ydensity = [densityGreenObjects densityRedObjects];

                b = bar(binCats, Ydensity, "stacked");

                b(1).FaceColor = [0 0.7 0];
                b(2).FaceColor = [1 0 0];

                xlabel("Distance to nearest " + largeColorNow + " large object edge, pixels")
                ylabel("Objects per density-scale pixels")
                title(pairBase + " " + objectSetNow + " by distance to " + largeColorNow + " large objects")
                legend(["Green objects", "Red objects"], "Location", "best")

                xtickangle(45)

                outObjectSet = regexprep(objectSetNow, '[^\w\-]', '_');

                saveas(gcf, pairBase + "_DistanceToLarge" + largeColorNow + "_" + outObjectSet + "_StackedDensity_GreenRed.png")

                if saveFigFiles
                    saveas(gcf, pairBase + "_DistanceToLarge" + largeColorNow + "_" + outObjectSet + "_StackedDensity_GreenRed.fig")
                end

                close(gcf)

            end

        end

    end

    if savePerFileTables
        writetable(pairDistanceBinTable, pairBase + "_DistanceToLarge_BinnedGreenRedDensity.csv")
    end

    compiledDistanceBinTableAll = [compiledDistanceBinTableAll; pairDistanceBinTable];

end

%% SAVE COMPILED TABLES

if height(compiledObjectTable) > 0
    writetable(compiledObjectTable, "Compiled_TwoChannel_Objects_Classified.csv")
end

if height(compiledSummaryTable) > 0
    writetable(compiledSummaryTable, "Compiled_TwoChannel_ObjectSummary.csv")
end

if height(compiledDistanceBinTableAll) > 0
    writetable(compiledDistanceBinTableAll, "Compiled_DistanceToLarge_BinnedGreenRedDensity_AllReplicates.csv")
end

%% COMPILE DISTANCE-BIN GREEN/RED DENSITY ACROSS REPLICATES
% For every (condition, large-object color, object set) combination,
% averages each distance bin's density/proportion values across every
% replicate that contributed data, then plots the mean +/- SEM stacked
% green/red density as one compiled figure per combination.

if height(compiledDistanceBinTableAll) > 0

    conditionsDensity = unique(compiledDistanceBinTableAll.Condition, "stable");
    largeColorNamesCompiled = unique(compiledDistanceBinTableAll.LargeObjectColor, "stable");

    compiledDensityStats = table();

    for c = 1:length(conditionsDensity)

        conditionNow = conditionsDensity(c);

        for largeColorNum = 1:length(largeColorNamesCompiled)

            largeColorNow = largeColorNamesCompiled(largeColorNum);

            for setNum = 1:length(objectSetNames)

                objectSetNow = objectSetNames(setNum);

                setRows = compiledDistanceBinTableAll.Condition == conditionNow & ...
                          compiledDistanceBinTableAll.LargeObjectColor == largeColorNow & ...
                          compiledDistanceBinTableAll.ObjectSet == objectSetNow;

                binsNow = unique(compiledDistanceBinTableAll.DistanceBin(setRows), "stable");

                for bnum = 1:length(binsNow)

                    binNow = binsNow(bnum);

                    rowsNow = setRows & compiledDistanceBinTableAll.DistanceBin == binNow;

                    valsGreen = compiledDistanceBinTableAll.DensityGreenObjects_perDensityScalePixels(rowsNow);
                    valsRed = compiledDistanceBinTableAll.DensityRedObjects_perDensityScalePixels(rowsNow);
                    valsTotal = compiledDistanceBinTableAll.DensityTotalObjects_perDensityScalePixels(rowsNow);

                    valsPropGreen = compiledDistanceBinTableAll.ProportionGreen(rowsNow);
                    valsPropRed = compiledDistanceBinTableAll.ProportionRed(rowsNow);

                    nReps = sum(~isnan(valsTotal));

                    if nReps > 0

                        meanGreenDensity = mean(valsGreen, "omitnan");
                        semGreenDensity = std(valsGreen, 0, "omitnan") / sqrt(sum(~isnan(valsGreen)));

                        meanRedDensity = mean(valsRed, "omitnan");
                        semRedDensity = std(valsRed, 0, "omitnan") / sqrt(sum(~isnan(valsRed)));

                        meanTotalDensity = mean(valsTotal, "omitnan");
                        semTotalDensity = std(valsTotal, 0, "omitnan") / sqrt(sum(~isnan(valsTotal)));

                        meanProportionGreen = mean(valsPropGreen, "omitnan");
                        semProportionGreen = std(valsPropGreen, 0, "omitnan") / sqrt(sum(~isnan(valsPropGreen)));

                        meanProportionRed = mean(valsPropRed, "omitnan");
                        semProportionRed = std(valsPropRed, 0, "omitnan") / sqrt(sum(~isnan(valsPropRed)));

                    else

                        meanGreenDensity = NaN;
                        semGreenDensity = NaN;

                        meanRedDensity = NaN;
                        semRedDensity = NaN;

                        meanTotalDensity = NaN;
                        semTotalDensity = NaN;

                        meanProportionGreen = NaN;
                        semProportionGreen = NaN;

                        meanProportionRed = NaN;
                        semProportionRed = NaN;

                    end

                    binStartNow = compiledDistanceBinTableAll.BinStart_pixels(find(rowsNow, 1));
                    binEndNow = compiledDistanceBinTableAll.BinEnd_pixels(find(rowsNow, 1));

                    oneStatsRow = table(conditionNow, largeColorNow, objectSetNow, binNow, ...
                        binStartNow, binEndNow, nReps, densityScale, ...
                        meanGreenDensity, semGreenDensity, ...
                        meanRedDensity, semRedDensity, ...
                        meanTotalDensity, semTotalDensity, ...
                        meanProportionGreen, semProportionGreen, ...
                        meanProportionRed, semProportionRed, ...
                        'VariableNames', {'Condition','LargeObjectColor','ObjectSet','DistanceBin', ...
                        'BinStart_pixels','BinEnd_pixels','NReplicates', ...
                        'DensityScale_pixels', ...
                        'MeanDensityGreenObjects_perDensityScalePixels', ...
                        'SEMDensityGreenObjects_perDensityScalePixels', ...
                        'MeanDensityRedObjects_perDensityScalePixels', ...
                        'SEMDensityRedObjects_perDensityScalePixels', ...
                        'MeanDensityTotalObjects_perDensityScalePixels', ...
                        'SEMDensityTotalObjects_perDensityScalePixels', ...
                        'MeanProportionGreen','SEMProportionGreen', ...
                        'MeanProportionRed','SEMProportionRed'});

                    compiledDensityStats = [compiledDensityStats; oneStatsRow];

                end

%% SAVE COMPILED PLOT FOR THIS CONDITION, LARGE-OBJECT COLOR, AND OBJECT SET

                rowsPlot = compiledDensityStats.Condition == conditionNow & ...
                           compiledDensityStats.LargeObjectColor == largeColorNow & ...
                           compiledDensityStats.ObjectSet == objectSetNow;

                if any(rowsPlot)

                    plotTable = compiledDensityStats(rowsPlot,:);

                    [~,ord] = sort(plotTable.BinStart_pixels);
                    plotTable = plotTable(ord,:);

                    figure

                    binCats = categorical(plotTable.DistanceBin, plotTable.DistanceBin, "Ordinal", true);

                    Y = [plotTable.MeanDensityGreenObjects_perDensityScalePixels, ...
                         plotTable.MeanDensityRedObjects_perDensityScalePixels];

                    b = bar(binCats, Y, "stacked");

                    b(1).FaceColor = [0 0.7 0];
                    b(2).FaceColor = [1 0 0];

                    hold on

                    x = 1:height(plotTable);

                    errorbar(x, plotTable.MeanDensityGreenObjects_perDensityScalePixels, ...
                        plotTable.SEMDensityGreenObjects_perDensityScalePixels, ...
                        "k", "LineStyle", "none", "LineWidth", 1)

                    errorbar(x, plotTable.MeanDensityTotalObjects_perDensityScalePixels, ...
                        plotTable.SEMDensityTotalObjects_perDensityScalePixels, ...
                        "k", "LineStyle", "none", "LineWidth", 1)

                    hold off

                    xlabel("Distance to nearest " + largeColorNow + " large object edge, pixels")
                    ylabel("Mean objects per density-scale pixels")
                    title(conditionNow + " " + objectSetNow + " by distance to " + largeColorNow + " large objects")
                    legend(["Green objects", "Red objects"], "Location", "best")

                    xtickangle(45)

                    outCondition = regexprep(conditionNow, '[^\w\-]', '_');
                    outObjectSet = regexprep(objectSetNow, '[^\w\-]', '_');

                    saveas(gcf, outCondition + "_Compiled_DistanceToLarge" + largeColorNow + "_" + outObjectSet + "_StackedDensity_GreenRed.png")

                    if saveFigFiles
                        saveas(gcf, outCondition + "_Compiled_DistanceToLarge" + largeColorNow + "_" + outObjectSet + "_StackedDensity_GreenRed.fig")
                    end

                    close(gcf)

                end

            end

        end

    end

    writetable(compiledDensityStats, "Compiled_DistanceToLarge_BinnedGreenRedDensity_MeanSEM.csv")

end
