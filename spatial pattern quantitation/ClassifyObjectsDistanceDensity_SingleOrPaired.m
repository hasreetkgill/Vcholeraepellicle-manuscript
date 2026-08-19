%% ClassifyObjectsDistanceDensity_SingleOrPaired.m
% Author: Hasreet Gill, Bassler Lab, HHMI/Princeton University
% Revised with Claude (Sonnet 5, Anthropic)

% For everything OTHER than symmetric red/green coculture imaging:
% single-channel "*_Denoised*.tif" volumes, or paired 488/561 volumes
% where the two channels play ASYMMETRIC roles rather than being treated
% as two competing populations. One channel (smallChannel, 488) always
% supplies the "small objects" being classified/measured; the other
% (largeChannel, 561) always supplies the "large objects" used purely as
% a distance-reference surface -- e.g. spike-in particles (small)
% relative to a host colony (large), or a monoculture image where the
% same channel serves as both (large blobs found and excluded/referenced
% within the very same image the small objects come from).
%
% Unlike the coculture script, a pair is optional here: if a 561 file
% has no matching 488 file (or vice versa), it's still processed alone,
% with that one image acting as both the small- and large-object source.
%
% Per input (single file or pair):
%   1. Thresholds each image (99th percentile by default; several
%      per-experiment-type alternatives are left as comments below --
%      edit the active branch/threshold before running on a different
%      experiment type than whatever this was last tuned for).
%   2. Picks each 3D object's single best Z-slice the same way as the
%      coculture script (drops <3-slice or over-spread objects as likely
%      motion artifacts).
%   3. Finds large objects (>largeAreaThresh) in the large-object image,
%      and classifies small objects (SingleCell vs. GreaterThanSingleCell,
%      by area/aspect-ratio/watershed-split) in the small-object image --
%      excluding any small-labeled object that sits entirely inside a
%      large object's footprint (not a standalone object), and, in the
%      single-channel case specifically, also excluding any object big
%      enough to itself qualify as "large" from being double-counted as
%      a small object too.
%   4. Bins small objects by distance to the nearest large object edge,
%      with both raw counts and area-normalized density, then compiles
%      counts/density/summary statistics (mean +/- SEM) across
%      replicates and conditions.

clear
clc
close all

%% PARAMETERS

binWidth = 45.5;                  % distance bin width in pixels
densityScale = 45455;           % report density as objects per this many available pixels

smallChannel = "488";           % in paired case, small objects come from this channel
largeChannel = "561";           % in paired case, large objects come from this channel

minObjectArea = 50;             % remove 3D objects if largest single-slice area is below this
projRatioThresh = 5;          % WT delmshA % put object in less-than-3 group if max projection area is too spread out
% projRatioThresh = 1.5;          % Protein addition % put object in less-than-3 group if max projection area is too spread out

singleMinArea = 50;             % exclude objects below this from small-object output
singleMaxArea = 500;            % single-cell upper area cutoff
singleMaxAspectRatio = 4;       % single-cell aspect ratio cutoff
largeAreaThresh = 1500;         % large objects must be > this area in the large-object image

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
          ~contains(allNames, "_Compiled_");

inputNames = allNames(isInput);

startFiles = strings(0);

% Unlike the coculture script's pairing (which drives the loop from the
% 488 side and requires a 561 partner to exist), this walks every input
% file and drops only the 561 half of any 488/561 pair that's actually
% present -- so the loop below is driven by either a lone file (single-
% channel or an unpaired 561) or the 488 half of a genuine pair, and each
% pair is only ever processed once.
for i = 1:length(inputNames)

    thisFile = inputNames(i);
    useThisFile = true;

    if contains(thisFile, "561")

        test488 = replace(thisFile, "561", "488");

        if any(inputNames == test488)
            useThisFile = false;
        end

    end

    if useThisFile
        startFiles = [startFiles thisFile];
    end

end

compiledObjectTable = table();
compiledSummaryTable = table();
compiledBinTableAll = table();

for s = 1:length(startFiles)

    fname1 = startFiles(s);
    fname2 = "";

    tok = regexp(char(fname1), '^(.*)_(\d+)_Denoised', 'tokens', 'once');

    if ~isempty(tok)
        conditionName = string(tok{1});
        replicateNum = str2double(tok{2});
    else
        conditionName = erase(fname1, ".tif");
        replicateNum = s;
    end

    % Looks for the OTHER half of a pair regardless of which channel
    % fname1 happens to be; if none exists, fname2 stays "" and this
    % replicate is processed as a single image below.
    if contains(fname1, "488")

        testPair = replace(fname1, "488", "561");

        if exist(testPair, "file")
            fname2 = testPair;
        end

    end

    if contains(fname1, "561")

        testPair = replace(fname1, "561", "488");

        if exist(testPair, "file")
            fname2 = testPair;
        end

    end

    files = fname1;

    if fname2 ~= ""
        files = [fname1 fname2];
    end

    nFiles = length(files);

    selectedBinCell = cell(nFiles,1);
    selectedGrayCell = cell(nFiles,1);

%% PROCESS EACH IMAGE
% Same threshold -> 3D-object triage -> single-best-slice selection
% pipeline as the coculture script, run independently on every file in
% this replicate (one or two).

    for f = 1:nFiles

        im = tiffreadVolume(files(f));
        imDouble = double(im);

        [H,W,Z] = size(im);

        thresh = prctile(imDouble(:), 99); % Protein addition = 98, WTdelmshA = 99

        % Alternative per-experiment-type threshold rules used on past
        % datasets -- swap in whichever applies before running on a new
        % experiment type (only one of these, or the line above, should
        % be active at a time):

        % Spike in QS mutants
        % if f == 1
        %     thresh = prctile(imDouble(:), 99.8);
        % else
        %     thresh = prctile(imDouble(:), 98.5);
        % end

        % Coculture
        % if f == 1
        %     thresh = prctile(imDouble(:), 99.8);
        % else
        %     thresh = prctile(imDouble(:), 99);
        % end

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
        bg = median(imDouble(~bin));

        selectedGray = zeros(H,W, "like", im);
        selectedGray(:) = bg;

        % Same per-3D-object triage as the coculture script: keep only
        % objects clearing minObjectArea, drop likely motion artifacts
        % (<3 slices, or over-spread in max-projection), and for
        % everything kept, project just its single best (largest-area)
        % slice into the flat selectedBin/selectedGray images.
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

%% CHOOSE SMALL-OBJECT IMAGE AND LARGE-OBJECT IMAGE
% For a genuine pair, small/large are picked by channel name
% (smallChannel/largeChannel). For a single-file replicate, smallIdx ==
% largeIdx == 1 -- the same image supplies both roles, so "large objects"
% become whatever big blobs exist within that same channel (used as an
% internal distance-reference surface, e.g. microcolonies within a
% monoculture image).

    smallIdx = 1;
    largeIdx = 1;

    if nFiles == 2

        smallIdx = find(contains(files, smallChannel), 1);
        largeIdx = find(contains(files, largeChannel), 1);

    end

    smallBinImage = selectedBinCell{smallIdx};
    largeBinImage = selectedBinCell{largeIdx};

%% FIND LARGE OBJECTS IN LARGE-OBJECT IMAGE

    CClarge = bwconncomp(largeBinImage, 8);
    propsLarge = regionprops(CClarge, "Area", "PixelIdxList");

    largeMask = false(size(largeBinImage));

    for j = 1:length(propsLarge)

        if propsLarge(j).Area > largeAreaThresh
            largeMask(propsLarge(j).PixelIdxList) = true;
        end

    end

    if any(largeMask(:))
        distToLarge = bwdist(largeMask);
    else
        distToLarge = nan(size(largeMask));
    end

%% WATERSHED SPLIT TEST ON SMALL-OBJECT IMAGE
% Same distance-transform + h-maxima + watershed recipe as elsewhere in
% this pipeline family; ridge flags which pixels were the dividing line
% between two objects that had to be separated (WasSplit below).

    D = bwdist(~smallBinImage);
    D = imhmax(D, 1);
    Lw = watershed(-D);
    ridge = Lw == 0;
    ridge = ridge & smallBinImage;

%% CLASSIFY OBJECTS IN SMALL-OBJECT IMAGE

    CCsmall = bwconncomp(smallBinImage, 8);
    propsSmall = regionprops(CCsmall, "Area", "Centroid", "PixelIdxList", ...
        "MajorAxisLength", "MinorAxisLength");

    objectID = [];
    objectArea = [];
    centroidX = [];
    centroidY = [];
    aspectRatio = [];
    wasSplit = [];
    watershedResult = strings(0,1);
    objectClass = strings(0,1);
    nearestLargeDist_pixels = [];

    singleMask = false(size(smallBinImage));
    otherSmallMask = false(size(smallBinImage));

    isSingleChannelCase = smallIdx == largeIdx;

    for j = 1:length(propsSmall)

        includeThisSmallObject = propsSmall(j).Area >= singleMinArea;

        % Single-channel case only: the same connected-component pass
        % that finds small objects would also pick up the large blobs
        % themselves (since small/large images are literally the same
        % image here) -- drop anything big enough to already be counted
        % as "large" so it isn't also double-counted as a small object.
        if isSingleChannelCase && propsSmall(j).Area > largeAreaThresh
            includeThisSmallObject = false;
        end

        % Exclude small objects that lie completely within a large (red)
        % object - they are pixels enclosed by a large object's footprint,
        % not standalone small objects, even though they passed connected-
        % component detection in the small-object channel.
        if includeThisSmallObject && all(largeMask(propsSmall(j).PixelIdxList))
            includeThisSmallObject = false;
        end

        if includeThisSmallObject

            c = propsSmall(j).Centroid;

            if propsSmall(j).MinorAxisLength == 0
                ar = inf;
            else
                ar = propsSmall(j).MajorAxisLength / propsSmall(j).MinorAxisLength;
            end

            splitNow = any(ridge(propsSmall(j).PixelIdxList));

            isSingleCell = propsSmall(j).Area > singleMinArea && ...
                           propsSmall(j).Area < singleMaxArea && ...
                           ar < singleMaxAspectRatio && ...
                           splitNow == false;

            if isSingleCell
                classNow = "SingleCell";
                singleMask(propsSmall(j).PixelIdxList) = true;
            else
                classNow = "GreaterThanSingleCell";
                otherSmallMask(propsSmall(j).PixelIdxList) = true;
            end

            if splitNow
                watershedNow = "Split";
            else
                watershedNow = "NotSplit";
            end

            d = interp2(distToLarge, c(1), c(2), "linear");

            objectID = [objectID; j];
            objectArea = [objectArea; propsSmall(j).Area];
            centroidX = [centroidX; c(1)];
            centroidY = [centroidY; c(2)];
            aspectRatio = [aspectRatio; ar];
            wasSplit = [wasSplit; splitNow];
            watershedResult = [watershedResult; watershedNow];
            objectClass = [objectClass; classNow];
            nearestLargeDist_pixels = [nearestLargeDist_pixels; d];

        end

    end

%% BIN DISTANCES FOR THIS REPLICATE

    if any(largeMask(:))
        availableDistVals = distToLarge(~largeMask);
        availableDistVals(isinf(availableDistVals)) = nan;
        maxDist = max(availableDistVals, [], "omitnan");
    else
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

    distanceBin = strings(length(nearestLargeDist_pixels),1);

    for j = 1:length(nearestLargeDist_pixels)

        b = find(nearestLargeDist_pixels(j) >= binStart & nearestLargeDist_pixels(j) < binEnd, 1);

        if ~isempty(b)
            distanceBin(j) = binLabel(b);
        else
            distanceBin(j) = "NoLargeObject";
        end

    end

    numSingleCells = zeros(length(binStart),1);
    numGreaterThanSingleCell = zeros(length(binStart),1);
    totalObjects = zeros(length(binStart),1);
    propSingleCells = zeros(length(binStart),1);
    propGreaterThanSingleCell = zeros(length(binStart),1);

    pixelsInBin = zeros(length(binStart),1);
    densitySingleCells = zeros(length(binStart),1);
    densityGreaterThanSingleCell = zeros(length(binStart),1);
    densityTotalObjects = zeros(length(binStart),1);

    for j = 1:length(binStart)

        inBin = nearestLargeDist_pixels >= binStart(j) & nearestLargeDist_pixels < binEnd(j);

        numSingleCells(j) = sum(inBin & objectClass == "SingleCell");
        numGreaterThanSingleCell(j) = sum(inBin & objectClass == "GreaterThanSingleCell");
        totalObjects(j) = numSingleCells(j) + numGreaterThanSingleCell(j);

        availablePixels = distToLarge >= binStart(j) & ...
                          distToLarge < binEnd(j) & ...
                          ~largeMask;

        pixelsInBin(j) = sum(availablePixels(:));

        if totalObjects(j) > 0
            propSingleCells(j) = numSingleCells(j) / totalObjects(j);
            propGreaterThanSingleCell(j) = numGreaterThanSingleCell(j) / totalObjects(j);
        else
            propSingleCells(j) = NaN;
            propGreaterThanSingleCell(j) = NaN;
        end

        if pixelsInBin(j) > 0
            densitySingleCells(j) = numSingleCells(j) / pixelsInBin(j) * densityScale;
            densityGreaterThanSingleCell(j) = numGreaterThanSingleCell(j) / pixelsInBin(j) * densityScale;
            densityTotalObjects(j) = totalObjects(j) / pixelsInBin(j) * densityScale;
        else
            densitySingleCells(j) = NaN;
            densityGreaterThanSingleCell(j) = NaN;
            densityTotalObjects(j) = NaN;
        end

    end

%% SAVE CSVs FOR THIS REPLICATE

    rawTable = table(objectID, objectArea, centroidX, centroidY, aspectRatio, ...
        wasSplit, watershedResult, objectClass, nearestLargeDist_pixels, distanceBin);

    rawTable.Condition = repmat(conditionName, height(rawTable), 1);
    rawTable.Replicate = repmat(replicateNum, height(rawTable), 1);
    rawTable.SourceFile = repmat(files(smallIdx), height(rawTable), 1);
    rawTable = movevars(rawTable, {'Condition','Replicate','SourceFile'}, 'Before', 1);

    binTable = table(binLabel, binStart, binEnd, pixelsInBin, ...
        numSingleCells, numGreaterThanSingleCell, totalObjects, ...
        propSingleCells, propGreaterThanSingleCell, ...
        densitySingleCells, densityGreaterThanSingleCell, densityTotalObjects, ...
        'VariableNames', {'DistanceBin','BinStart_pixels','BinEnd_pixels', ...
        'PixelsInBin','NumSingleCells','NumGreaterThanSingleCellObjects', ...
        'TotalObjects','ProportionSingleCells', ...
        'ProportionGreaterThanSingleCellObjects', ...
        'DensitySingleCells_per10000px', ...
        'DensityGreaterThanSingleCellObjects_per10000px', ...
        'DensityTotalObjects_per10000px'});

    binTable.Condition = repmat(conditionName, height(binTable), 1);
    binTable.Replicate = repmat(replicateNum, height(binTable), 1);
    binTable.SourceFile = repmat(files(smallIdx), height(binTable), 1);
    binTable = movevars(binTable, {'Condition','Replicate','SourceFile'}, 'Before', 1);

    totalIncludedObjects = length(objectID);
    totalSingleCells = sum(objectClass == "SingleCell");
    totalGreaterThanSingleCell = sum(objectClass == "GreaterThanSingleCell");

    if totalIncludedObjects > 0
        proportionSingleCells = totalSingleCells / totalIncludedObjects;
        proportionGreaterThanSingleCell = totalGreaterThanSingleCell / totalIncludedObjects;
    else
        proportionSingleCells = NaN;
        proportionGreaterThanSingleCell = NaN;
    end

    summaryTable = table(conditionName, replicateNum, files(smallIdx), ...
        totalIncludedObjects, totalSingleCells, totalGreaterThanSingleCell, ...
        proportionSingleCells, proportionGreaterThanSingleCell, ...
        'VariableNames', {'Condition','Replicate','SourceFile', ...
        'TotalIncludedObjects','TotalSingleCells', ...
        'TotalGreaterThanSingleCellObjects','ProportionSingleCells', ...
        'ProportionGreaterThanSingleCellObjects'});

    outBase = erase(files(smallIdx), ".tif");

    writetable(rawTable, outBase + "_Objects_RawDistances_Classified.csv")
    writetable(binTable, outBase + "_Objects_BinnedDistances_Classified.csv")

    compiledObjectTable = [compiledObjectTable; rawTable];
    compiledSummaryTable = [compiledSummaryTable; summaryTable];
    compiledBinTableAll = [compiledBinTableAll; binTable];

%% SAVE STACKED COUNT PLOT FOR THIS REPLICATE

    figure

    binCats = categorical(binLabel, binLabel, "Ordinal", true);

    Y = [numSingleCells numGreaterThanSingleCell];

    b = bar(binCats, Y, "stacked");

    b(1).FaceColor = [0 0.7 0];
    b(2).FaceColor = [0 0.2 1];

    xlabel("Distance to nearest large object edge, pixels")
    ylabel("Number of objects")
    title("Objects by distance bin")
    legend(["Single cells", ">single cell objects"], "Location", "best")

    xtickangle(45)

    saveas(gcf, outBase + "_Objects_DistanceBins_StackedCounts.png")
    saveas(gcf, outBase + "_Objects_DistanceBins_StackedCounts.fig")

    close(gcf)

%% SAVE DENSITY-NORMALIZED STACKED PLOT FOR THIS REPLICATE

    figure

    binCats = categorical(binLabel, binLabel, "Ordinal", true);

    Ydensity = [densitySingleCells densityGreaterThanSingleCell];

    b = bar(binCats, Ydensity, "stacked");

    b(1).FaceColor = [0 0.7 0];
    b(2).FaceColor = [0 0.2 1];

    xlabel("Distance to nearest large object edge, pixels")
    ylabel("Objects per 10,000 available pixels")
    title("Objects by distance bin, area-normalized")
    legend(["Single cells", ">single cell objects"], "Location", "best")

    xtickangle(45)

    saveas(gcf, outBase + "_Objects_DistanceBins_StackedDensity.png")
    saveas(gcf, outBase + "_Objects_DistanceBins_StackedDensity.fig")

    close(gcf)

%% SAVE OUTLINE OVERLAYS
% Red outline = large objects, green = SingleCell, blue =
% GreaterThanSingleCell small objects -- drawn on EVERY file in this
% replicate (both channels, if paired), so the classification derived
% from the small-object image is also visible on the large-object image.

    largeBoundary = boundarymask(largeMask);
    singleBoundary = boundarymask(singleMask);
    otherBoundary = boundarymask(otherSmallMask);

    for f = 1:nFiles

        [~,baseName,~] = fileparts(files(f));

        baseGray = mat2gray(selectedGrayCell{f});
        rgb = repmat(baseGray, [1 1 3]);

        red = rgb(:,:,1);
        green = rgb(:,:,2);
        blue = rgb(:,:,3);

        red(largeBoundary) = 1;
        green(largeBoundary) = 0;
        blue(largeBoundary) = 0;

        red(singleBoundary) = 0;
        green(singleBoundary) = 1;
        blue(singleBoundary) = 0;

        red(otherBoundary) = 0;
        green(otherBoundary) = 0;
        blue(otherBoundary) = 1;

        rgb(:,:,1) = red;
        rgb(:,:,2) = green;
        rgb(:,:,3) = blue;

        imwrite(im2uint8(rgb), baseName + "_SelectedSlice_Overlay.tif", "tif")

    end

end

%% COMPILE REPLICATES ACROSS FOLDER: RAW COUNTS BY DISTANCE BIN
% Re-bins every replicate's raw per-object distances onto one shared set
% of bin edges (spanning the largest distance seen across the whole
% compiled dataset, not just one replicate), then averages counts across
% replicates within each condition for a mean +/- SEM stacked-count plot.
% Also writes two Prism-friendly long-format tables alongside the
% MATLAB-oriented ones.

if height(compiledObjectTable) > 0

    writetable(compiledObjectTable, "Compiled_AllObjects_RawDistances_Classified.csv")

    compiledDist = compiledObjectTable.nearestLargeDist_pixels;
    compiledDist(isinf(compiledDist)) = nan;

    maxDistCompiled = max(compiledDist, [], "omitnan");

    if isempty(maxDistCompiled) || isnan(maxDistCompiled)
        maxDistCompiled = 0;
    end

    compiledEdges = 0:binWidth:(ceil(maxDistCompiled/binWidth)*binWidth + binWidth);

    if length(compiledEdges) < 2
        compiledEdges = [0 binWidth];
    end

    compiledBinStart = compiledEdges(1:end-1)';
    compiledBinEnd = compiledEdges(2:end)';
    compiledBinLabel = strings(length(compiledBinStart),1);

    for j = 1:length(compiledBinStart)
        compiledBinLabel(j) = string(compiledBinStart(j)) + "-" + string(compiledBinEnd(j));
    end

    conditions = unique(compiledObjectTable.Condition, "stable");

    compiledBinTable = table();
    compiledPrismLong = table();
    compiledStatsTable = table();

    for c = 1:length(conditions)

        thisCondition = conditions(c);
        conditionRows = compiledObjectTable.Condition == thisCondition;
        reps = unique(compiledObjectTable.Replicate(conditionRows), "stable");

        repBinSingle = zeros(length(reps), length(compiledBinStart));
        repBinGreater = zeros(length(reps), length(compiledBinStart));
        repBinTotal = zeros(length(reps), length(compiledBinStart));

        for r = 1:length(reps)

            thisRep = reps(r);
            repRows = conditionRows & compiledObjectTable.Replicate == thisRep;

            repDistances = compiledObjectTable.nearestLargeDist_pixels(repRows);
            repClasses = compiledObjectTable.objectClass(repRows);

            for bnum = 1:length(compiledBinStart)

                inBin = repDistances >= compiledBinStart(bnum) & ...
                        repDistances < compiledBinEnd(bnum);

                repBinSingle(r,bnum) = sum(inBin & repClasses == "SingleCell");
                repBinGreater(r,bnum) = sum(inBin & repClasses == "GreaterThanSingleCell");
                repBinTotal(r,bnum) = repBinSingle(r,bnum) + repBinGreater(r,bnum);

                if repBinTotal(r,bnum) > 0
                    pSingle = repBinSingle(r,bnum) / repBinTotal(r,bnum);
                    pGreater = repBinGreater(r,bnum) / repBinTotal(r,bnum);
                else
                    pSingle = NaN;
                    pGreater = NaN;
                end

                oneBinRow = table(thisCondition, thisRep, compiledBinLabel(bnum), ...
                    compiledBinStart(bnum), compiledBinEnd(bnum), ...
                    repBinSingle(r,bnum), repBinGreater(r,bnum), ...
                    repBinTotal(r,bnum), pSingle, pGreater, ...
                    'VariableNames', {'Condition','Replicate','DistanceBin', ...
                    'BinStart_pixels','BinEnd_pixels','NumSingleCells', ...
                    'NumGreaterThanSingleCellObjects','TotalObjects', ...
                    'ProportionSingleCells','ProportionGreaterThanSingleCellObjects'});

                compiledBinTable = [compiledBinTable; oneBinRow];

                prismRow = table(thisCondition, compiledBinLabel(bnum), thisRep, ...
                    repBinSingle(r,bnum), repBinGreater(r,bnum), ...
                    'VariableNames', {'Condition','DistanceBin','Replicate', ...
                    'SingleCells','GreaterThanSingleCellObjects'});

                compiledPrismLong = [compiledPrismLong; prismRow];

            end

        end

        meanSingle = mean(repBinSingle, 1, "omitnan")';
        meanGreater = mean(repBinGreater, 1, "omitnan")';
        meanTotal = mean(repBinTotal, 1, "omitnan")';

        semSingle = std(repBinSingle, 0, 1, "omitnan")' ./ sqrt(size(repBinSingle,1));
        semGreater = std(repBinGreater, 0, 1, "omitnan")' ./ sqrt(size(repBinGreater,1));
        semTotal = std(repBinTotal, 0, 1, "omitnan")' ./ sqrt(size(repBinTotal,1));

        statsCondition = repmat(thisCondition, length(compiledBinLabel), 1);

        statsTableThis = table(statsCondition, compiledBinLabel, compiledBinStart, ...
            compiledBinEnd, meanSingle, semSingle, meanGreater, semGreater, ...
            meanTotal, semTotal, ...
            'VariableNames', {'Condition','DistanceBin','BinStart_pixels', ...
            'BinEnd_pixels','MeanSingleCells','SEMSingleCells', ...
            'MeanGreaterThanSingleCellObjects','SEMGreaterThanSingleCellObjects', ...
            'MeanTotalObjects','SEMTotalObjects'});

        compiledStatsTable = [compiledStatsTable; statsTableThis];

        figure

        binCats = categorical(compiledBinLabel, compiledBinLabel, "Ordinal", true);

        Ymean = [meanSingle meanGreater];

        b = bar(binCats, Ymean, "stacked");

        b(1).FaceColor = [0 0.7 0];
        b(2).FaceColor = [0 0.2 1];

        hold on

        x = 1:length(compiledBinLabel);

        errorbar(x, meanSingle, semSingle, "k", "LineStyle", "none", "LineWidth", 1)
        errorbar(x, meanTotal, semTotal, "k", "LineStyle", "none", "LineWidth", 1)

        hold off

        xlabel("Distance to nearest large object edge, pixels")
        ylabel("Mean number of objects")
        title(thisCondition + " compiled replicate summary")
        legend(["Single cells", ">single cell objects"], "Location", "best")

        xtickangle(45)

        outCondition = regexprep(thisCondition, '[^\w\-]', '_');

        saveas(gcf, outCondition + "_Compiled_Objects_DistanceBins_StackedCounts_MeanSEM.png")
        saveas(gcf, outCondition + "_Compiled_Objects_DistanceBins_StackedCounts_MeanSEM.fig")

        close(gcf)

    end

    writetable(compiledBinTable, "Compiled_ReplicateBins_ForPrism_Long.csv")
    writetable(compiledPrismLong, "Compiled_Prism_StackedBar_Long.csv")
    writetable(compiledStatsTable, "Compiled_Stats_MeanSEM_ByDistanceBin.csv")

end

%% COMPILE DENSITY-NORMALIZED DISTANCE BINS ACROSS REPLICATES
% Same cross-replicate averaging idea as above, but on each replicate's
% OWN already-binned density table (compiledBinTableAll) rather than
% re-binning raw distances -- density values are only comparable when
% pooled at matching bin labels, which each replicate already used a
% shared binWidth to produce.

if height(compiledBinTableAll) > 0

    writetable(compiledBinTableAll, "Compiled_ReplicateBins_CountsAndDensity.csv")

    conditionsDensity = unique(compiledBinTableAll.Condition, "stable");

    compiledDensityStats = table();

    for c = 1:length(conditionsDensity)

        thisCondition = conditionsDensity(c);
        conditionRows = compiledBinTableAll.Condition == thisCondition;

        binsNow = unique(compiledBinTableAll.DistanceBin(conditionRows), "stable");

        for bnum = 1:length(binsNow)

            thisBin = binsNow(bnum);
            rowsNow = conditionRows & compiledBinTableAll.DistanceBin == thisBin;

            valsSingle = compiledBinTableAll.DensitySingleCells_per10000px(rowsNow);
            valsGreater = compiledBinTableAll.DensityGreaterThanSingleCellObjects_per10000px(rowsNow);
            valsTotal = compiledBinTableAll.DensityTotalObjects_per10000px(rowsNow);

            nReps = sum(~isnan(valsTotal));

            meanSingleDensity = mean(valsSingle, "omitnan");
            meanGreaterDensity = mean(valsGreater, "omitnan");
            meanTotalDensity = mean(valsTotal, "omitnan");

            semSingleDensity = std(valsSingle, 0, "omitnan") / sqrt(nReps);
            semGreaterDensity = std(valsGreater, 0, "omitnan") / sqrt(nReps);
            semTotalDensity = std(valsTotal, 0, "omitnan") / sqrt(nReps);

            binStartNow = compiledBinTableAll.BinStart_pixels(find(rowsNow, 1));
            binEndNow = compiledBinTableAll.BinEnd_pixels(find(rowsNow, 1));

            oneRow = table(thisCondition, thisBin, binStartNow, binEndNow, nReps, ...
                meanSingleDensity, semSingleDensity, ...
                meanGreaterDensity, semGreaterDensity, ...
                meanTotalDensity, semTotalDensity, ...
                'VariableNames', {'Condition','DistanceBin','BinStart_pixels', ...
                'BinEnd_pixels','NReplicates', ...
                'MeanDensitySingleCells_per10000px', ...
                'SEMDensitySingleCells_per10000px', ...
                'MeanDensityGreaterThanSingleCellObjects_per10000px', ...
                'SEMDensityGreaterThanSingleCellObjects_per10000px', ...
                'MeanDensityTotalObjects_per10000px', ...
                'SEMDensityTotalObjects_per10000px'});

            compiledDensityStats = [compiledDensityStats; oneRow];

        end

        rowsPlot = compiledDensityStats.Condition == thisCondition;
        plotTable = compiledDensityStats(rowsPlot,:);

        [~,ord] = sort(plotTable.BinStart_pixels);
        plotTable = plotTable(ord,:);

        figure

        binCats = categorical(plotTable.DistanceBin, plotTable.DistanceBin, "Ordinal", true);

        Y = [plotTable.MeanDensitySingleCells_per10000px, ...
             plotTable.MeanDensityGreaterThanSingleCellObjects_per10000px];

        b = bar(binCats, Y, "stacked");

        b(1).FaceColor = [0 0.7 0];
        b(2).FaceColor = [0 0.2 1];

        hold on

        x = 1:height(plotTable);

        errorbar(x, plotTable.MeanDensitySingleCells_per10000px, ...
            plotTable.SEMDensitySingleCells_per10000px, ...
            "k", "LineStyle", "none", "LineWidth", 1)

        errorbar(x, plotTable.MeanDensityTotalObjects_per10000px, ...
            plotTable.SEMDensityTotalObjects_per10000px, ...
            "k", "LineStyle", "none", "LineWidth", 1)

        hold off

        xlabel("Distance to nearest large object edge, pixels")
        ylabel("Mean objects per 10,000 available pixels")
        title(thisCondition + " area-normalized distance-bin summary")
        legend(["Single cells", ">single cell objects"], "Location", "best")

        xtickangle(45)

        outCondition = regexprep(thisCondition, '[^\w\-]', '_');

        saveas(gcf, outCondition + "_Compiled_Objects_DistanceBins_StackedDensity_MeanSEM.png")
        saveas(gcf, outCondition + "_Compiled_Objects_DistanceBins_StackedDensity_MeanSEM.fig")

        close(gcf)

    end

    writetable(compiledDensityStats, "Compiled_Stats_MeanSEM_ByDistanceBin_DensityNormalized.csv")

end

%% COMPILE TOTAL SUMMARY ACROSS REPLICATES
% Collapses every replicate down to one number per condition (no
% distance-bin breakdown at all) -- overall mean +/- SEM single-cell vs.
% greater-than-single-cell counts, for a quick top-level comparison.

if height(compiledSummaryTable) > 0

    writetable(compiledSummaryTable, "Compiled_ReplicateSummary_AllReplicates.csv")

    conditionsSummary = unique(compiledSummaryTable.Condition, "stable");

    nReps = zeros(length(conditionsSummary),1);

    meanSingle = zeros(length(conditionsSummary),1);
    semSingle = zeros(length(conditionsSummary),1);

    meanGreater = zeros(length(conditionsSummary),1);
    semGreater = zeros(length(conditionsSummary),1);

    meanTotal = zeros(length(conditionsSummary),1);
    semTotal = zeros(length(conditionsSummary),1);

    for c = 1:length(conditionsSummary)

        thisCondition = conditionsSummary(c);
        rows = compiledSummaryTable.Condition == thisCondition;

        singleVals = compiledSummaryTable.TotalSingleCells(rows);
        greaterVals = compiledSummaryTable.TotalGreaterThanSingleCellObjects(rows);
        totalVals = compiledSummaryTable.TotalIncludedObjects(rows);

        nReps(c) = sum(rows);

        meanSingle(c) = mean(singleVals, "omitnan");
        meanGreater(c) = mean(greaterVals, "omitnan");
        meanTotal(c) = mean(totalVals, "omitnan");

        semSingle(c) = std(singleVals, 0, "omitnan") / sqrt(nReps(c));
        semGreater(c) = std(greaterVals, 0, "omitnan") / sqrt(nReps(c));
        semTotal(c) = std(totalVals, 0, "omitnan") / sqrt(nReps(c));

    end

    totalStatsTable = table(conditionsSummary, nReps, ...
        meanSingle, semSingle, meanGreater, semGreater, meanTotal, semTotal, ...
        'VariableNames', {'Condition','NReplicates', ...
        'MeanSingleCells','SEMSingleCells', ...
        'MeanGreaterThanSingleCellObjects','SEMGreaterThanSingleCellObjects', ...
        'MeanTotalObjects','SEMTotalObjects'});

    writetable(totalStatsTable, "Compiled_TotalSummary_MeanSEM_ByCondition.csv")

    prismTotalTable = table(compiledSummaryTable.Condition, ...
        compiledSummaryTable.Replicate, ...
        compiledSummaryTable.TotalSingleCells, ...
        compiledSummaryTable.TotalGreaterThanSingleCellObjects, ...
        compiledSummaryTable.TotalIncludedObjects, ...
        'VariableNames', {'Condition','Replicate','SingleCells', ...
        'GreaterThanSingleCellObjects','TotalObjects'});

    writetable(prismTotalTable, "Compiled_Prism_TotalSummary_Replicates.csv")

    figure

    conditionCats = categorical(conditionsSummary, conditionsSummary, "Ordinal", true);

    Y = [meanSingle meanGreater];

    b = bar(conditionCats, Y, "stacked");

    b(1).FaceColor = [0 0.7 0];
    b(2).FaceColor = [0 0.2 1];

    hold on

    x = 1:length(conditionsSummary);

    errorbar(x, meanSingle, semSingle, "k", "LineStyle", "none", "LineWidth", 1)
    errorbar(x, meanTotal, semTotal, "k", "LineStyle", "none", "LineWidth", 1)

    hold off

    ylabel("Mean number of objects")
    title("Single-cell classification summary")
    legend(["Single cells", ">single cell objects"], "Location", "best")

    xtickangle(45)

    saveas(gcf, "Compiled_TotalSummary_StackedBar_MeanSEM.png")
    saveas(gcf, "Compiled_TotalSummary_StackedBar_MeanSEM.fig")

    close(gcf)

end
