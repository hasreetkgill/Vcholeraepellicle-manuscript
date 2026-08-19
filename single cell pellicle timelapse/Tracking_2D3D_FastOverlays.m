%% Tracking_2D3D_FastOverlays.m
% Author: Hasreet Gill, Bassler Lab, HHMI/Princeton University
% Developed/revised with Claude (Sonnet 5, Anthropic) - AI-assisted code review and comments.

% Per-tile object tracking across a timelapse, built on top of the MAX
% projection movies and per-timepoint 3D binary stacks written by the
% earlier segmentation/registration pipeline (ProcessSurfaceTiles_*).
%
% For each tile:
%   1. Loads that tile's time-lapse MAX-projection binary and does a
%      first-pass frame-to-frame link of every 2D object, using 3D
%      connectivity across TIME (stacking two consecutive frames and
%      checking if an object's footprint touches/overlaps between them)
%      plus a "didn't shrink too much" area rule.
%   2. Lets the user pick a track-length range (after inspecting one
%      frame's track lengths) to flag the long-lived tracks -- these are
%      trusted enough to use as local motion references, the same role
%      beads play in the acquisition-time drift correction, but derived
%      from the biology itself and spatially local rather than one global
%      shift.
%   3. Computes a per-tile "anchor" frame (the median start frame of the
%      selected long-lived tracks) and rebuilds frame-to-frame links
%      BEFORE that anchor from scratch via a more permissive backward
%      chase (nearest sufficiently-large previous-frame object), since
%      first-pass tracking is least reliable before enough of these
%      objects are established and distinguishable.
%   4. Builds a dense per-frame XY displacement field from the confirmed
%      long-track object motions (scattered interpolation with a global
%      median-shift fallback when there are too few references), then
%      uses it to cumulatively register every frame into the anchor
%      frame's coordinate system (pre-anchor frames chained forward,
%      post-anchor frames chained backward).
%   5. Re-runs the same 3D-connectivity linking a second time on these
%      registered coordinates for a cleaner final track set, additionally
%      detecting "merge" events (multiple objects converging into one) for
%      lineage bookkeeping.
%   6. Matches every final 2D (MAX-projection) tracked object back to its
%      corresponding object in that timepoint's full 3D stack (nearest
%      centroid within maxXYMatchDist) to pull in real calibrated
%      Volume/Z-centroid measurements.
%   7. Writes, per tile: an Excel workbook (Objects + TrackLineage sheets)
%      plus diagnostic overlay/vector-field TIFFs for visually checking
%      every stage of the correction.

%% Setup

clear
clc
close all

path = uigetdir();

if path == 0
    error("No folder selected.")
end

maxFolder = fullfile(path, "MAX");

if ~exist(maxFolder, "dir")
    error("Could not find MAX subfolder in selected folder.")
end

outputRoot = fullfile(path, "Tracking_Output");

if ~exist(outputRoot, "dir")
    mkdir(outputRoot)
end

% Frame-to-frame linking / vector-field parameters
areaRatioMin = 0.5;      % next object must be at least 50% of current object area
minRefsForInterp = 3;    % interpolate vector field if >= this many refs; 1-2 refs use global median shift
maxDispForDisplay = 20;  % pixels; display scaling for vector maps

lineWidth = 2;
spotRadius = 4;

% Calibration
pixelSizeXY = 0.11;   % um/pixel
zStep = 0.5;          % um/slice

pixelArea_um2 = pixelSizeXY^2;
voxelVolume_um3 = pixelSizeXY^2 * zStep;

maxXYMatchDist = 10;  % pixels; 2D-to-3D centroid matching threshold

% Special early-frame correction
% Anchor is computed per tile after long-lived tracks are selected:
%   earlyAnchorFrame = rounded median first frame of selected long-lived tracks
earlyFirstFrame = 1;          % first frame included in early backward extension
earlyAreaRatioMin = 0.5;      % previous-frame object must be > 50% area of target object
earlyMinObjectArea = 200;     % previous-frame object itself must be > this area, in pixels
earlyMaxMatchDist = 100;      % pixels; previous-frame object must be within this distance

%% Find MAX projection input files
% These are the running per-tile time-lapse max-projection TIFFs written
% by the earlier registration pipeline (one frame appended per timepoint).

maxFiles = dir(fullfile(maxFolder, "*Tile*_ShBin_NoNoiseSwim_Resc_StDrCorr_NoShaky_ZAl_MAX*.tif"));

% Excludes this script's own output naming patterns, in case an old run's
% outputs ended up alongside the MAX files (e.g. MAX subfolder reused) --
% otherwise a stale overlay/vector-field file could get picked up as if
% it were a real MAX-projection input.
excludeTerms = ["VectorField", "VectorRefs", "ShiftedBinary", ...
                "ShiftedNext", "TrackedObjects", "Overlay"];

keep = true(length(maxFiles), 1);

for k = 1:length(maxFiles)

    fname = string(maxFiles(k).name);

    for e = 1:length(excludeTerms)
        if contains(fname, excludeTerms(e))
            keep(k) = false;
        end
    end
end

maxFiles = maxFiles(keep);

if isempty(maxFiles)
    error("No valid MAX projection input files found.")
end

maxTileNums = nan(length(maxFiles), 1);

for k = 1:length(maxFiles)

    tok = regexp(maxFiles(k).name, "Tile(\d+)", "tokens", "once");

    if isempty(tok)
        error("Could not parse tile number from MAX file: " + string(maxFiles(k).name))
    end

    maxTileNums(k) = str2double(tok{1});
end

[~, tileOrder] = sort(maxTileNums);
maxFiles = maxFiles(tileOrder);
maxTileNums = maxTileNums(tileOrder);

%% Find all full 3D stack files in parent folder
% These are the final per-timepoint, per-tile aligned 3D binary stacks
% (also from the earlier registration pipeline) -- used later to look up
% each tracked 2D object's real 3D volume and Z position.

surfFilesAll = dir(fullfile(path, "Surf_*_Export_*_ShBin_NoNoiseSwim_Resc_StDrCorr_NoShaky_ZAl*.tif"));

if isempty(surfFilesAll)
    error("No Surf_*_Export_* full 3D stack files found in parent folder.")
end

surfNums = [];
surfTileNums = [];
surfNames = strings(0,1);
surfFullPaths = strings(0,1);

for k = 1:length(surfFilesAll)

    fname = string(surfFilesAll(k).name);

    tok = regexp(fname, ...
    "^Surf_(\d+)_Export_(\d+)_ShBin_NoNoiseSwim_Resc_StDrCorr_NoShaky_ZAl(?:_Shift)?\.tif$", ...
    "tokens", "once");

    if isempty(tok)
        continue
    end

    surfNums(end+1,1) = str2double(tok{1});
    surfTileNums(end+1,1) = str2double(tok{2});
    surfNames(end+1,1) = fname;
    surfFullPaths(end+1,1) = string(fullfile(path, surfFilesAll(k).name));
end

%% Prefer _Shift full 3D files if duplicate SurfNumber/Tile pairs exist
% If both:
%   Surf_X_Export_Y_..._ZAl.tif
%   Surf_X_Export_Y_..._ZAl_Shift.tif
% exist for the same SurfNumber and tile, keep the _Shift version only.

if ~isempty(surfNums)

    surfKeys = string(surfNums) + "_" + string(surfTileNums);
    uniqueKeys = unique(surfKeys, "stable");
    keepSurf = false(length(surfNums), 1);

    for k = 1:length(uniqueKeys)

        idx = find(surfKeys == uniqueKeys(k));

        shiftIdx = idx(contains(surfNames(idx), "_Shift.tif"));

        if ~isempty(shiftIdx)
            keepSurf(shiftIdx(end)) = true;
        else
            keepSurf(idx(1)) = true;
        end
    end

    surfNums = surfNums(keepSurf);
    surfTileNums = surfTileNums(keepSurf);
    surfNames = surfNames(keepSurf);
    surfFullPaths = surfFullPaths(keepSurf);
end

%% Loop through tiles

for fileIdx = 1:length(maxFiles)

    tileNum = maxTileNums(fileIdx);
    tileName = "Tile" + tileNum;

    disp("==============================================")
    disp("Processing " + tileName)
    disp("==============================================")

    inFile = fullfile(maxFolder, maxFiles(fileIdx).name);
    [~, baseName, ~] = fileparts(inFile);
    name = string(baseName);

    tileOutFolder = fullfile(outputRoot, tileName);

    if ~exist(tileOutFolder, "dir")
        mkdir(tileOutFolder)
    end

    %% Match full 3D files for this tile

    tileMask = surfTileNums == tileNum;

    tileSurfNums = surfNums(tileMask);
    tileSurfNames = surfNames(tileMask);
    tileSurfFullPaths = surfFullPaths(tileMask);

    [tileSurfNums, sortIdx] = sort(tileSurfNums);
    tileSurfNames = tileSurfNames(sortIdx);
    tileSurfFullPaths = tileSurfFullPaths(sortIdx);

    if isempty(tileSurfFullPaths)
        error("No full 3D stack files found for " + tileName)
    end

    %% Output files for this tile

    outRGB  = fullfile(tileOutFolder, name + "_VectorField_RGB.tif");
    outMag  = fullfile(tileOutFolder, name + "_VectorField_Magnitude.tif");
    outDx   = fullfile(tileOutFolder, name + "_VectorField_Dx.tif");
    outDy   = fullfile(tileOutFolder, name + "_VectorField_Dy.tif");
    outRefs = fullfile(tileOutFolder, name + "_VectorRefs_Overlay.tif");

    outShiftedBinary = fullfile(tileOutFolder, name + "_ShiftedBinary_Check.tif");
    outShiftedTracks = fullfile(tileOutFolder, name + "_ShiftedNext_Tracks_Overlay.tif");
    outFirstPassLongTracks = fullfile(tileOutFolder, name + "_FirstPassLongTracks_Overlay.tif");
    outExcel = fullfile(tileOutFolder, name + "_TrackedObjects_2D3D.xlsx");

    deleteThese = [outRGB outMag outDx outDy outRefs ...
                   outShiftedBinary outShiftedTracks outFirstPassLongTracks outExcel];

    for q = 1:length(deleteThese)
        if exist(deleteThese(q), "file")
            delete(deleteThese(q))
        end
    end

    %% Read MAX image and find objects in each frame
    % objPixels/objCentroids/objAreas are indexed [frame]{object} and
    % reused throughout the rest of this tile's processing as the
    % canonical per-frame 2D object list.

    im = tiffreadVolume(inFile);
    im = im > 0;

    [h, w, nFrames] = size(im);

    if length(tileSurfFullPaths) ~= nFrames

        disp("WARNING for " + tileName + ": number of full 3D files does not match MAX frames.")
        disp("3D files: " + length(tileSurfFullPaths) + ", MAX frames: " + nFrames)

        nUseFrames = min(nFrames, length(tileSurfFullPaths));

        disp("Limiting analysis to first " + nUseFrames + " matched frames.")

        im = im(:,:,1:nUseFrames);
        tileSurfNums = tileSurfNums(1:nUseFrames);
        tileSurfNames = tileSurfNames(1:nUseFrames);
        tileSurfFullPaths = tileSurfFullPaths(1:nUseFrames);

        nFrames = nUseFrames;
    end

    objPixels = cell(nFrames, 1);
    objCentroids = cell(nFrames, 1);
    objAreas = cell(nFrames, 1);

    for i = 1:nFrames

        slice = im(:,:,i);

        props = regionprops(slice, "PixelIdxList", "Centroid", "Area");

        nObj = length(props);

        objPixels{i} = cell(nObj, 1);
        objCentroids{i} = zeros(nObj, 2);
        objAreas{i} = zeros(nObj, 1);

        for j = 1:nObj
            objPixels{i}{j} = props(j).PixelIdxList;
            objCentroids{i}(j,:) = props(j).Centroid;
            objAreas{i}(j) = props(j).Area;
        end
    end

    %% Track using same 3D/area rules, and record only successful links for vectors
    % Links object j in frame i to object k in frame i+1 by stacking the two
    % frames' binary masks along a third (time) dimension and running a 3D
    % connected-components pass across them -- the same trick used to link
    % objects across Z-slices in the acquisition pipeline, just applied
    % across consecutive TIME frames instead. Any current/next object pair
    % that is 26-connected within that 2-frame stack is a link candidate;
    % among candidates in the same connected group, current and next
    % objects are matched biggest-to-biggest (curObjs/nxtObjs sorted by
    % area) so a small fragment doesn't steal a match meant for a bigger
    % object, and a match is only accepted if the next object didn't shrink
    % below areaRatioMin of the current object's area (guards against
    % matching onto an unrelated small neighboring object).
    % objectLinked/prevXY record ONLY successful links -- later used to
    % pick trustworthy motion samples for the vector field, since an
    % object that failed to link (a new track) has no real displacement
    % to offer.

    trackIDs = cell(nFrames, 1);
    objectAge = cell(nFrames, 1);
    objectLinked = cell(nFrames, 1);
    prevXY = cell(nFrames, 1);

    nextTrackID = 0;
    trackAge = [];

    curCC = bwconncomp(im(:,:,1));
    nObj = curCC.NumObjects;

    trackIDs{1} = zeros(nObj, 1);
    objectAge{1} = ones(nObj, 1);
    objectLinked{1} = false(nObj, 1);
    prevXY{1} = nan(nObj, 2);

    for j = 1:nObj
        nextTrackID = nextTrackID + 1;
        trackIDs{1}(j) = nextTrackID;
        trackAge(nextTrackID,1) = 1;
    end

    for i = 1:nFrames-1

        disp("Tracking frame " + i + " / " + (nFrames-1))
        tic

        curSlice = im(:,:,i);
        nxtSlice = im(:,:,i+1);

        curCC = bwconncomp(curSlice);
        nxtCC = bwconncomp(nxtSlice);

        curLabel = labelmatrix(curCC);
        nxtLabel = labelmatrix(nxtCC);

        curProps = regionprops(curCC, "Area");
        nxtProps = regionprops(nxtCC, "Area");

        nCur = curCC.NumObjects;
        nNxt = nxtCC.NumObjects;

        idsCur = trackIDs{i};
        idsNext = zeros(nNxt, 1);

        objectAge{i+1} = zeros(nNxt, 1);
        objectLinked{i+1} = false(nNxt, 1);
        prevXY{i+1} = nan(nNxt, 2);

        curArea = zeros(nCur, 1);
        nxtArea = zeros(nNxt, 1);

        for j = 1:nCur
            curArea(j) = curProps(j).Area;
        end

        for j = 1:nNxt
            nxtArea(j) = nxtProps(j).Area;
        end

        curUsed = false(nCur, 1);

        twoFrameStack = cat(3, curSlice, nxtSlice);
        CC3D = bwconncomp(twoFrameStack, 26);

        for g = 1:CC3D.NumObjects

            pix3D = CC3D.PixelIdxList{g};
            [y, x, z] = ind2sub(size(twoFrameStack), pix3D);

            idxCurPix = sub2ind(size(curSlice), y(z == 1), x(z == 1));
            curObjs = unique(curLabel(idxCurPix));
            curObjs(curObjs == 0) = [];

            idxNxtPix = sub2ind(size(nxtSlice), y(z == 2), x(z == 2));
            nxtObjs = unique(nxtLabel(idxNxtPix));
            nxtObjs(nxtObjs == 0) = [];

            if isempty(curObjs) || isempty(nxtObjs)
                continue
            end

            [~, curOrder] = sort(curArea(curObjs), "descend");
            [~, nxtOrder] = sort(nxtArea(nxtObjs), "descend");

            curObjs = curObjs(curOrder);
            nxtObjs = nxtObjs(nxtOrder);

            nPairs = min(length(curObjs), length(nxtObjs));

            for p = 1:nPairs

                cObj = curObjs(p);
                nObj = nxtObjs(p);

                if curUsed(cObj) || idsNext(nObj) ~= 0
                    continue
                end

                pairStack = false(size(twoFrameStack));
                pairStack(:,:,1) = curLabel == cObj;
                pairStack(:,:,2) = nxtLabel == nObj;

                pairCC = bwconncomp(pairStack, 26);
                same3DObject = pairCC.NumObjects == 1;

                areaOK = nxtArea(nObj) >= areaRatioMin * curArea(cObj);

                if same3DObject && areaOK

                    tid = idsCur(cObj);

                    idsNext(nObj) = tid;
                    curUsed(cObj) = true;

                    trackAge(tid) = trackAge(tid) + 1;

                    objectAge{i+1}(nObj) = trackAge(tid);
                    objectLinked{i+1}(nObj) = true;
                    prevXY{i+1}(nObj,:) = objCentroids{i}(cObj,:);
                end
            end
        end

        for j = 1:nNxt
            if idsNext(j) == 0
                nextTrackID = nextTrackID + 1;
                idsNext(j) = nextTrackID;
                trackAge(nextTrackID,1) = 1;
                objectAge{i+1}(j) = 1;
            end
        end

        trackIDs{i+1} = idsNext;

        toc
    end

    %% Inspect track lengths and choose track-length range
    % Shows one representative frame with each visible track's total
    % first-pass length labeled, so the user can pick a length range that
    % isolates the reliable, long-lived objects (e.g. established
    % microcolonies) worth trusting as motion references below -- short
    % tracks are more likely noise, fragments, or transient objects.

    inspectFrame = min(14, nFrames);
    nLabelsToShow = 50;

    nTracks = nextTrackID;

    trackLength = zeros(nTracks, 1);

    for i = 1:nFrames

        ids = trackIDs{i};

        for j = 1:length(ids)
            tid = ids(j);
            trackLength(tid) = trackLength(tid) + 1;
        end
    end

    idsInspect = trackIDs{inspectFrame};
    xyInspect = objCentroids{inspectFrame};

    if isempty(idsInspect)
        defaultMin = max(trackLength);
        defaultMax = max(trackLength);
    else

        lengthsInspect = trackLength(idsInspect);

        [sortedLengths, order] = sort(lengthsInspect, "descend");

        keepN = min(nLabelsToShow, length(order));
        keep = order(1:keepN);

        figure
        imshow(im(:,:,inspectFrame), [])
        hold on

        title(tileName + " frame " + inspectFrame + ": top " + keepN + " track lengths")

        for k = 1:keepN

            objIdx = keep(k);

            x = xyInspect(objIdx,1);
            y = xyInspect(objIdx,2);

            tid = idsInspect(objIdx);
            len = trackLength(tid);

            text(x, y, string(len), ...
                "Color", "blue", ...
                "FontSize", 12, ...
                "FontWeight", "bold", ...
                "HorizontalAlignment", "center")
        end

        hold off

        disp("Top track lengths represented in " + tileName + " frame " + inspectFrame + ":")
        disp(sortedLengths(1:keepN))

        defaultMin = sortedLengths(min(keepN, length(sortedLengths)));
        defaultMax = sortedLengths(1);
    end

    answer = inputdlg( ...
        {"Enter minimum track length to use:", ...
         "Enter maximum track length to use:"}, ...
        "Choose track length range for " + tileName, ...
        [1 50], ...
        {num2str(defaultMin), num2str(defaultMax)});

    if isempty(answer)
        error("No track length range entered.")
    end

    minTrackLength = str2double(answer{1});
    maxTrackLength = str2double(answer{2});

    if isnan(minTrackLength) || isnan(maxTrackLength)
        error("Entered track length range contains a non-number.")
    end

    if minTrackLength > maxTrackLength
        temp = minTrackLength;
        minTrackLength = maxTrackLength;
        maxTrackLength = temp;
    end

    disp("Using track length range for " + tileName + ":")
    disp([minTrackLength maxTrackLength])

    longTrackIDs = find(trackLength >= minTrackLength & trackLength <= maxTrackLength);

    isLongTrack = false(nTracks, 1);
    isLongTrack(longTrackIDs) = true;

    disp("Number of tracks in this length range:")
    disp(length(longTrackIDs))

    %% Compute dynamic early anchor frame from selected long-lived tracks
    % The anchor frame is the rounded median of the first frame in which each
    % selected long-lived track appears. This is computed per tile.

    longTrackFirstFrames = nan(length(longTrackIDs), 1);

    for kk = 1:length(longTrackIDs)

        tid = longTrackIDs(kk);

        for i = 1:nFrames

            if any(trackIDs{i} == tid)
                longTrackFirstFrames(kk) = i;
                break
            end
        end
    end

    validFirstFrames = longTrackFirstFrames(~isnan(longTrackFirstFrames));

    if isempty(validFirstFrames)
        earlyAnchorFrame = min(7, nFrames);
        disp("No valid long-track first frames found; using fallback early anchor frame:")
        disp(earlyAnchorFrame)
    else
        earlyAnchorFrame = round(median(validFirstFrames));
        earlyAnchorFrame = min(max(earlyAnchorFrame, earlyFirstFrame), nFrames);
        disp("Using early anchor frame from median first frame of selected long-lived tracks:")
        disp(earlyAnchorFrame)
    end

    %% Special early-frame backward extension into selected long-lived tracks
    % Completely ignore original first-pass tracking in frames before the anchor.
    %
    % The anchor frame keeps the selected long-lived first-pass TrackIDs.
    % Frames before the anchor are rebuilt from scratch by walking backward:
    %
    %   anchor frame object -> previous frame object -> previous frame object ...
    %
    % Previous-frame candidate rules:
    %   1. not already used by another early chain in that previous frame
    %   2. area > earlyAreaRatioMin * current target area
    %   3. area > earlyMinObjectArea
    %   4. distance <= earlyMaxMatchDist
    %
    % The chosen previous-frame object is assigned the SAME long-lived TrackID.
    % These clean early links are inserted into objectLinked / prevXY so the
    % vector-field section can use them.

    earlyMatch_TargetFrame = [];
    earlyMatch_PrevFrame = [];
    earlyMatch_TrackID = [];
    earlyMatch_TargetObject = [];
    earlyMatch_PrevObject = [];
    earlyMatch_Distance_pixels = [];
    earlyMatch_TargetArea = [];
    earlyMatch_PrevArea = [];
    earlyMatch_Mode = strings(0,1);

    if earlyAnchorFrame > earlyFirstFrame

        disp("Adding CLEAN special early-frame correction links for frames " + ...
             earlyFirstFrame + " to " + earlyAnchorFrame)
        disp("Ignoring original first-pass TrackIDs before anchor frame " + earlyAnchorFrame)

        %% Clear first-pass tracking information before the anchor frame
        % Do NOT clear the anchor frame, because that is where the selected
        % long-lived TrackIDs are used as the starting points.

        for f = earlyFirstFrame:(earlyAnchorFrame - 1)

            nObjF = length(objAreas{f});

            trackIDs{f} = zeros(nObjF, 1);
            objectLinked{f} = false(nObjF, 1);
            prevXY{f} = nan(nObjF, 2);

            if f == 1
                objectAge{f} = ones(nObjF, 1);
            else
                objectAge{f} = zeros(nObjF, 1);
            end
        end

        %% Also clear objectLinked / prevXY for anchor-frame targets
        % Their early link will be rebuilt cleanly from anchor-1.
        % This removes any old first-pass link into the anchor frame.

        nObjAnchor = length(objAreas{earlyAnchorFrame});
        objectLinked{earlyAnchorFrame} = false(nObjAnchor, 1);
        prevXY{earlyAnchorFrame} = nan(nObjAnchor, 2);

        %% Walk backward from the anchor frame

        for targetFrame = earlyAnchorFrame:-1:(earlyFirstFrame + 1)

            prevFrame = targetFrame - 1;

            disp("Clean early backward extension: frame " + prevFrame + ...
                 " -> frame " + targetFrame)

            targetIDs = trackIDs{targetFrame};

            targetAreas = objAreas{targetFrame};
            prevAreas = objAreas{prevFrame};

            targetXY = objCentroids{targetFrame};
            prevXYall = objCentroids{prevFrame};

            nPrev = length(prevAreas);
            usedPrev = false(nPrev, 1);

            % Extend only objects that currently belong to selected long-lived tracks.
            % At the anchor frame, these are original first-pass long-lived tracks.
            % In earlier target frames, these are ONLY objects assigned by this clean
            % early backward extension, because we cleared trackIDs before the anchor.
            validTargetID = targetIDs > 0 & targetIDs <= length(isLongTrack);

            targetIsLong = false(size(targetIDs));
            targetIsLong(validTargetID) = isLongTrack(targetIDs(validTargetID));

            targetObjs = find(targetIsLong);

            if isempty(targetObjs)
                disp("No clean selected long-lived target objects found in frame " + targetFrame)
                continue
            end

            % Process larger target objects first.
            [~, areaOrder] = sort(targetAreas(targetObjs), "descend");
            targetObjs = targetObjs(areaOrder);

            for a = 1:length(targetObjs)

                tgtObj = targetObjs(a);
                tgtTid = targetIDs(tgtObj);

                tgtArea = targetAreas(tgtObj);
                tgtCent = targetXY(tgtObj,:);

                % Candidate previous-frame objects are raw objects only.
                % No previous-frame first-pass TrackID is considered.
                eligiblePrev = find( ...
                    ~usedPrev & ...
                    prevAreas > earlyAreaRatioMin * tgtArea & ...
                    prevAreas > earlyMinObjectArea);

                if isempty(eligiblePrev)
                    continue
                end

                d = sqrt((prevXYall(eligiblePrev,1) - tgtCent(1)).^2 + ...
                         (prevXYall(eligiblePrev,2) - tgtCent(2)).^2);

                withinDist = d <= earlyMaxMatchDist;

                if ~any(withinDist)
                    continue
                end

                eligiblePrev = eligiblePrev(withinDist);
                d = d(withinDist);

                [dMin, bestIdx] = min(d);
                prevObj = eligiblePrev(bestIdx);

                usedPrev(prevObj) = true;

                % Assign the previous-frame object to the SAME long-lived TrackID.
                % This overwrites the cleared early-frame TrackID.
                trackIDs{prevFrame}(prevObj) = tgtTid;

                % Record the clean link from previous frame into target frame.
                % This is what the vector-field code will use.
                objectLinked{targetFrame}(tgtObj) = true;
                prevXY{targetFrame}(tgtObj,:) = prevXYall(prevObj,:);

                earlyMatch_TargetFrame(end+1,1) = targetFrame;
                earlyMatch_PrevFrame(end+1,1) = prevFrame;
                earlyMatch_TrackID(end+1,1) = tgtTid;
                earlyMatch_TargetObject(end+1,1) = tgtObj;
                earlyMatch_PrevObject(end+1,1) = prevObj;
                earlyMatch_Distance_pixels(end+1,1) = dMin;
                earlyMatch_TargetArea(end+1,1) = tgtArea;
                earlyMatch_PrevArea(end+1,1) = prevAreas(prevObj);
                earlyMatch_Mode(end+1,1) = "CleanEarlyBackwardLink";
            end
        end

    else

        disp("Skipping special early-frame correction because earlyAnchorFrame <= earlyFirstFrame.")

    end

    earlyMatchTable = table( ...
        earlyMatch_TargetFrame, ...
        earlyMatch_PrevFrame, ...
        earlyMatch_TrackID, ...
        earlyMatch_TargetObject, ...
        earlyMatch_PrevObject, ...
        earlyMatch_Distance_pixels, ...
        earlyMatch_TargetArea, ...
        earlyMatch_PrevArea, ...
        earlyMatch_Mode);

    disp("Number of clean early-frame links added:")
    disp(height(earlyMatchTable))

    %% Assign fresh TrackIDs to unlinked early-frame objects
    % Pre-anchor objects not claimed by a backward chain get their own new
    % single-appearance TrackIDs so they stay on record (e.g. for counting
    % unlinked one-frame objects later). They are not long-lived tracks,
    % so isLongTrack is padded with false for them.

    nTracksBeforeEarly = nTracks;

    for f = earlyFirstFrame:(earlyAnchorFrame - 1)

        ids = trackIDs{f};

        for j = 1:length(ids)

            if ids(j) == 0

                nextTrackID = nextTrackID + 1;
                trackIDs{f}(j) = nextTrackID;
                trackAge(nextTrackID,1) = 1;
                objectAge{f}(j) = 1;
            end
        end
    end

    nTracks = nextTrackID;

    if nTracks > nTracksBeforeEarly
        isLongTrack(nTracksBeforeEarly+1:nTracks, 1) = false;
    end

    %% Recompute first-pass track lengths after early correction
    % The early clearing changed which frames belong to which track, so the
    % old trackLength is stale. This fresh version is what you'd use for
    % plots like "frame vs. number of one-frame unlinked objects".

    trackLength = zeros(nTracks, 1);

    for i = 1:nFrames

        ids = trackIDs{i};

        for j = 1:length(ids)
            trackLength(ids(j)) = trackLength(ids(j)) + 1;
        end
    end

    %% Write overlay of first-pass selected long-lived tracks after early extension
    % This is only the first-pass selected long-lived tracks, including any early
    % nearest-smaller-object extension into frames before the anchor.

    disp("Writing first-pass long-lived track overlay after early-frame extension.")

    firstPassTrackFrames = cell(nTracks, 1);
    firstPassTrackXY = cell(nTracks, 1);

    for i = 1:nFrames

        ids = trackIDs{i};
        xy = objCentroids{i};

        for j = 1:length(ids)

            tid = ids(j);

            if tid > 0 && tid <= length(isLongTrack) && isLongTrack(tid)
                firstPassTrackFrames{tid}(end+1,1) = i;
                firstPassTrackXY{tid}(end+1,:) = xy(j,:);
            end
        end
    end

    colorsFirstPass = uint8(255 * hsv(max(length(longTrackIDs),1)));

    % Precompute all segments and dots for the selected long-lived tracks,
    % tagged with the frame at which each appears. A segment appears at the
    % frame of its LATER endpoint (same rule as the original fr <= i logic).

    nPtsFP = zeros(length(longTrackIDs), 1);

    for kk = 1:length(longTrackIDs)
        nPtsFP(kk) = size(firstPassTrackFrames{longTrackIDs(kk)}, 1);
    end

    totalSegsFP = sum(max(nPtsFP - 1, 0));
    totalDotsFP = sum(nPtsFP);

    segXYFP = zeros(totalSegsFP, 4);
    segFrameFP = zeros(totalSegsFP, 1);
    segColorFP = zeros(totalSegsFP, 3);

    dotXYFP = zeros(totalDotsFP, 2);
    dotFrameFP = zeros(totalDotsFP, 1);
    dotColorFP = zeros(totalDotsFP, 3);

    segRowFP = 0;
    dotRowFP = 0;

    for kk = 1:length(longTrackIDs)

        tid = longTrackIDs(kk);

        fr = firstPassTrackFrames{tid};
        xy = firstPassTrackXY{tid};

        if isempty(fr)
            continue
        end

        % Make sure points are in frame order.
        [fr, ord] = sort(fr);
        xy = xy(ord,:);

        color = double(colorsFirstPass(kk,:));
        n = length(fr);

        if n >= 2
            idx = segRowFP + (1:n-1);
            segXYFP(idx,:) = [xy(1:n-1,:), xy(2:n,:)];
            segFrameFP(idx) = fr(2:n);
            segColorFP(idx,:) = repmat(color, n-1, 1);
            segRowFP = segRowFP + n - 1;
        end

        idx = dotRowFP + (1:n);
        dotXYFP(idx,:) = xy;
        dotFrameFP(idx) = fr;
        dotColorFP(idx,:) = repmat(color, n, 1);
        dotRowFP = dotRowFP + n;
    end

    for i = 1:nFrames

        rgb = uint8(255 * repmat(im(:,:,i), 1, 1, 3));

        sel = segFrameFP <= i;

        if any(sel)
            rgb = insertShape(rgb, "Line", segXYFP(sel,:), ...
                "Color", segColorFP(sel,:), ...
                "LineWidth", lineWidth);
        end

        dsel = dotFrameFP == i;

        if any(dsel)
            rgb = insertShape(rgb, "FilledCircle", ...
                [dotXYFP(dsel,:), repmat(spotRadius, nnz(dsel), 1)], ...
                "Color", dotColorFP(dsel,:), ...
                "Opacity", 1);
        end

        if i == 1
            imwrite(rgb, outFirstPassLongTracks, "tif", "Compression", "none");
        else
            imwrite(rgb, outFirstPassLongTracks, "tif", ...
                "WriteMode", "append", ...
                "Compression", "none");
        end
    end

    disp("Saved first-pass long-lived track overlay:")
    disp(outFirstPassLongTracks)

    %% Build vector fields from linked objects in chosen track-length range
    % Treats every successfully-linked long-lived-track object as a local
    % motion reference (refX/refY at its position, refDX/refDY its
    % displacement from the previous frame -- including any early-frame
    % backward-extension links added above). With enough references
    % (>= minRefsForInterp), scatteredInterpolant builds a smooth dense
    % displacement field over the whole tile from those sample points
    % (computed on a coarse "scale"-downsampled grid, then upsampled by
    % imresize for speed); with too few points to interpolate reliably,
    % it falls back to one global median shift for the whole frame
    % (mirroring the bead median-shift approach elsewhere, and important
    % for early frames where few long tracks have appeared yet). With no
    % references at all, the field is left at zero (no correction).

    scale = 8;
    [Xsmall, Ysmall] = meshgrid(1:scale:w, 1:scale:h);

    dxMaps = cell(nFrames, 1);
    dyMaps = cell(nFrames, 1);
    nVectorRefs = zeros(nFrames, 1);
    useVectorShift = false(nFrames, 1);

    for i = 1:nFrames

        disp("Writing vector field frame " + i + " / " + nFrames)

        dxMap = zeros(h, w);
        dyMap = zeros(h, w);

        refOverlay = uint8(255 * repmat(im(:,:,i), 1, 1, 3));

        if i > 1

            linkedNow = objectLinked{i};
            idsNow = trackIDs{i};

            xyNow = objCentroids{i};
            xyPrev = prevXY{i};

            refIdx = find(linkedNow & isLongTrack(idsNow));

            if ~isempty(refIdx)

                refX = xyNow(refIdx,1);
                refY = xyNow(refIdx,2);

                refDX = xyNow(refIdx,1) - xyPrev(refIdx,1);
                refDY = xyNow(refIdx,2) - xyPrev(refIdx,2);

                refPts = round([refX refY]);
                [~, uniqueIdx] = unique(refPts, "rows");

                refX = refX(uniqueIdx);
                refY = refY(uniqueIdx);
                refDX = refDX(uniqueIdx);
                refDY = refDY(uniqueIdx);

                % Draw all reference vectors in two batched calls:
                % red lines from previous to current position, then green dots.
                refOverlay = insertShape(refOverlay, "Line", ...
                    [refX - refDX, refY - refDY, refX, refY], ...
                    "Color", [255 0 0], ...
                    "LineWidth", lineWidth);

                refOverlay = insertShape(refOverlay, "FilledCircle", ...
                    [refX, refY, repmat(spotRadius, length(refX), 1)], ...
                    "Color", [0 255 0], ...
                    "Opacity", 1);

                nVectorRefs(i) = length(refX);

                if nVectorRefs(i) >= minRefsForInterp

                    useVectorShift(i) = true;

                    try
                        Fdx = scatteredInterpolant(refX, refY, refDX, "natural", "nearest");
                        Fdy = scatteredInterpolant(refX, refY, refDY, "natural", "nearest");

                        dxSmall = Fdx(Xsmall, Ysmall);
                        dySmall = Fdy(Xsmall, Ysmall);

                        dxMap = imresize(dxSmall, [h w], "bilinear");
                        dyMap = imresize(dySmall, [h w], "bilinear");

                    catch
                        useVectorShift(i) = false;
                        dxMap(:) = 0;
                        dyMap(:) = 0;
                    end

                elseif nVectorRefs(i) > 0

                    % If there are too few references for interpolation, still use
                    % the available information as a global translation for this frame.
                    % This is especially important for the early-frame correction.
                    useVectorShift(i) = true;

                    dxGlobal = median(refDX, "omitnan");
                    dyGlobal = median(refDY, "omitnan");

                    dxMap(:) = dxGlobal;
                    dyMap(:) = dyGlobal;

                else
                    useVectorShift(i) = false;
                    dxMap(:) = 0;
                    dyMap(:) = 0;
                end
            end
        end

        magMap = sqrt(dxMap.^2 + dyMap.^2);
        dxMaps{i} = dxMap;
        dyMaps{i} = dyMap;

        dxShow = (dxMap + maxDispForDisplay) ./ (2 * maxDispForDisplay);
        dxShow = min(max(dxShow, 0), 1);
        dxShow = uint8(255 * dxShow);

        dyShow = (dyMap + maxDispForDisplay) ./ (2 * maxDispForDisplay);
        dyShow = min(max(dyShow, 0), 1);
        dyShow = uint8(255 * dyShow);

        magShow = magMap ./ maxDispForDisplay;
        magShow = min(max(magShow, 0), 1);
        magShow = uint8(255 * magShow);

        angleMap = atan2(dyMap, dxMap);
        hueMap = mod(angleMap, 2*pi) ./ (2*pi);
        satMap = ones(h, w);
        valMap = min(magMap ./ maxDispForDisplay, 1);

        rgb = hsv2rgb(cat(3, hueMap, satMap, valMap));
        rgb = uint8(255 * rgb);

        if i == 1
            imwrite(rgb, outRGB, "tif", "Compression", "none");
            imwrite(magShow, outMag, "tif", "Compression", "none");
            imwrite(dxShow, outDx, "tif", "Compression", "none");
            imwrite(dyShow, outDy, "tif", "Compression", "none");
            imwrite(refOverlay, outRefs, "tif", "Compression", "none");
        else
            imwrite(rgb, outRGB, "tif", "WriteMode", "append", "Compression", "none");
            imwrite(magShow, outMag, "tif", "WriteMode", "append", "Compression", "none");
            imwrite(dxShow, outDx, "tif", "WriteMode", "append", "Compression", "none");
            imwrite(dyShow, outDy, "tif", "WriteMode", "append", "Compression", "none");
            imwrite(refOverlay, outRefs, "tif", "WriteMode", "append", "Compression", "none");
        end
    end

    %% Build anchor-based cumulative shifted binary stack and label maps
    % The anchor frame is earlyAnchorFrame, computed from long-lived tracks.
    %
    % Frames before the anchor are shifted FORWARD toward the anchor.
    % Example if anchor = 7:
    %   frame 6 -> frame 7
    %   frame 5 -> frame 6 -> frame 7
    %
    % Frames after the anchor are shifted BACKWARD toward the anchor.
    % Example if anchor = 7:
    %   frame 8 -> frame 7
    %   frame 9 -> frame 8 -> frame 7
    %
    % These anchor-shifted frames are used for:
    %   1. the shifted binary visual check stack
    %   2. the final cumulative shifted tracking pass

    anchorFrame = min(max(earlyAnchorFrame, 1), nFrames);

    disp("Using cumulative registration anchor frame:")
    disp(anchorFrame)

    cumShiftedSlices = false(h, w, nFrames);
    cumShiftedLabels = cell(nFrames, 1);
    cumObjCentroids = cell(nFrames, 1);

    for i = 1:nFrames

        disp("Building anchor-shifted binary frame " + i + " / " + nFrames)

        slice = im(:,:,i);

        CC = bwconncomp(slice);
        props = regionprops(CC, "Centroid");

        shiftedSlice = false(h, w);
        shiftedLabel = zeros(h, w, "uint32");
        shiftedCentroids = nan(CC.NumObjects, 2);

        for j = 1:CC.NumObjects

            pix = CC.PixelIdxList{j};
            [yy, xx] = ind2sub([h w], pix);

            xxShift = double(xx);
            yyShift = double(yy);

            c = props(j).Centroid;
            cxShift = c(1);
            cyShift = c(2);

            if i > anchorFrame

                % Shift frame i backward:
                % i -> i-1 -> ... -> anchorFrame
                for t = i:-1:(anchorFrame + 1)

                    if isempty(dxMaps{t}) || isempty(dyMaps{t}) || ~useVectorShift(t)
                        continue
                    end

                    [xxShift, yyShift, cxShift, cyShift] = applyVectorMapToObject( ...
                        xxShift, yyShift, cxShift, cyShift, ...
                        dxMaps{t}, dyMaps{t}, h, w, -1);
                end

            elseif i < anchorFrame

                % Shift frame i forward:
                % i -> i+1 -> ... -> anchorFrame
                for t = (i + 1):anchorFrame

                    if isempty(dxMaps{t}) || isempty(dyMaps{t}) || ~useVectorShift(t)
                        continue
                    end

                    [xxShift, yyShift, cxShift, cyShift] = applyVectorMapToObject( ...
                        xxShift, yyShift, cxShift, cyShift, ...
                        dxMaps{t}, dyMaps{t}, h, w, +1);
                end

            else

                % Anchor frame is unchanged.
            end

            xxShift = round(xxShift);
            yyShift = round(yyShift);

            valid = xxShift >= 1 & xxShift <= w & ...
                    yyShift >= 1 & yyShift <= h;

            xxShift = xxShift(valid);
            yyShift = yyShift(valid);

            if isempty(xxShift)
                continue
            end

            shiftedPix = sub2ind([h w], yyShift, xxShift);

            shiftedSlice(shiftedPix) = true;

            % Label value is the ORIGINAL object number in frame i.
            shiftedLabel(shiftedPix) = uint32(j);

            shiftedCentroids(j,:) = [mean(xxShift), mean(yyShift)];
        end

        cumShiftedSlices(:,:,i) = shiftedSlice;
        cumShiftedLabels{i} = shiftedLabel;
        cumObjCentroids{i} = shiftedCentroids;

        shiftedOut = uint8(255 * shiftedSlice);

        if i == 1
            imwrite(shiftedOut, outShiftedBinary, "tif", "Compression", "none");
        else
            imwrite(shiftedOut, outShiftedBinary, "tif", ...
                "WriteMode", "append", ...
                "Compression", "none");
        end
    end

    %% Track again using cumulative-shifted frame coordinates
    % This final tracking pass compares frame i and frame i+1 after BOTH have
    % been cumulatively shifted into the anchor-frame coordinate system.
    %
    % Track IDs are still assigned to the original object numbers from each frame,
    % because cumShiftedLabels stores original object indices as labels.

    trackIDsShift = cell(nFrames, 1);
    nextTrackIDShift = 0;
    trackAgeShift = [];

    % merge/lineage info for final TrackLineage sheet
    mergedIntoTrackID = nan(0, 1);
    mergeFrame = nan(0, 1);
    mergeSourceObject = nan(0, 1);
    mergeTargetObject = nan(0, 1);

    nObj = length(objAreas{1});

    trackIDsShift{1} = zeros(nObj, 1);

    for j = 1:nObj
        nextTrackIDShift = nextTrackIDShift + 1;
        trackIDsShift{1}(j) = nextTrackIDShift;
        trackAgeShift(nextTrackIDShift,1) = 1;

        mergedIntoTrackID(nextTrackIDShift,1) = nan;
        mergeFrame(nextTrackIDShift,1) = nan;
        mergeSourceObject(nextTrackIDShift,1) = nan;
        mergeTargetObject(nextTrackIDShift,1) = nan;
    end

    for i = 1:nFrames-1

        disp("Cumulative shifted linking transition " + i + " / " + (nFrames-1) + ...
             "  (frame " + i + " -> " + (i+1) + ")")
        tic

        curSliceShift = cumShiftedSlices(:,:,i);
        nxtSliceShift = cumShiftedSlices(:,:,i+1);

        curLabelShift = cumShiftedLabels{i};
        nxtLabelShift = cumShiftedLabels{i+1};

        nCur = length(objAreas{i});
        nNxt = length(objAreas{i+1});

        idsCur = trackIDsShift{i};
        idsNext = zeros(nNxt, 1);

        curArea = objAreas{i};
        nxtArea = objAreas{i+1};

        curUsed = false(nCur, 1);

        twoFrameStack = cat(3, curSliceShift, nxtSliceShift);
        CC3D = bwconncomp(twoFrameStack, 26);

        for g = 1:CC3D.NumObjects

            pix3D = CC3D.PixelIdxList{g};
            [y, x, z] = ind2sub(size(twoFrameStack), pix3D);

            idxCurPix = sub2ind(size(curSliceShift), y(z == 1), x(z == 1));
            curObjs = unique(curLabelShift(idxCurPix));
            curObjs(curObjs == 0) = [];
            curObjs = double(curObjs);

            idxNxtPix = sub2ind(size(nxtSliceShift), y(z == 2), x(z == 2));
            nxtObjs = unique(nxtLabelShift(idxNxtPix));
            nxtObjs(nxtObjs == 0) = [];
            nxtObjs = double(nxtObjs);

            if isempty(curObjs) || isempty(nxtObjs)
                continue
            end

            [~, curOrder] = sort(curArea(curObjs), "descend");
            [~, nxtOrder] = sort(nxtArea(nxtObjs), "descend");

            curObjs = curObjs(curOrder);
            nxtObjs = nxtObjs(nxtOrder);

            nPairs = min(length(curObjs), length(nxtObjs));

            for p = 1:nPairs

                cObj = curObjs(p);
                nObj = nxtObjs(p);

                if curUsed(cObj) || idsNext(nObj) ~= 0
                    continue
                end

                pairStack = false(size(twoFrameStack));
                pairStack(:,:,1) = curLabelShift == cObj;
                pairStack(:,:,2) = nxtLabelShift == nObj;

                pairCC = bwconncomp(pairStack, 26);
                same3DObject = pairCC.NumObjects == 1;

                areaOK = nxtArea(nObj) >= areaRatioMin * curArea(cObj);

                if same3DObject && areaOK

                    tid = idsCur(cObj);

                    idsNext(nObj) = tid;
                    curUsed(cObj) = true;

                    trackAgeShift(tid) = trackAgeShift(tid) + 1;
                end
            end

            %% Record merge-like events in this connected group
            % If multiple current objects touch/connect to a next object, but only one
            % current object continues, mark the other current tracks as merged.
            %
            % This is now evaluated in cumulative-shifted coordinates.

            linkedNxtObjs = nxtObjs(idsNext(nxtObjs) ~= 0);
            unlinkedCurObjs = curObjs(~curUsed(curObjs));

            if ~isempty(linkedNxtObjs) && ~isempty(unlinkedCurObjs)

                for m = 1:length(unlinkedCurObjs)

                    srcObj = unlinkedCurObjs(m);
                    srcTid = idsCur(srcObj);

                    % only record first merge/end call for this source track
                    if ~isnan(mergeFrame(srcTid))
                        continue
                    end

                    % Choose target linked-next object.
                    % If more than one linked next object is in the group, use nearest
                    % next centroid to the source current centroid in cumulative-shifted coordinates.

                    srcCent = cumObjCentroids{i}(srcObj,:);

                    if any(isnan(srcCent))
                        srcCent = objCentroids{i}(srcObj,:);
                    end

                    dTarget = zeros(length(linkedNxtObjs), 1);

                    for q = 1:length(linkedNxtObjs)

                        tgtObjTest = linkedNxtObjs(q);
                        tgtCent = cumObjCentroids{i+1}(tgtObjTest,:);

                        if any(isnan(tgtCent))
                            tgtCent = objCentroids{i+1}(tgtObjTest,:);
                        end

                        dTarget(q) = sqrt((tgtCent(1) - srcCent(1)).^2 + ...
                                          (tgtCent(2) - srcCent(2)).^2);
                    end

                    [~, bestTargetIdx] = min(dTarget);

                    tgtObj = linkedNxtObjs(bestTargetIdx);
                    tgtTid = idsNext(tgtObj);

                    % Confirm source current object directly connects to target next object
                    % in cumulative-shifted coordinates.
                    mergePairStack = false(size(twoFrameStack));
                    mergePairStack(:,:,1) = curLabelShift == srcObj;
                    mergePairStack(:,:,2) = nxtLabelShift == tgtObj;

                    mergePairCC = bwconncomp(mergePairStack, 26);
                    directMergeConnection = mergePairCC.NumObjects == 1;

                    if directMergeConnection
                        mergedIntoTrackID(srcTid,1) = tgtTid;
                        mergeFrame(srcTid,1) = i + 1;
                        mergeSourceObject(srcTid,1) = srcObj;
                        mergeTargetObject(srcTid,1) = tgtObj;
                    end
                end
            end
        end

        for j = 1:nNxt
            if idsNext(j) == 0
                nextTrackIDShift = nextTrackIDShift + 1;
                idsNext(j) = nextTrackIDShift;
                trackAgeShift(nextTrackIDShift,1) = 1;

                mergedIntoTrackID(nextTrackIDShift,1) = nan;
                mergeFrame(nextTrackIDShift,1) = nan;
                mergeSourceObject(nextTrackIDShift,1) = nan;
                mergeTargetObject(nextTrackIDShift,1) = nan;
            end
        end

        trackIDsShift{i+1} = idsNext;

        toc
    end

    nTracksShift = nextTrackIDShift;

    %% Export calibrated object table using per-frame full 3D stacks
    % For each frame, reloads that timepoint's full 3D binary stack (the
    % Surf_*_ZAl file matched to this tile earlier), finds its 3D
    % connected objects, and matches each final tracked 2D (MAX-
    % projection) object to the nearest 3D object centroid in XY (within
    % maxXYMatchDist). A successful match supplies the real calibrated
    % Volume and Z-centroid for that object/frame; 2D objects with no
    % close-enough 3D match are left with NaN volume/Z (e.g. if the object
    % dropped out of the 3D segmentation that frame).

    trackDuration_frames = zeros(nTracksShift, 1);

    for i = 1:nFrames

        ids = trackIDsShift{i};

        for j = 1:length(ids)
            tid = ids(j);
            trackDuration_frames(tid) = trackDuration_frames(tid) + 1;
        end
    end

    totalObjects = 0;

    for i = 1:nFrames
        totalObjects = totalObjects + length(trackIDsShift{i});
    end

    Tile = repmat(tileNum, totalObjects, 1);
    Object = zeros(totalObjects, 1);
    Frame = zeros(totalObjects, 1);
    SourceSurfNumber = zeros(totalObjects, 1);
    Source3DFile = strings(totalObjects, 1);

    Area_pixels = zeros(totalObjects, 1);
    Area_um2 = zeros(totalObjects, 1);

    Volume_voxels = nan(totalObjects, 1);
    Volume_um3 = nan(totalObjects, 1);

    CentroidX_pixels = nan(totalObjects, 1);
    CentroidY_pixels = nan(totalObjects, 1);
    CentroidZ_slices = nan(totalObjects, 1);

    CentroidX_um = nan(totalObjects, 1);
    CentroidY_um = nan(totalObjects, 1);
    CentroidZ_um = nan(totalObjects, 1);

    TrackID = zeros(totalObjects, 1);
    TrackDuration_frames = zeros(totalObjects, 1);

    Matched3DObject = nan(totalObjects, 1);
    MatchDistanceXY_pixels = nan(totalObjects, 1);

    row = 0;

    for i = 1:nFrames

        disp("Matching 2D objects to 3D objects for frame " + i + " / " + nFrames)
        tic

        this3DFile = tileSurfFullPaths(i);

        vol3D = tiffreadVolume(this3DFile) > 0;

        [vh, vw, ~] = size(vol3D);

        if vh ~= h || vw ~= w
            error("XY size mismatch for " + tileName + " frame " + i + ...
                ". MAX size is " + h + "x" + w + ...
                ", 3D stack size is " + vh + "x" + vw)
        end

        CC3Dfull = bwconncomp(vol3D, 26);
        props3D = regionprops3(CC3Dfull, "Volume", "Centroid");

        if CC3Dfull.NumObjects > 0
            cent3D = props3D.Centroid;
            volVals = props3D.Volume;
        else
            cent3D = zeros(0,3);
            volVals = zeros(0,1);
        end

        ids2D = trackIDsShift{i};
        xy2D = objCentroids{i};
        area2D = objAreas{i};

        for j = 1:length(ids2D)

            row = row + 1;

            tid = ids2D(j);

            Object(row) = j;
            Frame(row) = i;
            SourceSurfNumber(row) = tileSurfNums(i);
            Source3DFile(row) = tileSurfNames(i);

            Area_pixels(row) = area2D(j);
            Area_um2(row) = area2D(j) * pixelArea_um2;

            CentroidX_pixels(row) = xy2D(j,1);
            CentroidY_pixels(row) = xy2D(j,2);

            CentroidX_um(row) = xy2D(j,1) * pixelSizeXY;
            CentroidY_um(row) = xy2D(j,2) * pixelSizeXY;

            TrackID(row) = tid;
            TrackDuration_frames(row) = trackDuration_frames(tid);

            if ~isempty(cent3D)

                dXY = sqrt((cent3D(:,1) - xy2D(j,1)).^2 + ...
                           (cent3D(:,2) - xy2D(j,2)).^2);

                [dMin, matchIdx] = min(dXY);

                if dMin <= maxXYMatchDist

                    Matched3DObject(row) = matchIdx;
                    MatchDistanceXY_pixels(row) = dMin;

                    Volume_voxels(row) = volVals(matchIdx);
                    Volume_um3(row) = volVals(matchIdx) * voxelVolume_um3;

                    CentroidZ_slices(row) = cent3D(matchIdx,3);
                    CentroidZ_um(row) = (cent3D(matchIdx,3) - 1) * zStep;
                end
            end
        end

        toc
    end

    objectTable = table( ...
        Tile, ...
        Object, ...
        Frame, ...
        SourceSurfNumber, ...
        Source3DFile, ...
        Area_pixels, ...
        Area_um2, ...
        Volume_voxels, ...
        Volume_um3, ...
        CentroidX_pixels, ...
        CentroidY_pixels, ...
        CentroidZ_slices, ...
        CentroidX_um, ...
        CentroidY_um, ...
        CentroidZ_um, ...
        TrackID, ...
        TrackDuration_frames, ...
        Matched3DObject, ...
        MatchDistanceXY_pixels);

    %% Build track lineage table
    % One row per final (post-cumulative-shift) track: its start/end frame,
    % total duration, and how it ended -- merged into another track
    % (EndReason = "Merged", with the merge target/frame/objects recorded
    % from the merge-detection step above), still present when the movie
    % ended ("StillPresentAtEnd"), or simply lost/disappeared.

    LineageTile = repmat(tileNum, nTracksShift, 1);
    LineageTrackID = (1:nTracksShift)';

    StartFrame = nan(nTracksShift, 1);
    EndFrame = nan(nTracksShift, 1);
    Duration_frames = nan(nTracksShift, 1);

    StartSourceSurfNumber = nan(nTracksShift, 1);
    EndSourceSurfNumber = nan(nTracksShift, 1);
    MergeSourceSurfNumber = nan(nTracksShift, 1);

    StartObject = nan(nTracksShift, 1);
    EndObject = nan(nTracksShift, 1);

    EndReason = strings(nTracksShift, 1);

    MergedIntoTrackID = mergedIntoTrackID(1:nTracksShift);
    MergeFrame = mergeFrame(1:nTracksShift);
    MergeSourceObject = mergeSourceObject(1:nTracksShift);
    MergeTargetObject = mergeTargetObject(1:nTracksShift);

    for tid = 1:nTracksShift

        framesFound = [];
        objectsFound = [];

        for i = 1:nFrames

            objIdx = find(trackIDsShift{i} == tid, 1);

            if ~isempty(objIdx)
                framesFound(end+1,1) = i;
                objectsFound(end+1,1) = objIdx;
            end
        end

        if isempty(framesFound)
            continue
        end

        StartFrame(tid) = framesFound(1);
        EndFrame(tid) = framesFound(end);
        Duration_frames(tid) = length(framesFound);

        StartObject(tid) = objectsFound(1);
        EndObject(tid) = objectsFound(end);

        StartSourceSurfNumber(tid) = tileSurfNums(StartFrame(tid));
        EndSourceSurfNumber(tid) = tileSurfNums(EndFrame(tid));

        if ~isnan(MergeFrame(tid))
            EndReason(tid) = "Merged";
            MergeSourceSurfNumber(tid) = tileSurfNums(MergeFrame(tid));
        elseif EndFrame(tid) == nFrames
            EndReason(tid) = "StillPresentAtEnd";
        else
            EndReason(tid) = "Disappeared";
        end
    end

    trackLineageTable = table( ...
        LineageTile, ...
        LineageTrackID, ...
        StartFrame, ...
        EndFrame, ...
        Duration_frames, ...
        StartSourceSurfNumber, ...
        EndSourceSurfNumber, ...
        StartObject, ...
        EndObject, ...
        EndReason, ...
        MergedIntoTrackID, ...
        MergeFrame, ...
        MergeSourceSurfNumber, ...
        MergeSourceObject, ...
        MergeTargetObject);

    %% Write exactly two Excel sheets

    writetable(objectTable, outExcel, "Sheet", "Objects")
    writetable(trackLineageTable, outExcel, "Sheet", "TrackLineage")

    disp("Saved Excel file with Objects and TrackLineage sheets:")
    disp(outExcel)

    %% Collect shifted-linking tracks using original centroid positions
    % Re-groups the final (cumulative-shift-derived) TrackIDs back onto
    % each object's ORIGINAL (unshifted) centroid position, purely so the
    % track overlay below can be drawn on the real, un-registered movie
    % rather than the anchor-registered coordinate space used internally
    % for linking.

    trackFramesShift = cell(nTracksShift, 1);
    trackXYShift = cell(nTracksShift, 1);

    for i = 1:nFrames

        ids = trackIDsShift{i};
        xy = objCentroids{i};

        for j = 1:length(ids)

            tid = ids(j);

            trackFramesShift{tid}(end+1,1) = i;
            trackXYShift{tid}(end+1,:) = xy(j,:);
        end
    end

    %% Write shifted-linking TrackMate-style overlay on original movie

    colors = uint8(255 * hsv(max(nTracksShift,1)));

    % Precompute every line segment and dot once, tagged with the frame
    % at which it appears. A segment appears at the frame of its LATER
    % endpoint (same rule as the original fr <= i logic).

    nPts = cellfun(@(v) size(v,1), trackFramesShift);

    totalSegs = sum(max(nPts - 1, 0));
    totalDots = sum(nPts);

    segXY = zeros(totalSegs, 4);
    segFrame = zeros(totalSegs, 1);
    segColor = zeros(totalSegs, 3);

    dotXY = zeros(totalDots, 2);
    dotFrame = zeros(totalDots, 1);
    dotColor = zeros(totalDots, 3);

    segRow = 0;
    dotRow = 0;

    for t = 1:nTracksShift

        fr = trackFramesShift{t};
        xy = trackXYShift{t};

        color = double(colors(t,:));
        n = nPts(t);

        if n >= 2
            idx = segRow + (1:n-1);
            segXY(idx,:) = [xy(1:n-1,:), xy(2:n,:)];
            segFrame(idx) = fr(2:n);
            segColor(idx,:) = repmat(color, n-1, 1);
            segRow = segRow + n - 1;
        end

        idx = dotRow + (1:n);
        dotXY(idx,:) = xy;
        dotFrame(idx) = fr;
        dotColor(idx,:) = repmat(color, n, 1);
        dotRow = dotRow + n;
    end

    for i = 1:nFrames

        rgb = uint8(255 * repmat(im(:,:,i), 1, 1, 3));

        sel = segFrame <= i;

        if any(sel)
            rgb = insertShape(rgb, "Line", segXY(sel,:), ...
                "Color", segColor(sel,:), ...
                "LineWidth", lineWidth);
        end

        dsel = dotFrame == i;

        if any(dsel)
            rgb = insertShape(rgb, "FilledCircle", ...
                [dotXY(dsel,:), repmat(spotRadius, nnz(dsel), 1)], ...
                "Color", dotColor(dsel,:), ...
                "Opacity", 1);
        end

        if i == 1
            imwrite(rgb, outShiftedTracks, "tif", "Compression", "none");
        else
            imwrite(rgb, outShiftedTracks, "tif", "WriteMode", "append", "Compression", "none");
        end
    end

    disp("Done shifted linking for " + tileName)
    disp("Shifted-linking tracks:")
    disp(nTracksShift)
    disp("Saved:")
    disp(outShiftedTracks)

end

disp("All tiles done.")

function [xxShift, yyShift, cxShift, cyShift] = applyVectorMapToObject( ...
    xxShift, yyShift, cxShift, cyShift, dxMap, dyMap, h, w, directionSign)
% Shifts one object's pixel coordinates (and its running centroid) by the
% displacement field value AT ITS OWN CURRENT CENTROID -- i.e. treats the
% whole object as moving by one local vector rather than sampling the
% field per-pixel, since the field varies smoothly and slowly at the
% scale of a single object. directionSign flips the sign depending on
% whether this call is chaining the object forward or backward toward the
% anchor frame (see the two loops in the caller). The centroid (cxShift/
% cyShift) is updated alongside the pixels so the NEXT chained step in a
% multi-frame walk samples the field at the object's new, already-shifted
% position rather than its original one.

    mapH = size(dxMap, 1);
    mapW = size(dxMap, 2);

    cxNow = round(cxShift);
    cyNow = round(cyShift);

    cxNow = min(max(cxNow, 1), w);
    cyNow = min(max(cyNow, 1), h);

    cxMap = round((cxNow - 1) * (mapW - 1) / (w - 1) + 1);
    cyMap = round((cyNow - 1) * (mapH - 1) / (h - 1) + 1);

    cxMap = min(max(cxMap, 1), mapW);
    cyMap = min(max(cyMap, 1), mapH);

    dx = dxMap(cyMap, cxMap);
    dy = dyMap(cyMap, cxMap);

    xxShift = xxShift + directionSign * dx;
    yyShift = yyShift + directionSign * dy;

    cxShift = cxShift + directionSign * dx;
    cyShift = cyShift + directionSign * dy;
end
