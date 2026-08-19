%% ProcessSurfaceTiles_ZDriftStitchCorrection.m
% Author: Hasreet Gill, Bassler Lab, HHMI/Princeton University
% Developed/revised with Claude (Sonnet 5, Anthropic) - AI-assisted code review and comments.

% Per-experiment-folder, per-timepoint, per-spatial-tile, per-channel
% object segmentation and drift/artifact correction pipeline for the
% surface (Surf*) Z-stacks.
%
% Assumes each experiment folder already has:
%   - Surf_<n>_...488...tif (and matching ...561... and/or ...640... tifs)
%     Z-stack tiffs, one triplet/pair per imaging timepoint.
%   - Shift_Stack_TotalShifts.csv (written by the manual+auto XY drift
%     correction script) giving the whole-frame TotalX/TotalY shift to
%     apply at each timepoint.
%
% For every timepoint, the field of view is split into q spatial tiles
% (so a whole Z-stack doesn't have to be held in memory at once) and,
% per tile per channel-geometry job, the pipeline:
%   1. Loads the raw Z-stack for that tile, applying the timepoint's
%      global XY drift shift while loading (loadAndGlobalShiftTile).
%   2. Thresholds to a binary mask, removes small 2D noise, and drops
%      "swimming" objects that only exist for a couple of Z-slices --
%      unless they're large, in which case they're rescued and
%      reconnected to their full extent via proximity dilation
%      (removeSwimmingCells / proximityLinkRescuedObjects).
%   3. Detects two distinct kinds of Z-stack artifacts and computes a
%      per-object correction for each:
%        - Stitching oscillation: an object's XY centroid jumps then
%          snaps back between adjacent Z-slices, along one axis --
%          characteristic of a tile-stitching wobble, not real 3D shape
%          (identifyStitchOscillations).
%        - Bead-based Z drift: small, circular, many-slice objects are
%          treated as fiducial beads, and the median bead-to-bead
%          centroid displacement between consecutive Z-slices gives a
%          smooth cumulative XY-drift-vs-Z curve (identifyZDrift).
%   4. Applies both corrections together, per object per Z-slice, moving
%      each object's (hole-filled) footprint and backfilling the vacated
%      area with the local background intensity (applyObjectCorrections).
%   5. Runs a second cleanup pass dropping objects that are still
%      unstable after correction -- either their projected XY footprint
%      is way bigger than any single slice (i.e. they wander) or their
%      per-slice centroid strays far from their own median position
%      (removeShakyObjects).
%   6. Aligns every surviving object's top Z-slice to the top of the
%      stack, so all objects share a common Z reference (presumably the
%      coverslip/surface) regardless of where the focus happened to sit
%      that timepoint (alignObjectsToTopZ).
%   7. Writes the final binary mask, the final masked raw-intensity
%      stack per channel, a table of per-object volume/shape/intensity
%      measurements, and appends this timepoint's max-Z-projection onto
%      a running per-tile time-lapse TIFF (with ImageJ hyperstack
%      metadata so it opens as a proper T-series).
%
% Channel handling per timepoint:
%   - If a 561 image exists alongside 488, they are treated as two
%     INDEPENDENT geometry channels (each segmented and corrected using
%     its own mask).
%   - Otherwise, if a 640 image exists alongside 488, 488 supplies the
%     geometry/mask and 640 is just measured through that same mask
%     (paired mode) -- 640 has no mask of its own here.
%   - If both 561 and 640 exist for the same 488 image, that's treated
%     as an unexpected/ambiguous case: it's forced into 488/561
%     independent mode and 640 is ignored (with a warning), rather than
%     silently guessing which pairing was intended.
%
% Output filenames chain together suffixes describing which corrections
% have been applied so far (…_ShBin_NoNoiseSwim_Resc_StDrCorr_NoShaky_ZAl…),
% making it possible to inspect any intermediate stage of the pipeline.

clear; clc; close all

path = uigetdir; % High resolution imaging
cd(path);

folders = dir("2026*"); % folders with .tifs, XY correction computed from Shift_ files

for f = 1:length(folders)

    cd(fullfile(path, folders(f).name))

    % 488 is always present and drives the per-timepoint file list / count;
    % 561 and 640 filenames are derived from it further down.
    files488 = dir("Surf*488*.tif");
    names = {files488.name};
    nums = cellfun(@(x) str2double(regexp(x, '(?<=Surf_)\d+', 'match', 'once')), names);
    [~, idx] = sort(nums);
    files488 = files488(idx);
    info = imfinfo(files488(1).name);
    w = info(1).Width;
    h = info(1).Height;

    % Global (whole-frame, per-timepoint) XY drift correction computed by
    % the earlier bead-tracking script -- applied while loading each tile
    % below, before any of the per-object corrections in this script.
    shifts = readtable("Shift_Stack_TotalShifts.csv");
    q = 4;
    split = makeSplit(w, h, q);
    % Padding for the region read around each tile so imtranslate has real
    % image data to shift in from near the tile edges instead of zero-fill;
    % must be at least as large as the biggest global shift seen this
    % experiment, plus a small margin.
    pad = ceil(max(abs([shifts.TotalX; shifts.TotalY]))) + 5;

    count = 0; % running timepoint counter, used as the frame index when appending to the time-lapse max-projection tif

    for i = 1:length(files488)

        count = count + 1;

        name488 = files488(i).name;
        name561 = char(strrep(string(name488), "488", "561"));
        name640 = char(strrep(string(name488), "488", "640"));
        has561 = isfile(name561);
        has640 = isfile(name640);

        dx = shifts.TotalX(i);
        dy = shifts.TotalY(i);

        fprintf('\n==============================\n');
        fprintf('Folder: %s\n', folders(f).name);
        fprintf('Image: %s\n', name488);
        fprintf('Global time shift: X = %.3f, Y = %.3f\n', dx, dy);

        if has561 && has640
            warning('Both 561 and 640 found for %s. Processing as 488/561 independent mode and ignoring 640.', name488);
            has640 = false;
        end

        for j = 1:q

            % channelJobs enumerates the geometry-channel(s) to segment for
            % this timepoint: 488+561 run as two independent masks, while
            % 488+640 runs as a single job where 640 is only ever measured
            % through the 488-derived mask (see channel-handling note above).
            if has561
                channelJobs = {
                    struct("geomName", name488, "geomChan", "488", "rawNames", {{name488}}, "rawChans", {{"488"}})
                    struct("geomName", name561, "geomChan", "561", "rawNames", {{name561}}, "rawChans", {{"561"}})
                    };
            elseif has640
                channelJobs = {
                    struct("geomName", name488, "geomChan", "488", "rawNames", {{name488, name640}}, "rawChans", {{"488", "640"}})
                    };
            end

            for jobIdx = 1:numel(channelJobs)

                job = channelJobs{jobIdx};
                fprintf('\nProcessing tile %d, geometry channel %s: %s\n', ...
                    j, job.geomChan, job.geomName);
                rawTileByChannel = struct();

                for rr = 1:numel(job.rawNames)
                    rawName = job.rawNames{rr};
                    rawChan = job.rawChans{rr};
                    stepName = sprintf('Load + global shift raw channel %s: %s', rawChan, rawName);
                    tStep = startTimer(stepName);
                    rawTile = loadAndGlobalShiftTile(rawName, split, j, pad, h, w, dx, dy);
                    stopTimer(stepName, tStep);

                    [~, ~, nZ] = size(rawTile);
                    rawTileByChannel.("ch" + rawChan) = rawTile;
                    clear rawTile
                end

                geomField = "ch" + job.geomChan;
                rawGeomTile = rawTileByChannel.(geomField);

                [Y, X, Z] = size(rawGeomTile);
                nZ = Z;

               %% make binary
                % Threshold each Z-slice independently (channel-dependent
                % intensity threshold) and drop small 2D speckle noise.
                stepName = 'Make shifted binary and remove swimmers';
                tStep = startTimer(stepName);

                minObjectSlices = 2; %%%%%%%%%%%%%
                rescueVolume = 500; %%%%%%%%%%%%%
                if job.geomChan == "561"
                    binThresh = 115; %%%%%%%%%%%%%
                else
                    binThresh = 140; %%%%%%%%%%%%%
                end
                noise2D = 5;
                noise3D = 40; %%%%%%%%%%%%%

                bin = false(size(rawGeomTile));
                for k = 1:nZ
                    slice = rawGeomTile(:,:,k);
                    binslice = slice >= binThresh;
                    binslice = bwareaopen(binslice, noise2D);
                    bin(:,:,k) = binslice;
                end
                % bin = bwareaopen(bin, noise3D);

                % remove likely swimming cells / short-Z objects
                [binTile, rescuedMask, CC3D] = removeSwimmingCells(bin, minObjectSlices, rescueVolume);
                stopTimer(stepName, tStep);

                % proximity-link rescued objects across adjacent slices
                stepName = 'Proximity-link rescued objects';
                tStep = startTimer(stepName);

                linkDist = 20; %%%%%%%%%%%%%
                if any(rescuedMask(:))
                    [binTile, CC3D] = proximityLinkRescuedObjects(bin, binTile, rescuedMask, minObjectSlices, linkDist);
                end

                outBin = [extractBefore(job.geomName, '1X'), num2str(j), '_ShBin_NoNoiseSwim_Resc.tif'];
                writeTiffStack(binTile, outBin);
                L3D = labelmatrix(CC3D);
                stopTimer(stepName, tStep);

                %% stitching oscillation correction
                % Flags objects whose XY centroid jumps then snaps back
                % between adjacent Z-slices along a single axis (a
                % tile-stitching wobble, not real object motion) and
                % records the per-object per-Z correction needed to pull
                % each affected slice back toward that object's own
                % median centroid.
                stepName = 'Detect stitching oscillations';
                tStep = startTimer(stepName);

                jumpThresh = 15; %%%%%%%%%%%%%
                returnThresh = 8; %%%%%%%%%%%%%
                axisRatioThresh = 0.30; %%%%%%%%%%%%%

                [stitchObjID, stitchDx, stitchDy, oscTable] = identifyStitchOscillations(CC3D, [Y X Z], nZ, jumpThresh, returnThresh, axisRatioThresh);
                oscComponentID = oscTable.ComponentID;
                fprintf('Tile %d: found %d stitching oscillating objects \n', j, numel(unique(oscComponentID)));
                outOscCSV = [extractBefore(job.geomName, '1X'), num2str(j), '_StitchOscillationCorrections.csv'];
                if isfile(outOscCSV)
                    delete(outOscCSV);
                end
                writetable(oscTable, outOscCSV);
                stopTimer(stepName, tStep);

                %% find beads and calculate smooth XY drift in Z
                % Fiducial beads are identified as small, highly circular
                % objects present across many consecutive Z-slices, then
                % the median XY displacement of matched beads between
                % adjacent slices is cumulatively summed to give a smooth
                % XY-drift-vs-Z curve for this tile/channel.
                stepName = 'Identify bead-based Z drift and save shift CSV';
                tStep = startTimer(stepName);

                maxBeadArea = 60; %100; %%%%%%%%%%%%%
                minBeadSlices = 5; %%%%%%%%%%%%%
                minMedianCircularity = 0.8;

                [zDx, zDy, zStepDx, zStepDy, nBeads, beadCount] = identifyZDrift(CC3D, L3D, [Y X Z], ...
                    nZ, maxBeadArea, minBeadSlices, minMedianCircularity);
                if beadCount == 0
                    % No beads found with the primary threshold/size --
                    % rebinarize from scratch with a looser max area and
                    % retry once before giving up on Z-drift correction
                    % for this tile.
                    fprintf('No bead candidates found. Retrying.\n');
                    binThreshredo = 115;
                    maxBeadArea = 200; %%%%%%%%%%%%%
                    minBeadSlices = 5; %%%%%%%%%%%%%
                    minMedianCircularity = 0.8;

                    binredo = false(size(rawGeomTile));
                    for k = 1:nZ
                        sliceredo = rawGeomTile(:,:,k);
                        binsliceredo = sliceredo >= binThreshredo;
                        binsliceredo = bwareaopen(binsliceredo, noise2D);
                        binredo(:,:,k) = binsliceredo;
                    end
                    % binredo = bwareaopen(binredo, noise3D);

                    CC3D_redo = bwconncomp(binredo, 26);
                    L3D_redo = labelmatrix(CC3D_redo);
                    [zDx, zDy, zStepDx, zStepDy, nBeads, beadCount] = identifyZDrift(CC3D_redo, L3D_redo, [Y X Z], ...
                        nZ, maxBeadArea, minBeadSlices, minMedianCircularity);
                    clear binredo CC3D_redo L3D_redo
                end

                fprintf('Tile %d channel %s: found %d bead candidates\n', ...
                    j, job.geomChan, beadCount);
                shiftTable = table((1:nZ)', zDx, zDy, zStepDx, zStepDy, nBeads, ...
                    'VariableNames', {'Z','DriftX','DriftY','StepX','StepY','NumBeads'});
                outShiftCSV = [extractBefore(job.geomName, '1X'), ...
                    num2str(j), '_ZDriftShifts.csv'];
                if isfile(outShiftCSV)
                    delete(outShiftCSV);
                end
                writetable(shiftTable, outShiftCSV);
                stopTimer(stepName, tStep);

               %% XY correction by object in Z
                % Apply the bead-based Z-drift shift and any stitching
                % correction together, per object per Z-slice (see
                % applyObjectCorrections), rather than shifting whole
                % slices -- different objects in the same slice can need
                % different stitching corrections.
                stepName = 'Apply object corrections';
                tStep = startTimer(stepName);
                rawFields = cellstr(strcat("ch", string(job.rawChans)));
                outBinCorr = [extractBefore(job.geomName, '1X'), num2str(j), '_ShBin_NoNoiseSwim_Resc_StDrCorr.tif'];

                [binTileCorr, rawTileCorrByChannel] = applyObjectCorrections(binTile, L3D, rawTileByChannel, rawFields, zDx, zDy, stitchObjID, stitchDx, stitchDy);
                writeTiffStack(binTileCorr, outBinCorr);

                stopTimer(stepName, tStep);
                clear binTile rawTileByChannel rawGeomTile L3D CC3D rescuedMask

                %% remove shaky objects after Z correction, before top-Z alignment
                % Second QC pass: drop any remaining short-lived fragments,
                % then remove objects that are still unstable even after
                % correction (projected footprint much bigger than any
                % single slice, or centroid wandering far from its own
                % median position) -- see removeShakyObjects.
                stepName = 'Remove shaky objects after Z correction';
                tStep = startTimer(stepName);
                minObjectSlicesPost = 3;
                maxProjRatioThresh = 2.5; %%%%%%%%%%%%%
                areaThreshold = 2000; %%%%%%%%%%%%%
                smallArea = 100;
                shakyThresh = 20; %%%%%%%%%%%%%

                [binTileCorr, ~, ~] = removeSwimmingCells(binTileCorr, minObjectSlicesPost, []);
                binTileCorr = removeShakyObjects(binTileCorr, maxProjRatioThresh, shakyThresh, smallArea, areaThreshold);

                outBinClean = [extractBefore(job.geomName, '1X'), num2str(j), '_ShBin_NoNoiseSwim_Resc_StDrCorr_NoShaky.tif'];
                writeTiffStack(binTileCorr, outBinClean);

                stopTimer(stepName, tStep);

                %% align final kept objects to top Z and save final binary
                % Shift every surviving 3D object bodily along Z so its
                % topmost slice sits at the top of the stack -- normalizes
                % away timepoint-to-timepoint focus/Z-stage variation by
                % anchoring all objects to a shared reference plane
                % (presumably the coverslip/surface).
                stepName = 'Align final objects to top Z and write final binary';
                tStep = startTimer(stepName);

                outBinFinal = [extractBefore(job.geomName, '1X'), ...
                    num2str(j), '_ShBin_NoNoiseSwim_Resc_StDrCorr_NoShaky_ZAl.tif'];

                [binTileFinal, rawTileFinalByChannel] = alignObjectsToTopZ( binTileCorr, rawTileCorrByChannel, rawFields);

                clear binTileCorr rawTileCorrByChannel
                writeTiffStack(binTileFinal, outBinFinal);

                stopTimer(stepName, tStep);

                %% append final max projection over time
                % Builds a running time-lapse TIFF (one per tile/geometry
                % channel, spanning every timepoint processed so far) of
                % this frame's max-Z-projection, so drift/segmentation
                % quality can be checked over time without opening every
                % individual Z-stack.
                stepName = 'Append final max projection time frame';
                tStep = startTimer(stepName);

                maxFinal = max(binTileFinal,[],3);
                outMaxFinal = ['Geom' char(job.geomChan),'_Tile' num2str(j), '_ShBin_NoNoiseSwim_Resc_StDrCorr_NoShaky_ZAl_MAX.tif'];
                if count == 1
                    % Fresh run: don't append this timepoint's frame onto
                    % a leftover file from a previous run.
                    if isfile(outMaxFinal)
                        delete(outMaxFinal);
                    end
                end
                sliceOut = uint8(maxFinal)*255;
                if isfile(outMaxFinal)
                    imwrite(sliceOut, outMaxFinal, 'WriteMode','append','Compression','none');
                else
                    imwrite(sliceOut, outMaxFinal, 'Compression','none');
                end

                frameNum = count;
                addImageJTimeMetadata(outMaxFinal, frameNum);

                stopTimer(stepName, tStep);

                %% final object properties
                % Per-final-object geometry (volume, principal axis
                % length, centroid) plus, for every raw channel measured
                % in this job, the mean/std intensity over that object's
                % own voxels -- written one row per object, per tile per
                % timepoint.
                stepName = 'Measure and save final object properties';
                tStep = startTimer(stepName);

                CCtmp = bwconncomp(binTileFinal,26);

                props = regionprops3(CCtmp, ...
                    'Volume','PrincipalAxisLength','Centroid');

                for rf = 1:numel(rawFields)

                    fieldName = rawFields{rf};
                    chanName = erase(fieldName,"ch");

                    rawStack = rawTileFinalByChannel.(fieldName);

                    meanVals = zeros(CCtmp.NumObjects,1);
                    stdVals  = zeros(CCtmp.NumObjects,1);

                    for obj = 1:CCtmp.NumObjects

                        vox = CCtmp.PixelIdxList{obj};
                        vals = double(rawStack(vox));

                        meanVals(obj) = mean(vals);
                        stdVals(obj)  = std(vals);

                    end

                    props.("MeanIntensity_" + chanName) = meanVals;
                    props.("StdIntensity_" + chanName)  = stdVals;

                end

                outPropsFinal = [extractBefore(job.geomName,'1X'), ...
                    num2str(j), '_ShBin_NoNoiseSwim_Resc_StDrCorr_NoShaky_ZAl.csv'];

                writetable(props, outPropsFinal);

                stopTimer(stepName, tStep);

                clear CCtmp props meanVals stdVals vals vox

                %% write masked final raw/color channel stacks
                % The actual processed intensity images (not just the
                % binary mask) for every raw channel in this job, fully
                % corrected and Z-aligned, background-filled outside
                % objects -- one TIFF per raw channel.
                stepName = 'Write final masked raw/color output stacks';
                tStep = startTimer(stepName);

                for rf = 1:numel(rawFields)

                    fieldName = rawFields{rf};
                    chanName = erase(fieldName, "ch");

                    outRawMasked = [extractBefore(job.rawNames{rf}, '1X'), ...
                        num2str(j), '_Raw' char(chanName) ...
                        '_ShBin_NoNoiseSwim_Resc_StDrCorr_NoShaky_ZAl_Mskd.tif'];

                    rawMasked = rawTileFinalByChannel.(fieldName);

                    writeTiffStack(rawMasked, outRawMasked);

                    clear rawMasked

                end

                stopTimer(stepName, tStep);

                clear binTileFinal rawTileFinalByChannel maxFinal sliceOut
            end
        end
    end

    cd(path)

end

function rawTile = loadAndGlobalShiftTile(rawName, split, j, pad, h, w, dx, dy)
% Loads only the padded pixel region for tile j (via tiffreadVolume's
% PixelRegion, so the whole multi-GB Z-stack doesn't need to be read into
% memory), applies the timepoint's global XY shift to every Z-slice, then
% crops back down to exactly the tile's own footprint. The padding exists
% so imtranslate has real neighboring image data to pull in near the
% tile's edges instead of filling with zeros right where a real object
% might be shifted into view.

    rows0 = split{j,1};
    cols0 = split{j,2};

    r1 = rows0(1);
    r2 = rows0(3);
    c1 = cols0(1);
    c2 = cols0(3);

    rp1 = max(1, r1 - pad);
    rp2 = min(h, r2 + pad);
    cp1 = max(1, c1 - pad);
    cp2 = min(w, c2 + pad);

    paddedRegion = {[rp1 1 rp2], [cp1 1 cp2], [1 Inf]};

    cropRows = (r1 - rp1 + 1):(r2 - rp1 + 1);
    cropCols = (c1 - cp1 + 1):(c2 - cp1 + 1);

    imPad = tiffreadVolume(rawName, 'PixelRegion', paddedRegion);
    [~, ~, nZ] = size(imPad);

    rawTile = zeros(numel(cropRows), numel(cropCols), nZ, 'like', imPad);

    for k = 1:nZ
        slicePad = imPad(:,:,k);
        sliceShiftPad = imtranslate(slicePad, [dx dy], ...
            'OutputView', 'same', 'FillValues', 0);
        rawTile(:,:,k) = sliceShiftPad(cropRows, cropCols);
    end

end

function [binTile, rescuedMask, CC3D] = removeSwimmingCells(binTile, minObjectSlices, rescueVolume)
% Drops 3D objects that span fewer than minObjectSlices Z-slices (likely
% transient/motion-blurred "swimmers" rather than real static structures).
% If rescueVolume is given (non-empty) and such a short-lived object is
% still large by voxel count, it isn't deleted -- its voxels are instead
% marked in rescuedMask so a caller can try to recover its full extent
% via proximity linking (see proximityLinkRescuedObjects) instead of
% losing what might be a real, just poorly-thresholded, object.

    doRescue = nargin >= 3 && ~isempty(rescueVolume);
    rescuedMask = false(size(binTile));

    props = regionprops3(binTile, 'Volume', 'VoxelList', 'VoxelIdxList');
    for obj = 1:length(props.VoxelList)

        pix = props.VoxelIdxList{obj};
        spansFewSlices = numel(unique(props.VoxelList{obj}(:,3))) < minObjectSlices;

        if spansFewSlices
            if doRescue && numel(pix) > rescueVolume
                rescuedMask(pix) = true;
            else
                binTile(pix) = false;
            end
        end
    end

    CC3D = bwconncomp(binTile, 26);
end

function [binTile, CC3D] = proximityLinkRescuedObjects(bin, binTile, rescuedMask, minObjectSlices, linkDist)
% For each object rescued above (a short-Z-span but large fragment),
% checks whether it connects up to a bigger structure once the ORIGINAL
% (pre-noise-removal) binary is dilated by linkDist in XY on each slice --
% i.e. treats nearby fragments of the same underlying object, which
% thresholding/noise-removal split apart, as one linked object. Only
% objects that reach minObjectSlices once fully recovered are kept; the
% final binary and CC3D are rebuilt to include these recovered full
% objects as single components alongside the normally-kept ones.

    [Y, X, nZ] = size(binTile);

    binBase = binTile;
    binBase(rescuedMask) = false;

    % dilate rescued objects slice-by-slice in XY
    se = strel("disk", linkDist);
    binDil = false(size(binTile));
    for z = 1:nZ
        binDil(:,:,z) = imdilate(bin(:,:,z), se);
    end

    % label the dilated unfiltered binary
    CCdil = bwconncomp(binDil, 26);
    Ldil = labelmatrix(CCdil);

    % label rescued seed objects
    CCrescued = bwconncomp(rescuedMask, 26);
    matchedDilLabels = zeros(CCrescued.NumObjects, 1);
    for r = 1:CCrescued.NumObjects
        rescuedPix = CCrescued.PixelIdxList{r};

        % find which dilated unfiltered object this rescued seed overlaps
        dilLabels = Ldil(rescuedPix);
        dilLabels(dilLabels == 0) = [];
        if isempty(dilLabels)
            continue
        end

        % usually there should be one label; mode is safest if there are several
        matchedDilLabels(r) = mode(double(dilLabels));
    end

    matchedDilLabels = unique(matchedDilLabels);
    matchedDilLabels(matchedDilLabels == 0) = [];

    recoveredGroupPix = {};
    keepRecovered = false(size(binTile));
    for g = 1:numel(matchedDilLabels)

        dilID = matchedDilLabels(g);

        % pixels from the ORIGINAL undilated binary that fall inside
        % this matched dilated object
        recoveredPix = find((Ldil == dilID) & bin);

        if isempty(recoveredPix)
            continue
        end

        [~, ~, zz] = ind2sub([Y X nZ], recoveredPix);
        zSpan = numel(unique(zz));

        if zSpan < minObjectSlices
            continue
        end

        keepRecovered(recoveredPix) = true;
        recoveredGroupPix{end+1,1} = recoveredPix;

    end

% final binary = normal filtered objects + recovered full objects from bin
binTile = binBase | keepRecovered;

% build CC3D so each recovered proximity-linked object is treated as one object
CC3D = bwconncomp(binBase, 26);

for g = 1:numel(recoveredGroupPix)

    pixGroup = recoveredGroupPix{g};

    CC3D.NumObjects = CC3D.NumObjects + 1;
    CC3D.PixelIdxList{CC3D.NumObjects} = pixGroup;

end

end

function [stitchObjID, stitchDx, stitchDy, oscTable] = identifyStitchOscillations(CC3D, imSize, nZ, jumpThresh, returnThresh, axisRatioThresh)
% For every 3D object, walks its per-slice centroid across CONSECUTIVE
% Z-slices looking for a "jump then return" pattern: centroid moves by
% >= jumpThresh from one slice to the next, then moves back to within
% returnThresh of where it started two slices ago. That round-trip,
% restricted to slices with no Z gaps, is the signature of a
% tile-stitching seam wobble rather than genuine 3D object movement.
% axisRatioThresh additionally requires the jump/return to be
% predominantly along ONE axis (X or Y) -- a real stitching offset is
% axis-aligned, unlike organic 3D motion. Any object matching this
% pattern anywhere along its length is flagged as oscillating for its
% ENTIRE Z extent, and each of its slices gets a correction that pulls
% that slice's centroid back toward the object's own overall median
% position (its best estimate of where it "should" be).

    Y = imSize(1);
    X = imSize(2);
    Z = imSize(3);

    stitchObjID = cell(nZ,1);
    stitchDx = cell(nZ,1);
    stitchDy = cell(nZ,1);

    oscComponentID = [];
    oscZ = [];
    oscDx = [];
    oscDy = [];

    for obj = 1:CC3D.NumObjects

        pix = CC3D.PixelIdxList{obj};
        [yy, xx, zz] = ind2sub([Y X Z], pix);
        zList = unique(zz);

        if numel(zList) < 3
            continue
        end

        cx = zeros(numel(zList),1);
        cy = zeros(numel(zList),1);

        for zi = 1:numel(zList)

            ztmp = zList(zi);
            idxZ = zz == ztmp;

            cx(zi) = mean(xx(idxZ));
            cy(zi) = mean(yy(idxZ));

        end

        isOscObj = false;

        for zi = 2:numel(zList)-1

            zPrev = zList(zi-1);
            zCurr = zList(zi);
            zNext = zList(zi+1);

            if zCurr ~= zPrev + 1 || zNext ~= zCurr + 1
                continue
            end

            pPrev = [cx(zi-1), cy(zi-1)];
            pCurr = [cx(zi),   cy(zi)];
            pNext = [cx(zi+1), cy(zi+1)];

            vPrevCurr = pCurr - pPrev;
            vCurrNext = pNext - pCurr;
            vPrevNext = pNext - pPrev;

            dPrevCurr = norm(vPrevCurr);
            dCurrNext = norm(vCurrNext);
            dPrevNext = norm(vPrevNext);

            isJumpReturn = ...
                dPrevCurr >= jumpThresh && ...
                dCurrNext >= jumpThresh && ...
                dPrevNext <= returnThresh;

            if ~isJumpReturn
                continue
            end

            absDx1 = abs(vPrevCurr(1));
            absDy1 = abs(vPrevCurr(2));
            absDx2 = abs(vCurrNext(1));
            absDy2 = abs(vCurrNext(2));

            isXOnly = ...
                absDx1 >= jumpThresh && ...
                absDx2 >= jumpThresh && ...
                absDy1 <= axisRatioThresh * absDx1 && ...
                absDy2 <= axisRatioThresh * absDx2;

            isYOnly = ...
                absDy1 >= jumpThresh && ...
                absDy2 >= jumpThresh && ...
                absDx1 <= axisRatioThresh * absDy1 && ...
                absDx2 <= axisRatioThresh * absDy2;

            if isXOnly || isYOnly
                isOscObj = true;
                break
            end

        end

        if ~isOscObj
            continue
        end

        targetX = median(cx);
        targetY = median(cy);

        for zi = 1:numel(zList)

            zCurr = zList(zi);

            thisDx = targetX - cx(zi);
            thisDy = targetY - cy(zi);

            stitchObjID{zCurr}(end+1,1) = obj;
            stitchDx{zCurr}(end+1,1) = thisDx;
            stitchDy{zCurr}(end+1,1) = thisDy;

            oscComponentID(end+1,1) = obj;
            oscZ(end+1,1) = zCurr;
            oscDx(end+1,1) = thisDx;
            oscDy(end+1,1) = thisDy;

        end

    end

    oscTable = table(oscComponentID, oscZ, oscDx, oscDy, ...
        'VariableNames', {'ComponentID','Z','StitchCorrectionX','StitchCorrectionY'});

end

function [zDx, zDy, zStepDx, zStepDy, nBeads, beadCount] = identifyZDrift(CC3D, L3D, imSize, nZ, maxArea, minBeadSlices, minMedianCircularity)
% Selects bead candidates from the 3D objects: small (max single-slice
% area below maxArea), present across at least minBeadSlices, and highly
% circular in cross-section on the median slice (bounding-box
% short-side/long-side ratio above minMedianCircularity) -- real cells
% are rarely this compact and consistent, so this reliably isolates the
% fiducial beads. Then, slice by slice, matches each bead's label between
% consecutive Z-slices and takes the MEDIAN centroid displacement across
% all matched beads as that slice-step's drift estimate (median is
% robust to any one mismatched/lost bead); cumulatively summing those
% per-step estimates over Z gives a smooth absolute XY-drift-vs-Z curve.

Y = imSize(1);
X = imSize(2);
Z = imSize(3);
beadID = [];
beadCount = 0;

    for obj = 1:CC3D.NumObjects
        pix = CC3D.PixelIdxList{obj};
        [yy, xx, zz] = ind2sub([Y X Z], pix);
        zList = unique(zz);
        if numel(zList) < minBeadSlices
            continue
        end

        areaPerSlice = zeros(numel(zList),1);
        circPerSlice = nan(numel(zList),1);
        for zi = 1:numel(zList)
            zCurr = zList(zi);
            idxZ = zz == zCurr;
              areaPerSlice(zi) = nnz(idxZ);

                % bounding box around ALL pieces of this object in this slice
                xMin = min(xx(idxZ));
                xMax = max(xx(idxZ));
                yMin = min(yy(idxZ));
                yMax = max(yy(idxZ));

                boxWidth  = xMax - xMin + 1;
                boxHeight = yMax - yMin + 1;

                shortSide = min(boxWidth, boxHeight);
                longSide  = max(boxWidth, boxHeight);

                circPerSlice(zi) = shortSide / max(longSide, eps);
        end

        maxObjArea = max(areaPerSlice);

        % max 2D slice area must be less than designated bead max area
        if maxObjArea >= maxArea
            continue
        end

        medianCircularity = median(circPerSlice, 'omitnan');

        % median circularity across slices must be greater than 0.8
        if isnan(medianCircularity) || medianCircularity <= minMedianCircularity
            continue
        end

        beadCount = beadCount + 1;
        beadID(beadCount,1) = obj;

    end

    zStepDx = zeros(nZ,1);
    zStepDy = zeros(nZ,1);
    nBeads = zeros(nZ,1);
    for z = 2:nZ
        Lprev = L3D(:,:,z-1);
        Lcurr = L3D(:,:,z);
        dxList = [];
        dyList = [];

        for b = 1:beadCount
            id = beadID(b);
            maskPrev = Lprev == id;
            maskCurr = Lcurr == id;
            if ~any(maskPrev(:)) || ~any(maskCurr(:))
                continue
            end

            propsPrev = regionprops(maskPrev, 'Centroid');
            propsCurr = regionprops(maskCurr, 'Centroid');
            cPrev = propsPrev.Centroid;
            cCurr = propsCurr.Centroid;
            dxList(end+1,1) = cCurr(1) - cPrev(1);
            dyList(end+1,1) = cCurr(2) - cPrev(2);
        end

        nBeads(z) = numel(dxList);
        if nBeads(z) > 0
            zStepDx(z) = median(dxList);
            zStepDy(z) = median(dyList);
        end
    end
    zDx = cumsum(zStepDx);
    zDy = cumsum(zStepDy);

end

function [binTileCorr, rawTileCorrByChannel] = applyObjectCorrections( ...
    binTile, L3D, rawTileByChannel, rawFields, ...
    zDx, zDy, stitchObjID, stitchDx, stitchDy)
% Moves every 2D object slice-by-slice by the combination of the bead-
% based Z-drift shift (baseXShift/baseYShift = -zDx/-zDy, negated because
% the drift curve describes how much the beads moved, so objects need to
% move the opposite way to compensate) and, if that object was flagged as
% stitching-oscillating at this Z, its own extra stitching correction
% looked up by matching its most common label within the slice. Before
% shifting, small fully-enclosed holes in the object's local footprint
% are filled by hand (finding background-connected-components under
% maxHoleArea that don't touch the local crop border -- equivalent to
% imfill but limited to a small local crop for speed) so a partially
% thresholded object doesn't end up with holes punched through it by the
% shift. Every raw channel's output stack is pre-filled with the
% per-slice local background median BEFORE any object is placed, so that
% once objects are moved away from their original positions the vacated
% area reads as background rather than zero/black.

    [Y, X, nZ] = size(binTile);

    maxHoleArea = 10;
    boxPad = 1;

    binTileCorr = false(size(binTile));
    rawTileCorrByChannel = struct();

    for rf = 1:numel(rawFields)

        fieldName = rawFields{rf};
        rawStack = rawTileByChannel.(fieldName);

        rawTileCorrByChannel.(fieldName) = zeros(size(rawStack), ...
            'like', rawStack);

        for z = 1:nZ

            slice = rawStack(:,:,z);
            binSlice = binTile(:,:,z);

            bgPixels = slice(~binSlice);

            if isempty(bgPixels)
                bgVal = median(slice(:));
            else
                bgVal = median(bgPixels);
            end

            rawTileCorrByChannel.(fieldName)(:,:,z) = bgVal;

        end

    end

    for z = 1:nZ

        CC2D = bwconncomp(binTile(:,:,z));
        Lslice = L3D(:,:,z);

        baseXShift = round(-zDx(z));
        baseYShift = round(-zDy(z));

        for obj = 1:CC2D.NumObjects

            pix = CC2D.PixelIdxList{obj};
            [yy, xx] = ind2sub([Y X], pix);

            idVals = Lslice(pix);
            idVals(idVals == 0) = [];

            thisStitchX = 0;
            thisStitchY = 0;

            if ~isempty(idVals)

                compID = mode(double(idVals));
                stitchIdx = find(stitchObjID{z} == compID, 1);

                if ~isempty(stitchIdx)
                    thisStitchX = stitchDx{z}(stitchIdx);
                    thisStitchY = stitchDy{z}(stitchIdx);
                end

            end

            xShift = baseXShift + round(thisStitchX);
            yShift = baseYShift + round(thisStitchY);

            % local cropped object mask
            r1 = max(1, min(yy) - boxPad);
            r2 = min(Y, max(yy) + boxPad);
            c1 = max(1, min(xx) - boxPad);
            c2 = min(X, max(xx) + boxPad);

            localMask = false(r2-r1+1, c2-c1+1);

            yyLocal = yy - r1 + 1;
            xxLocal = xx - c1 + 1;
            localPix = sub2ind(size(localMask), yyLocal, xxLocal);
            localMask(localPix) = true;

            % no-imfill hole detection:
            % background components fully enclosed inside localMask are holes
            bgMask = ~localMask;
            CCbg = bwconncomp(bgMask, 4);

            fillLocal = false(size(localMask));
            [localH, localW] = size(localMask);

            for hObj = 1:CCbg.NumObjects

                bgPix = CCbg.PixelIdxList{hObj};

                if numel(bgPix) > maxHoleArea
                    continue
                end

                [bgY, bgX] = ind2sub([localH localW], bgPix);

                touchesBorder = any(bgY == 1 | bgY == localH | ...
                                    bgX == 1 | bgX == localW);

                if touchesBorder
                    continue
                end

                fillLocal(bgPix) = true;

            end

            localMask = localMask | fillLocal;

            [yyFilledLocal, xxFilledLocal] = find(localMask);

            yyFilled = yyFilledLocal + r1 - 1;
            xxFilled = xxFilledLocal + c1 - 1;

            pixFilled = sub2ind([Y X], yyFilled, xxFilled);

            yy2 = yyFilled + yShift;
            xx2 = xxFilled + xShift;

            keep = yy2 >= 1 & yy2 <= Y & xx2 >= 1 & xx2 <= X;

            if any(keep)

                pix2 = sub2ind([Y X], yy2(keep), xx2(keep));

                temp = binTileCorr(:,:,z);
                temp(pix2) = true;
                binTileCorr(:,:,z) = temp;

                for rf = 1:numel(rawFields)

                    fieldName = rawFields{rf};

                    rawSlice = rawTileByChannel.(fieldName)(:,:,z);
                    rawVals = rawSlice(pixFilled);

                    tempRaw = rawTileCorrByChannel.(fieldName)(:,:,z);
                    tempRaw(pix2) = rawVals(keep);
                    rawTileCorrByChannel.(fieldName)(:,:,z) = tempRaw;

                end

            end

        end

    end

end

function [binTileAligned, rawTileAlignedByChannel] = alignObjectsToTopZ(binTileCorr, rawTileCorrByChannel, rawFields)
% Shifts every final 3D object rigidly along Z so its own topmost slice
% lands on the top of the stack (Z = Z_max) -- removing timepoint-to-
% timepoint variation in where an object sits in Z (e.g. from refocusing)
% by anchoring every object to the same reference plane. As in
% applyObjectCorrections, each output raw channel is pre-filled with the
% per-slice background median before any object is placed, so the area
% an object moved away from reads as background.

    [Y, X, Z] = size(binTileCorr);

    binTileAligned = false(size(binTileCorr));

    rawTileAlignedByChannel = struct();

    for rf = 1:numel(rawFields)

        fieldName = rawFields{rf};
        rawStack = rawTileCorrByChannel.(fieldName);

        rawTileAlignedByChannel.(fieldName) = zeros(size(rawStack), 'like', rawStack);

        for z = 1:Z

            binSlice = binTileCorr(:,:,z);
            slice = rawStack(:,:,z);

            bgPixels = slice(~binSlice);

            if isempty(bgPixels)
                bgVal = median(slice(:));
            else
                bgVal = median(bgPixels);
            end

            rawTileAlignedByChannel.(fieldName)(:,:,z) = bgVal;

        end

    end

    CC3Dcorr = bwconncomp(binTileCorr, 26);

    for obj = 1:CC3Dcorr.NumObjects

        pix = CC3Dcorr.PixelIdxList{obj};
        [yy, xx, zz] = ind2sub([Y X Z], pix);

        topZ = max(zz);
        zShift = Z - topZ;

        zz2 = zz + zShift;

        pix2 = sub2ind([Y X Z], yy, xx, zz2);

        binTileAligned(pix2) = true;

        for rf = 1:numel(rawFields)

            fieldName = rawFields{rf};

            rawVals = rawTileCorrByChannel.(fieldName)(pix);
            rawTileAlignedByChannel.(fieldName)(pix2) = rawVals;

        end

    end

end

function binTileOut = removeShakyObjects(binTileIn, maxProjRatioThresh, shakyThresh, smallArea, areaThreshold)
% Final QC pass on objects that survived Z-drift/stitching correction.
% Trivially small objects (max-projection area below smallArea) are
% dropped outright without being judged for "shakiness" -- too little
% signal to trust either metric. Very large objects (above areaThreshold)
% are exempted/protected from removal, since both metrics below get
% noisier and the cost of wrongly discarding a big real structure is
% higher. For everything in between, an object is considered shaky (and
% removed) if EITHER:
%   - its max-Z-projection footprint is much bigger than its largest
%     single-slice area (maxProjRatioThresh) -- i.e. it wanders around in
%     XY across Z rather than sitting still, or
%   - its per-slice centroid strays more than shakyThresh from its own
%     median centroid position at some point in Z.

    binTileOut = binTileIn;
    [Y, X, Z] = size(binTileIn);
    CC3D = bwconncomp(binTileIn, 26);

    nRemovedSmall = 0;
    nRemovedRatio = 0;
    nRemovedCentroid = 0;
    nRemovedBoth = 0;

    for obj = 1:CC3D.NumObjects

        pix = CC3D.PixelIdxList{obj};
        [yy, xx, zz] = ind2sub([Y X Z], pix);

        % area of object in each Z slice
        sliceAreas = accumarray(zz, 1, [Z 1]);
        maxSliceArea = max(sliceAreas);

        % area of max projection of this object
        xyPix = sub2ind([Y X], yy, xx);
        maxProjectionArea = numel(unique(xyPix));

        % first remove tiny projected objects
        % do not count these as shaky/moving
        if maxProjectionArea < smallArea
            binTileOut(pix) = false;
            nRemovedSmall = nRemovedSmall + 1;
            continue
        end

        % protect large objects
        if maxProjectionArea > areaThreshold
            continue
        end

        % max-projection spread metric
        maxProjRatio = maxProjectionArea / maxSliceArea;

        % centroid movement metric
        zList = unique(zz);
        cx = zeros(numel(zList), 1);
        cy = zeros(numel(zList), 1);

        for zi = 1:numel(zList)

            ztmp = zList(zi);
            idxZ = zz == ztmp;

            cx(zi) = mean(xx(idxZ));
            cy(zi) = mean(yy(idxZ));

        end

        medX = median(cx);
        medY = median(cy);

        distFromMedian = hypot(cx - medX, cy - medY);
        maxCentroidDist = max(distFromMedian);

        isShakyByRatio = maxProjRatio > maxProjRatioThresh;
        isShakyByCentroid = maxCentroidDist > shakyThresh;

        if isShakyByRatio || isShakyByCentroid

            binTileOut(pix) = false;

            if isShakyByRatio && isShakyByCentroid
                nRemovedBoth = nRemovedBoth + 1;
            elseif isShakyByRatio
                nRemovedRatio = nRemovedRatio + 1;
            elseif isShakyByCentroid
                nRemovedCentroid = nRemovedCentroid + 1;
            end

        end

    end

    fprintf(['Removed %d small objects, %d by max-projection ratio, ', ...
        '%d by centroid movement, %d by both criteria\n'], ...
        nRemovedSmall, nRemovedRatio, nRemovedCentroid, nRemovedBoth);

end

function writeTiffStack(stack, outName)
% Writes a 3D array as a multi-page TIFF, one slice per page. Deletes any
% existing file of the same name first since imwrite would otherwise
% append onto stale pages left over from a previous run. Logical (binary
% mask) input is converted to 0/255 uint8 for viewing; anything else is
% written as-is.

    if isfile(outName)
        delete(outName);
    end

    for z = 1:size(stack,3)

        if islogical(stack)
            slice = uint8(stack(:,:,z)) * 255;
        else
            slice = stack(:,:,z);
        end

        if z == 1
            imwrite(slice, outName, 'tif', 'Compression', 'none');
        else
            imwrite(slice, outName, 'tif', ...
                'WriteMode', 'append', 'Compression', 'none');
        end

    end

end

function addImageJTimeMetadata(fname, nT)
% Stamps the given TIFF with an ImageJ hyperstack ImageDescription tag
% declaring it as a single-channel, single-slice, nT-frame time series.
% Needed because plain sequential imwrite appends don't carry any
% dimension metadata, so without this ImageJ would open the growing
% max-projection file as an unlabeled stack instead of a proper T-series.
% Must be re-written every time a frame is appended, since nT changes.
t = Tiff(fname, 'r+');
desc = sprintf(['ImageJ=1.53\n' ...
                'images=%d\n' ...
                'channels=1\n' ...
                'slices=1\n' ...
                'frames=%d\n' ...
                'hyperstack=true\n' ...
                'mode=grayscale\n' ...
                'loop=false\n'], nT, nT);
setTag(t, 'ImageDescription', desc);
close(t);
end

function split = makeSplit(w,h,q)
% Divides a w-by-h image into a q-tile grid (q must be a perfect square;
% n = sqrt(q) tiles per axis) so each tile's Z-stack can be loaded and
% processed independently -- keeps memory use bounded regardless of how
% large the full stitched field of view is. Returns one {rows, cols, all
% Z} pixel-region spec per tile, in the format tiffreadVolume's
% 'PixelRegion' option expects.

n = sqrt(q);              % number of tiles per dimension
x = round(linspace(1,w+1,n+1));
y = round(linspace(1,h+1,n+1));

split = cell(q,3);
k = 1;

for r = 1:n
    for c = 1:n
        split{k,1} = [y(r) 1 y(r+1)-1];   % rows (Y)
        split{k,2} = [x(c) 1 x(c+1)-1];   % cols (X)
        split{k,3} = [1 Inf];             % all slices
        k = k + 1;
    end
end
end

function tStart = startTimer(~)
% Starts a stopwatch for a named processing step; paired with stopTimer.
    tStart = tic;

end

function stopTimer(stepName, tStart)
% Prints elapsed time for a step started with startTimer, timestamped so
% the console log doubles as a timeline of a long batch run.

    elapsedSec = toc(tStart);

    fprintf('[%s] DONE:  %s | %.2f sec | %.2f min\n', ...
        char(datetime("now","Format","yyyy-MM-dd HH:mm:ss")), ...
        stepName, elapsedSec, elapsedSec/60);

end
