function SegmentSpikeInAttachment(path)
% Author: Hasreet Gill, Bassler Lab, HHMI/Princeton University
% Revised with Claude (Sonnet 5, Anthropic) 

% Cluster batch driver (SLURM, "della") for a spike-in attachment assay:
% for every (Date, XY-position) combination found in PATH, segments each
% timepoint's 561-channel colony image into Voronoi-separated colony
% "zones", cross-references externally-supplied TrackMate track/spot
% CSVs (tracks presumed already run on the corresponding images -- see
% the loading note below) to find short, low-motion "attached" particle
% tracks, assigns each attached track to whichever colony zone it falls
% inside, and writes per-position and dataset-wide summary spreadsheets
% plus a QC overlay PNG per timepoint.
%
% Inputs expected directly under PATH: raw "*_T<n>_XY<n>_..._20X...561.tif"
% Z-stacks (NOT the "Binary_"-prefixed ones -- unlike the culture
% pipeline, this script re-thresholds the raw volume itself rather than
% consuming any pre-made mask), matching "*tracks.csv"/"*spots.csv"
% TrackMate exports, and metadata.xlsx / genotypes.xlsx lookup tables
% (metadata.xlsx here has one row per date, a "Time step"/"Pixel"
% column, and one column per XY position whose value is an
% underscore-joined "StrainCode_Hours_Replicate_Arabinose[_Nspd]" string
% -- see getMetaInfo below).
%
% NOTE: ims/trks/spts (the image, tracks, and spots file lists) are each
% independently sorted by filename and then assumed to line up 1:1 in
% that order -- there's no join by matching filename stem. If any of the
% three globs picks up a different count or ordering than the others,
% an image would silently get paired with the wrong tracks/spots file.

path = char(path);
nCores = str2double(getenv('SLURM_CPUS_PER_TASK'));
if isnan(nCores) || nCores < 1
    nCores = 1;   % fall back to a single worker if the SLURM env var is missing/invalid
end

% Same size-to-allocation local parpool pattern as the culture pipeline's
% cluster driver: reuse a correctly-sized pool, or tear down and rebuild
% one that isn't.
pobj = gcp('nocreate');
if isempty(pobj)
    pc = parcluster('local');
    pc.JobStorageLocation = '/tmp/';
    parpool(pc, nCores);
elseif pobj.NumWorkers ~= nCores
    delete(pobj);
    pc = parcluster('local');
    pc.JobStorageLocation = '/tmp/';
    parpool(pc, nCores);
end

areaThresh = 150;  % spot-area cutoff (pixels) used to split attached particles into small/large counts

% load metadata
meta = readtable(fullfile(path, 'metadata.xlsx'), 'VariableNamingRule', 'preserve');
meta.Date = string(meta.Date);

geno = readtable(fullfile(path, 'genotypes.xlsx'), 'VariableNamingRule', 'preserve');
strainMap = containers.Map('KeyType','char','ValueType','char');
for i = 1:height(geno)
    code = strtrim(string(geno{i,1}));
    name = strtrim(string(geno{i,2}));
    strainMap(char(code)) = char(name);
end

% load files
ims = dir(fullfile(path, '*561.tif'));
ims(contains(string({ims.name}), "Binary")) = [];
trks = dir(fullfile(path, '*tracks.csv'));
spts = dir(fullfile(path, '*spots.csv'));
[~, ord] = sort(string({ims.name}));
ims = ims(ord);
[~, ord] = sort(string({trks.name}));
trks = trks(ord);
[~, ord] = sort(string({spts.name}));
spts = spts(ord);

nIm = numel(ims);
imDate = strings(nIm,1);
imT    = nan(nIm,1);
imXY   = strings(nIm,1);
for i = 1:nIm
    % Filenames are parsed as "<Date>_..._T<n>_XY<n>_...": Date is
    % everything before the first underscore, T and XY are pulled out by
    % their surrounding literal tags.
    nm = string(ims(i).name);
    imDate(i) = extractBefore(nm, "_");
    tStr = extractBetween(nm, "_T", "_XY");
    imT(i) = str2double(tStr);
    xyStr = extractBetween(nm, "_XY", "_");
    imXY(i) = "XY" + xyStr;
end

% groups for workers
% Each (Date, XY) pair is one continuous timelapse position, processed as
% one independent unit -- these are what get distributed across workers.
groupTab = unique(table(imDate, imXY, 'VariableNames', {'Date','XY'}), 'rows');

%% main loop
SummaryCell = cell(height(groupTab),1);
parfor g = 1:height(groupTab)

    dateStr = groupTab.Date(g);
    xyStr   = groupTab.XY(g);
    fprintf('Processing %s %s\n', dateStr, xyStr);

    metaInfo = getMetaInfo(meta, strainMap, dateStr, xyStr);
    idxGroup = find(imDate == dateStr & imXY == xyStr);

    SpotsXY = table();
    ObjectsXY = table();
    Summary_g = table();

    for k = 1:length(idxGroup)

        i = idxGroup(k);
        Tval = imT(i);
        fprintf('  T = %d\n', Tval);

        imFile  = fullfile(path, ims(i).name);
        trkFile = fullfile(path, trks(i).name);
        sptFile = fullfile(path, spts(i).name);

        result = analyzeOneTimepoint( ...
            imFile, trkFile, sptFile, ...
            dateStr, Tval, xyStr, ...
            metaInfo, areaThresh, path);

        if isempty(SpotsXY)
            SpotsXY = result.SpotsSummary;
        else
            SpotsXY = [SpotsXY; result.SpotsSummary];
        end

        if isempty(ObjectsXY)
            ObjectsXY = result.ObjectSummary;
        else
            % ObjectSummary carries a variable NUMBER of "TrackJ_..."
            % columns per timepoint (as many as that timepoint's busiest
            % object has tracks -- see analyzeOneTimepoint), so table
            % vertcat can't be used directly. This pads whichever side is
            % missing a column with NaN (numeric) or blank strings
            % (text) before aligning both tables to the same column set
            % and column order, then concatenates.
            vars1 = string(ObjectsXY.Properties.VariableNames);
            vars2 = string(result.ObjectSummary.Properties.VariableNames);
            allVars = unique([vars1 vars2], 'stable');

            for v = allVars
                if ~ismember(v, vars1)
                    if isstring(result.ObjectSummary.(v)) || iscell(result.ObjectSummary.(v))
                        ObjectsXY.(v) = strings(height(ObjectsXY),1);
                    else
                        ObjectsXY.(v) = nan(height(ObjectsXY),1);
                    end
                end

                if ~ismember(v, vars2)
                    if isstring(ObjectsXY.(v)) || iscell(ObjectsXY.(v))
                        result.ObjectSummary.(v) = strings(height(result.ObjectSummary),1);
                    else
                        result.ObjectSummary.(v) = nan(height(result.ObjectSummary),1);
                    end
                end
            end

            ObjectsXY = ObjectsXY(:, cellstr(allVars));
            result.ObjectSummary = result.ObjectSummary(:, cellstr(allVars));
            ObjectsXY = [ObjectsXY; result.ObjectSummary];
        end

        Summary_g = [Summary_g; result.ImageSummaryRow];
    end

    SpotsXY   = sortrows(SpotsXY,   'T');
    ObjectsXY = sortrows(ObjectsXY, 'T');

    writetable(SpotsXY, fullfile(path, "SpotsSummary_" + dateStr + "_" + xyStr + ".xlsx"));
    writetable(ObjectsXY, fullfile(path, "ObjectSummary_" + dateStr + "_" + xyStr + ".xlsx"));

    SummaryCell{g} = Summary_g;
end

SummaryAll = vertcat(SummaryCell{:});
SummaryAll = sortrows(SummaryAll, {'Date','XY','T'});
writetable(SummaryAll, fullfile(path, 'ImageSummary_All.xlsx'));

disp('Done.');

end

function result = analyzeOneTimepoint(imFile, trkFile, sptFile, dateStr, Tval, xyStr, metaInfo, areaThresh, path)
% Segments one timepoint's colonies, filters and spatially assigns that
% timepoint's TrackMate tracks/spots to those colonies, and returns three
% tables: one row per attached spot (SpotsSummary), one row per colony
% zone (ObjectSummary), and one summary row for the whole image
% (ImageSummaryRow).

%% ----- image -----
im = tiffreadVolume(imFile);
maxim = max(im, [], 3);

% Adaptive (image-specific) intensity threshold, same percentile
% convention used elsewhere in this pipeline family.
thr = prctile(im(:), 99);
maxbin = maxim > thr;
maxbin = bwareafilt(maxbin, [50 Inf]);      % general noise floor
maxcol = bwareafilt(maxbin, [500 Inf]);     % stricter floor: only genuine colony-scale blobs

% Generous halo around each colony -- the "zone of interest" that a
% nearby attached particle is considered to belong to, even if it sits
% well outside the colony's own footprint.
dil = imdilate(maxcol, strel("disk",100));

mcols = regionprops(maxcol, "Centroid");
if isempty(mcols)
    ctrs = zeros(0,2);
else
    ctrs = cat(1, mcols.Centroid);
end

if size(ctrs,1) > 2
    % Voronoi-partition the colony centroids and carve those boundary
    % lines out of the dilated halo mask, so two nearby colonies' 100px
    % halos don't merge into one blob (and a spot near the midpoint
    % between them isn't ambiguously shared) -- each colony's halo is
    % capped at its own side of the Voronoi cell.
    [V, C] = voronoin(ctrs);
    vorBW = false(size(dil));
    H = size(vorBW,1);
    W = size(vorBW,2);

    for ci = 1:numel(C)
        vidx = C{ci};
        if any(vidx == 1), continue; end   % vertex 1 = point at infinity -> unbounded cell, skip
        poly = V(vidx, :);

        for kk = 1:size(poly,1)
            p1 = poly(kk,:);
            p2 = poly(mod(kk, size(poly,1)) + 1, :);

            n = max(abs(round(p2 - p1))) + 1;
            xs = round(linspace(p1(1), p2(1), n));
            ys = round(linspace(p1(2), p2(2), n));

            good = xs>=1 & xs<=W & ys>=1 & ys<=H;
            xs = xs(good);
            ys = ys(good);

            vorBW(sub2ind([H W], ys, xs)) = true;
        end
    end

    dil = dil & ~imdilate(vorBW, strel("disk",2));
end

% Final Voronoi-separated colony zones -- these labels are what spots
% get assigned to below, and boundaries drawn on the QC overlay.
L = bwlabel(dil);
b = bwboundaries(dil);

%% ----- tracking -----
trkdata = readtable(trkFile);
sptdata = readtable(sptFile);

% Keep only tracks that persisted >= 2 frames but barely moved (small
% end-to-end displacement AND small total path length) -- the signature
% of a particle that ATTACHED and stayed put, as opposed to one still
% freely diffusing (large displacement/distance) or a single-frame
% detection blip (< 2 frames).
nFrames = trkdata.TRACK_STOP - trkdata.TRACK_START + 1;
keepTrk = (nFrames >= 2) & (trkdata.TRACK_DISPLACEMENT <= 15) & (trkdata.TOTAL_DISTANCE_TRAVELED <= 20);
trk_filt = trkdata(keepTrk,:);

keepSpt = ismember(sptdata.TRACK_ID, trk_filt.TRACK_ID);
spt_filt = sptdata(keepSpt,:);

[~, baseName] = fileparts(imFile);
imageName = string(extractBefore(baseName, "_20X"));

%% ----- preliminary spots summary -----
% One row per surviving track, using that track's EARLIEST frame as its
% representative position/area (its position at first attachment).
spotTrackIDs = unique(spt_filt.TRACK_ID);
nSpotTracks = numel(spotTrackIDs);

ImageName = repmat(imageName, nSpotTracks, 1);
DateCol = repmat(string(dateStr), nSpotTracks, 1);
TCol = repmat(Tval, nSpotTracks, 1);
XYCol = repmat(string(xyStr), nSpotTracks, 1);
StrainCodeCol = repmat(metaInfo.StrainCode, nSpotTracks, 1);
GenotypeCol = repmat(metaInfo.Genotype, nSpotTracks, 1);
HoursCol = repmat(metaInfo.Hours, nSpotTracks, 1);
ReplicateCol = repmat(metaInfo.Replicate, nSpotTracks, 1);
ArabinoseCol = repmat(metaInfo.Arabinose, nSpotTracks, 1);
NspdCol = repmat(metaInfo.Nspd, nSpotTracks, 1);   % NEW
TimeStepCol = repmat(metaInfo.TimeStep_m, nSpotTracks, 1);
PixelCol = repmat(metaInfo.Pixel_um, nSpotTracks, 1);

TrackID = nan(nSpotTracks,1);
SpotArea = nan(nSpotTracks,1);
X = nan(nSpotTracks,1);
Y = nan(nSpotTracks,1);
DistToNearestDilCentroid = nan(nSpotTracks,1);
TrackDisplacement = nan(nSpotTracks,1);
TrackDuration = nan(nSpotTracks,1);

for ss = 1:nSpotTracks
    tid = spotTrackIDs(ss);
    idx = (spt_filt.TRACK_ID == tid);

    frameVals = spt_filt.FRAME(idx);
    [~, j] = min(frameVals);

    xtmp = spt_filt.POSITION_X(idx);
    ytmp = spt_filt.POSITION_Y(idx);

    x0 = xtmp(j);
    y0 = ytmp(j);

    TrackID(ss) = tid;
    X(ss) = x0;
    Y(ss) = y0;

    if ismember("AREA", string(spt_filt.Properties.VariableNames))
        atmp = spt_filt.AREA(idx);
        SpotArea(ss) = atmp(j);
    end

    if ~isempty(ctrs)
        d = sqrt((ctrs(:,1) - x0).^2 + (ctrs(:,2) - y0).^2);
        DistToNearestDilCentroid(ss) = min(d);
    end

    trkRow = find(trk_filt.TRACK_ID == tid, 1);
    if ~isempty(trkRow)
        TrackDisplacement(ss) = trk_filt.TRACK_DISPLACEMENT(trkRow);
        TrackDuration(ss) = trk_filt.TRACK_STOP(trkRow) - trk_filt.TRACK_START(trkRow) + 1;
    end
end

SpotsSummary = table(DateCol, TCol, XYCol, ImageName, ...
    StrainCodeCol, GenotypeCol, HoursCol, ReplicateCol, ArabinoseCol, NspdCol, TimeStepCol, PixelCol, ...
    TrackID, SpotArea, X, Y, DistToNearestDilCentroid, TrackDisplacement, TrackDuration, ...
    'VariableNames', {'Date','T','XY','ImageName', ...
    'StrainCode','Genotype','Hours','Replicate','Arabinose', 'Nspd', 'TimeStep_m','Pixel_um', ...
    'TrackID','SpotArea','X','Y','DistToNearestDilCentroid', ...
    'TrackDisplacement','TrackDuration'});

%% ----- image-summary spot columns -----
totalAttachedSpots_All = height(SpotsSummary);

goodArea = ~isnan(SpotsSummary.SpotArea);
numAttachedSpotsLT150 = sum(SpotsSummary.SpotArea(goodArea) < areaThresh);
numAttachedSpotsGE150 = sum(SpotsSummary.SpotArea(goodArea) >= areaThresh);

goodDist = ~isnan(SpotsSummary.DistToNearestDilCentroid);
if any(goodDist)
    meanDistToNearestDilCentroid_All = mean(SpotsSummary.DistToNearestDilCentroid(goodDist), 'omitnan');
else
    meanDistToNearestDilCentroid_All = NaN;
end

goodLT = ~isnan(SpotsSummary.DistToNearestDilCentroid) & ...
         ~isnan(SpotsSummary.SpotArea) & ...
         (SpotsSummary.SpotArea < areaThresh);
if any(goodLT)
    meanDistToNearestDilCentroid_LT150 = mean(SpotsSummary.DistToNearestDilCentroid(goodLT), 'omitnan');
else
    meanDistToNearestDilCentroid_LT150 = NaN;
end

%% ----- assign tracks to objects -----
% Which colony zone (label in L) does each surviving spot physically sit
% in? Round/clamp to valid pixel coordinates first in case a spot sits
% exactly on/past the image edge.
x = round(spt_filt.POSITION_X);
y = round(spt_filt.POSITION_Y);
x = max(1, min(size(L,2), x));
y = max(1, min(size(L,1), y));

if isempty(x)
    spt_filt.OBJECT_ID = zeros(0,1);
else
    spt_filt.OBJECT_ID = L(sub2ind(size(L), y, x));
end

uTracks = unique(spt_filt.TRACK_ID);
trackObject = zeros(size(uTracks));

for t = 1:numel(uTracks)
    idx = (spt_filt.TRACK_ID == uTracks(t));
    objs = spt_filt.OBJECT_ID(idx);
    objs = objs(objs > 0);

    if isempty(objs)
        trackObject(t) = 0;
    else
        trackObject(t) = mode(objs);   % majority zone across this track's few frames
    end
end

% Tracks that never land inside any colony zone (trackObject == 0) are
% dropped -- a spike-in particle not near any colony isn't "attached".
keepTrack = trackObject > 0;
keepTrackIDs2 = uTracks(keepTrack);
keepTrk2 = ismember(trk_filt.TRACK_ID, keepTrackIDs2);

trk_filt = trk_filt(keepTrk2,:);
trk_filt.OBJECT_ID = zeros(height(trk_filt),1);

for r = 1:height(trk_filt)
    tid = trk_filt.TRACK_ID(r);
    trk_filt.OBJECT_ID(r) = trackObject(uTracks == tid);
end

% SpotsSummary is re-filtered to match (a second, narrower filter beyond
% the displacement/duration one above) and stamped with the same
% per-track OBJECT_ID for use below.
SpotsSummary = SpotsSummary(ismember(SpotsSummary.TrackID, trk_filt.TRACK_ID), :);
SpotsSummary.OBJECT_ID = zeros(height(SpotsSummary),1);

for r = 1:height(SpotsSummary)
    tid = SpotsSummary.TrackID(r);
    trkRow = find(trk_filt.TRACK_ID == tid, 1);
    if ~isempty(trkRow)
        SpotsSummary.OBJECT_ID(r) = trk_filt.OBJECT_ID(trkRow);
    end
end

%% ----- overlay image for cluster (no figures) -----
% Hand-built RGB QC image (cluster nodes are headless, so this avoids
% MATLAB's normal figure/plot functions entirely): grayscale background
% from the general mask, true colonies in blue, colony-zone boundaries in
% red, and every retained attached spot as a filled yellow disc.
overlay = zeros(size(maxbin,1), size(maxbin,2), 3, 'uint8');
overlay(:,:,1) = uint8(maxbin) * 255;
overlay(:,:,2) = uint8(maxbin) * 255;
overlay(:,:,3) = uint8(maxbin) * 255;

R = overlay(:,:,1);
G = overlay(:,:,2);
B = overlay(:,:,3);

R(maxcol) = 0;
G(maxcol) = 0;
B(maxcol) = 255;

overlay(:,:,1) = R;
overlay(:,:,2) = G;
overlay(:,:,3) = B;

boundaryMask = false(size(maxbin));
for kk = 1:length(b)
    boundary = b{kk};
    rr = boundary(:,1);
    cc = boundary(:,2);

    good = rr >= 1 & rr <= size(maxbin,1) & ...
           cc >= 1 & cc <= size(maxbin,2);
    rr = rr(good);
    cc = cc(good);

    boundaryMask(sub2ind(size(maxbin), rr, cc)) = true;
end

lineWidth = 3;
boundaryMask = imdilate(boundaryMask, strel("disk", lineWidth));

R = overlay(:,:,1);
G = overlay(:,:,2);
B = overlay(:,:,3);

R(boundaryMask) = 255;
G(boundaryMask) = 0;
B(boundaryMask) = 0;

overlay(:,:,1) = R;
overlay(:,:,2) = G;
overlay(:,:,3) = B;

for tt = 1:height(SpotsSummary)
    x0 = round(SpotsSummary.X(tt));
    y0 = round(SpotsSummary.Y(tt));

    rad = 16;
    r1 = max(1, y0-rad):min(size(maxbin,1), y0+rad);
    c1 = max(1, x0-rad):min(size(maxbin,2), x0+rad);

    for rr = r1
        for cc = c1
            if (rr-y0)^2 + (cc-x0)^2 <= rad^2
                overlay(rr,cc,1) = 255;
                overlay(rr,cc,2) = 255;
                overlay(rr,cc,3) = 0;
            end
        end
    end
end

outName = fullfile(path, ['FiltTracks_' char(imageName) '.png']);
imwrite(overlay, outName);

%% ----- object summary -----
% One row per colony zone: its own true colony area vs. the full
% dilated-zone area (their difference is the surrounding "moat" that
% attached particles are searched within), how many surviving
% tracks/spots landed in it, and what fraction of its moat area is
% covered by attached-particle area.
rp = regionprops(L, "Centroid", "PixelIdxList");
nObj = numel(rp);

objAreaMaxcol = zeros(nObj,1);
dilArea = zeros(nObj,1);

for oo = 1:nObj
    dilArea(oo) = numel(rp(oo).PixelIdxList);
    objAreaMaxcol(oo) = nnz(maxcol(rp(oo).PixelIdxList));
end

if isempty(trk_filt)
    tracksPerObj = zeros(nObj,1);
else
    tracksPerObj = accumarray(trk_filt.OBJECT_ID, 1, [nObj 1], @sum, 0);
end

numObjectsGT1Track = sum(tracksPerObj > 1);
numObjectsGT2Track = sum(tracksPerObj > 2);

NumSpots = zeros(nObj,1);
for oo = 1:nObj
    NumSpots(oo) = sum(SpotsSummary.OBJECT_ID == oo);
end

DateObj = repmat(string(dateStr), nObj, 1);
TObj = repmat(Tval, nObj, 1);
XYObj = repmat(string(xyStr), nObj, 1);
ImageObj = repmat(imageName, nObj, 1);
StrainObj = repmat(metaInfo.StrainCode, nObj, 1);
GenotypeObj = repmat(metaInfo.Genotype, nObj, 1);
HoursObj = repmat(metaInfo.Hours, nObj, 1);
RepObj = repmat(metaInfo.Replicate, nObj, 1);
ArabinoseObj = repmat(metaInfo.Arabinose, nObj, 1);
NspdObj = repmat(metaInfo.Nspd, nObj, 1);   % NEW
TimeStepObj = repmat(metaInfo.TimeStep_m, nObj, 1);
PixelObj = repmat(metaInfo.Pixel_um, nObj, 1);

Centroid_X = nan(nObj,1);
Centroid_Y = nan(nObj,1);
Area = objAreaMaxcol;
DilatedArea = dilArea;
AreaAroundObject = dilArea - objAreaMaxcol;
NumTracks = tracksPerObj;
FractionAreaFilledWithSpots = nan(nObj,1);

for oo = 1:nObj
    Centroid_X(oo) = rp(oo).Centroid(1);
    Centroid_Y(oo) = rp(oo).Centroid(2);

    idxObj = SpotsSummary.OBJECT_ID == oo;
    sumSpotArea = sum(SpotsSummary.SpotArea(idxObj), 'omitnan');

    denom = AreaAroundObject(oo);
    if denom > 0
        FractionAreaFilledWithSpots(oo) = sumSpotArea / denom;
    else
        FractionAreaFilledWithSpots(oo) = NaN;
    end
end

ObjectSummary = table(DateObj, TObj, XYObj, ImageObj, StrainObj, GenotypeObj, ...
    HoursObj, RepObj, ArabinoseObj, NspdObj, TimeStepObj, PixelObj, ...
    Centroid_X, Centroid_Y, Area, DilatedArea, AreaAroundObject, ...
    NumTracks, NumSpots, FractionAreaFilledWithSpots, ...
    'VariableNames', {'Date','T','XY','ImageName','StrainCode','Genotype', ...
    'Hours','Replicate','Arabinose','Nspd','TimeStep_m','Pixel_um', ...
    'Centroid_X','Centroid_Y','Area','DilatedArea','AreaAroundObject', ...
    'NumTracks','NumSpots','FractionAreaFilledWithSpots'});

% Per-object, per-track detail columns (TrackJ_SpotArea/Displacement/
% DurationFrames, J = 1..however many tracks this object has, ordered by
% ascending TRACK_ID). Different objects end up with different numbers
% of these columns -- that's what the caller's vertcat-alignment logic
% has to reconcile across timepoints.
for oo = 1:nObj
    if isempty(trk_filt), continue; end

    rows = find(trk_filt.OBJECT_ID == oo);
    if isempty(rows), continue; end

    [~, ord] = sort(trk_filt.TRACK_ID(rows), "ascend");
    rows = rows(ord);

    for j = 1:numel(rows)
        rtrk = rows(j);
        tid  = trk_filt.TRACK_ID(rtrk);

        idxSp = SpotsSummary.TrackID == tid;
        if any(idxSp)
            spotAreaJ = SpotsSummary.SpotArea(find(idxSp,1,'first'));
        else
            spotAreaJ = NaN;
        end

        dispj = trk_filt.TRACK_DISPLACEMENT(rtrk);
        durj  = trk_filt.TRACK_STOP(rtrk) - trk_filt.TRACK_START(rtrk) + 1;

        ObjectSummary.("Track" + j + "_SpotArea")(oo) = spotAreaJ;
        ObjectSummary.("Track" + j + "_Displacement")(oo) = dispj;
        ObjectSummary.("Track" + j + "_DurationFrames")(oo) = durj;
    end
end

%% ----- image summary row -----
totalNumObjects = nObj;

if isempty(trk_filt)
    numObjectsWithTracks = 0;
else
    numObjectsWithTracks = numel(unique(trk_filt.OBJECT_ID));
end

proportionWithTracks = numObjectsWithTracks / totalNumObjects;
proportionGT1Track = numObjectsGT1Track / totalNumObjects;
proportionGT2Track = numObjectsGT2Track / totalNumObjects;

meanFractionAreaFilledWithSpots = mean(FractionAreaFilledWithSpots, 'omitnan');
stdFractionAreaFilledWithSpots  = std(FractionAreaFilledWithSpots,  'omitnan');

meanSpotsPerObject = mean(NumSpots, 'omitnan');
stdSpotsPerObject  = std(NumSpots,  'omitnan');

ImageSummaryRow = table(imageName, string(dateStr), Tval, string(xyStr), ...
    metaInfo.StrainCode, metaInfo.Genotype, metaInfo.Hours, metaInfo.Replicate, metaInfo.Arabinose, ...
    metaInfo.Nspd, metaInfo.TimeStep_m, metaInfo.Pixel_um, ...
    totalNumObjects, numObjectsWithTracks, proportionWithTracks, ...
    numObjectsGT1Track, proportionGT1Track, ...
    numObjectsGT2Track, proportionGT2Track, ...
    totalAttachedSpots_All, numAttachedSpotsLT150, numAttachedSpotsGE150, ...
    meanDistToNearestDilCentroid_All, meanDistToNearestDilCentroid_LT150, ...
    meanFractionAreaFilledWithSpots, stdFractionAreaFilledWithSpots, ...
    meanSpotsPerObject, stdSpotsPerObject, ...
    'VariableNames', {'ImageName','Date','T','XY', ...
    'StrainCode','Genotype','Hours','Replicate','Arabinose','Nspd','TimeStep_m','Pixel_um', ...
    'TotalNumObjects','NumObjectsWithTracks','ProportionWithTracks', ...
    'NumObjectsGT1Track','ProportionGT1Track', ...
    'NumObjectsGT2Track','ProportionGT2Track', ...
    'TotalAttachedSpots_All','NumAttachedSpotsLT150','NumAttachedSpotsGE150', ...
    'MeanDistToNearestDilCentroid_All','MeanDistToNearestDilCentroid_LT150', ...
    'MeanFractionAreaFilledWithSpots','StdFractionAreaFilledWithSpots', ...
    'MeanSpotsPerObject','StdSpotsPerObject'});

result = struct;
% OBJECT_ID is dropped from the saved SpotsSummary: it's only meaningful
% within this one timepoint's own label map L, so keeping it in the
% multi-timepoint combined table would wrongly suggest object 5 at T=1
% and object 5 at T=8 are the same physical colony.
result.SpotsSummary = removevars(SpotsSummary, 'OBJECT_ID');
result.ObjectSummary = ObjectSummary;
result.ImageSummaryRow = ImageSummaryRow;

end

function metaInfo = getMetaInfo(meta, strainMap, dateStr, xyStr)
% Looks up this date's metadata row and decodes its per-XY-position
% column: an underscore-joined "StrainCode_Hours_Replicate_Arabinose[_Nspd]"
% string. Trailing parts are optional -- if fewer than 5 underscore-
% separated fields are present, the corresponding metaInfo fields are
% left as NaN. TimeStep_m/Pixel_um come from fixed per-date columns
% instead, not from this encoded string.

    row = find(meta.Date == dateStr, 1);
    cellVal = string(meta.(xyStr)(row));
    parts = split(strtrim(cellVal), "_");

    metaInfo = struct;
    metaInfo.StrainCode = "";
    metaInfo.Genotype = "";
    metaInfo.Hours = NaN;
    metaInfo.Replicate = NaN;
    metaInfo.Arabinose = NaN;
    metaInfo.Nspd = NaN;              % NEW
    metaInfo.TimeStep_m = meta.("Time step")(row);
    metaInfo.Pixel_um = meta.Pixel(row);

    if numel(parts) >= 1
        metaInfo.StrainCode = strtrim(parts(1));
        if isKey(strainMap, char(metaInfo.StrainCode))
            metaInfo.Genotype = string(strainMap(char(metaInfo.StrainCode)));
        end
    end

    if numel(parts) >= 2
        metaInfo.Hours = str2double(parts(2));
    end

    if numel(parts) >= 3
        metaInfo.Replicate = str2double(parts(3));
    end

    if numel(parts) >= 4
        metaInfo.Arabinose = str2double(parts(4));
    end

    if numel(parts) >= 5           % NEW
        metaInfo.Nspd = str2double(parts(5));
    end

end
