%% QuantifyPunctaCellsMicrocolonies.m
% Author: Hasreet Gill, Bassler Lab, HHMI/Princeton University
% Revised with Claude (Sonnet 5, Anthropic)

% Per-condition, per-replicate, per-field-of-view quantification of
% smFISH puncta (561/640 channels) relative to single cells (405 DAPI +
% 488 GFP) and microcolonies, for exported 4-channel max-projection
% TIFFs ("MAX_..._Export_..._<channel>.tif", one file per channel per
% field of view).
%
% For each experimental CONDITION (files grouped by stripping the
% trailing replicate-number suffix off their base filename) and each
% REPLICATE within it, every field of view (XY position) is:
%   1. Puncta-fit in the 561 and 640 channels independently (LoG blob
%      detection + 2D Gaussian fit per spot -- see extractpunctae).
%   2. Segmented into single-cell objects from DAPI + GFP together (see
%      the CELL_MASK_MODE toggle below -- this is the one step that
%      differs between the two experiment types this script supports),
%      then watershed-split to separate touching cells.
%   3. Segmented into microcolonies: large (>=3000 px) DAPI-positive
%      blobs that aren't already accounted for by the single-cell mask.
%   4. Every fitted puncta is assigned to whichever single cell and/or
%      microcolony (if any) it falls inside, and each cell/microcolony's
%      row records how many puncta (and their summed fit metrics)
%      landed in it, plus (for cells) distance to the nearest
%      microcolony and (for both) whether it came from DAPI, GFP, or a
%      watershed split.
%   5. QC overlay TIFFs (per-spot ROI/peak-call images, and per-FOV
%      green/non-green x split/non-split cell-outline overlays, and
%      microcolony-outline overlays) are appended one page per FOV.
% Per condition, all of this rolls up into one cell-level (Data_*.xlsx),
% one microcolony-level (DataMicro_*.xlsx), and per-channel puncta
% histogram/sigma spreadsheets, plus matching .mat structs.
%
% NOTE: redhistall/frhistall/areas/cellareas/microareas/redsigmas/
% frsigmas accumulate across every condition in the outer loop but are
% never written out or otherwise used after the loop finishes

clear
clc
close all

%% ================= CELL-MASK MODE (choose one before running) =================
% Two different rules exist for combining the DAPI (blue) single-cell
% mask with the GFP (green) single-cell mask into one "cells" mask, and
% they are NOT interchangeable -- pick whichever matches this dataset's
% experiment type. Only one runs; there is no dual-run option.
%
%   "Coculture": DAPI-first. A GFP-segmented object is added to the cell
%       mask only if its centroid does NOT already fall inside an
%       existing DAPI object -- i.e. it fills in cells DAPI missed, but
%       never overrides a DAPI object's own shape. greenOnlyMask marks
%       ONLY these DAPI-less additions.
%   "SpikeIn": GFP-first. A GFP-segmented object is always added, AND if
%       its centroid falls inside an existing DAPI object, that DAPI
%       object's shape is deleted and replaced by the GFP-derived shape
%       -- GFP segmentation takes priority over DAPI whenever both exist
%       for the same cell. greenOnlyMask marks EVERY GFP-derived pixel
%       inserted into the mask (matched or not), not just DAPI-less
%       additions -- so it excludes more area from the microcolony mask
%       (microcol(greenOnlyMask) = false, below) than "Coculture" mode
%       does for the same image. This mode also applies an extra
%       morphological close to the green mask before combining (sealing
%       small gaps the "Coculture" mode leaves alone).
CELL_MASK_MODE = "Coculture";   % "Coculture" or "SpikeIn"

%% ================= SETUP =================

path = uigetdir();
cd(path)

outDir = fullfile(path, "Output");
if ~exist(outDir, 'dir')
    mkdir(outDir);
end

ims = dir("MAX*.tif"); % Collect all images
ims = string({ims.name})';
baseAll = extractBefore(ims, "_Export");

% group all replicates of the same condition together
groupNames = regexprep(baseAll, '_(\d+)(?:_\d+)?$', '');
groupuni = unique(groupNames, 'stable')';

redhistall = zeros(0,4);
frhistall = zeros(0,4);
areas = zeros(0,1);
cellareas = zeros(0,1);
microareas = zeros(0,1);
redsigmas = zeros(0,1);
frsigmas = zeros(0,1);

for i = 1:length(groupuni)

    data = struct([]);
    datam = struct([]);

    redhist = table('Size',[0 6], ...
        'VariableTypes', {'double','string','double','double','double','double'}, ...
        'VariableNames', {'Replicate','XYPos','PeakIntensity','GaussianArea','BgSubSum','PeakIntensityFixedSigma'});

    frhist = table('Size',[0 6], ...
        'VariableTypes', {'double','string','double','double','double','double'}, ...
        'VariableNames', {'Replicate','XYPos','PeakIntensity','GaussianArea','BgSubSum','PeakIntensityFixedSigma'});

    redsig = table('Size',[0 3], ...
        'VariableTypes', {'double','string','double'}, ...
        'VariableNames', {'Replicate','XYPos','Sigma'});

    frsig = table('Size',[0 3], ...
        'VariableTypes', {'double','string','double'}, ...
        'VariableNames', {'Replicate','XYPos','Sigma'});

    thisGroup = groupuni(i);
    groupMask = groupNames == thisGroup;
    imsGroup = ims(groupMask);
    baseGroup = extractBefore(imsGroup, "_Export");
    baseUni = unique(baseGroup, 'stable')';

    % get replicate number
    repnums = nan(length(baseUni),1);
    for b = 1:length(baseUni)
        tok = regexp(baseUni(b), '_(\d+)(?:_\d+)?$', 'tokens', 'once');
        repnums(b) = str2double(tok{1});
    end
    reps = unique(repnums(~isnan(repnums)))';

    for r = 1:length(reps)

        repnum = reps(r);
        repBases = baseUni(repnums == repnum);
        repMask = ismember(baseGroup, repBases);
        xyc = imsGroup(repMask);
        xy = unique(extractBetween(xyc, "MAX_", "_2.8X"));
        repBaseName = thisGroup + "_Rep" + string(repnum);

        for j = 1:length(xy)

            % load each color
            blcolor = '405';
            grcolor = '488';
            redcolor = '561';
            frcolor = '640';

            blname = xyc(contains(xyc, xy(j)) & contains(xyc, blcolor));
            bl = imread(blname);
            grname = xyc(contains(xyc, xy(j)) & contains(xyc, grcolor));
            gr = imread(grname);
            redname = xyc(contains(xyc, xy(j)) & contains(xyc, redcolor));
            red = imread(redname);
            frname = xyc(contains(xyc, xy(j)) & contains(xyc, frcolor));
            fr = imread(frname);

            % measure smFISH punctae
            [redfits_xy, redfits_mat, sigmared] = extractpunctae(repBaseName, redcolor, red, outDir);
            [frfits_xy, frfits_mat, sigmafr] = extractpunctae(repBaseName, frcolor, fr, outDir);

            % store puncta outputs with replicate + XY bookkeeping
            if ~isempty(redfits_mat)
                Tred = array2table(redfits_mat, 'VariableNames', ...
                    {'PeakIntensity','GaussianArea','BgSubSum','PeakIntensityFixedSigma'});
                Tred.Replicate = repmat(repnum, height(Tred), 1);
                Tred.XYPos = repmat(string(xy(j)), height(Tred), 1);
                Tred = movevars(Tred, {'Replicate','XYPos'}, 'Before', 1);

                TsigRed = table(repmat(repnum, numel(sigmared), 1), ...
                                repmat(string(xy(j)), numel(sigmared), 1), ...
                                sigmared(:), ...
                                'VariableNames', {'Replicate','XYPos','Sigma'});

                redhist = [redhist; Tred];
                redsig = [redsig; TsigRed];
            end

            if ~isempty(frfits_mat)
                Tfr = array2table(frfits_mat, 'VariableNames', ...
                    {'PeakIntensity','GaussianArea','BgSubSum','PeakIntensityFixedSigma'});
                Tfr.Replicate = repmat(repnum, height(Tfr), 1);
                Tfr.XYPos = repmat(string(xy(j)), height(Tfr), 1);
                Tfr = movevars(Tfr, {'Replicate','XYPos'}, 'Before', 1);

                TsigFr = table(repmat(repnum, numel(sigmafr), 1), ...
                               repmat(string(xy(j)), numel(sigmafr), 1), ...
                               sigmafr(:), ...
                               'VariableNames', {'Replicate','XYPos','Sigma'});

                frhist = [frhist; Tfr];
                frsig = [frsig; TsigFr];
            end

            % DAPI mask: smooth, adaptive-threshold, clear anything
            % touching the image border, close small gaps, then keep
            % only single-cell-sized objects (150-2000 px).
            smbl = imgaussfilt(bl, 6);
            bin = imbinarize(smbl, "adaptive", Sensitivity = 0.45);
            bin = imclearborder(bin);
            se = strel("disk", 10);
            bin = imclose(bin, se);
            binbl = bwareafilt(bin, [150 2000]);

            % GFP mask: lighter smoothing (GFP objects are typically
            % smaller/sharper than the DAPI signal), adaptive threshold,
            % border-clear, same single-cell size range.
            smgr = imgaussfilt(gr, 2);
            bingr = imbinarize(gr, "adaptive", Sensitivity = 0.30);
            bingr = imclearborder(bingr);
            bingr = bwareafilt(bingr, [150 2000]);

            if CELL_MASK_MODE == "SpikeIn"
                % SpikeIn mode only: seal small gaps in the GFP mask with
                % the same structuring element used for the DAPI mask
                % above, before it's combined with DAPI below.
                bingr = imclose(bingr, se);
            end

            % Build cell mask by combining DAPI and GFP objects.
            % ccgr/grprops and ccbl/blprops (below) are shared by both
            % CELL_MASK_MODE branches; only the combination RULE differs.
            ccgr = bwconncomp(bingr);
            grprops = regionprops(ccgr, 'Centroid', 'PixelIdxList');
            ccbl = bwconncomp(binbl);
            blprops = regionprops(ccbl, 'PixelIdxList');

            cells = binbl;
            greenOnlyMask = false(size(binbl));

            switch CELL_MASK_MODE

                case "Coculture"
                    % DAPI-first: only add a GFP object if it has NO
                    % matching DAPI object (by centroid-inside-DAPI-
                    % object test). Matched GFP objects are ignored --
                    % their corresponding DAPI object's own shape is
                    % kept as that cell's boundary.
                    for ii = 1:numel(grprops)
                        cgr = round(grprops(ii).Centroid);
                        x = cgr(1);
                        y = cgr(2);
                        hasMatch = false;
                        if x >= 1 && x <= size(binbl,2) && y >= 1 && y <= size(binbl,1)
                            pix = sub2ind(size(binbl), y, x);
                            for jj = 1:numel(blprops)
                                if any(blprops(jj).PixelIdxList == pix)
                                    hasMatch = true;
                                    break
                                end
                            end
                        end
                        if ~hasMatch
                            cells(grprops(ii).PixelIdxList) = true;
                            greenOnlyMask(grprops(ii).PixelIdxList) = true;
                        end
                    end

                case "SpikeIn"
                    % GFP-first: always add the GFP object; if it also
                    % matches a DAPI object, delete that DAPI object's
                    % shape first so the GFP shape replaces it instead of
                    % just supplementing it.
                    for ii = 1:numel(grprops)

                        cgr = round(grprops(ii).Centroid);
                        x = cgr(1);
                        y = cgr(2);

                        hasMatch = false;
                        matchIdx = NaN;

                        if x >= 1 && x <= size(binbl,2) && y >= 1 && y <= size(binbl,1)

                            pix = sub2ind(size(binbl), y, x);

                            for jj = 1:numel(blprops)
                                if any(blprops(jj).PixelIdxList == pix)
                                    hasMatch = true;
                                    matchIdx = jj;
                                    break
                                end
                            end

                        end

                        if hasMatch
                            % Remove the DAPI-derived mask for this green cell.
                            cells(blprops(matchIdx).PixelIdxList) = false;
                        end

                        % Add the green-derived mask for this green cell.
                        cells(grprops(ii).PixelIdxList) = true;
                        greenOnlyMask(grprops(ii).PixelIdxList) = true;

                    end

                otherwise
                    error("Unknown CELL_MASK_MODE: %s (expected \"Coculture\" or \"SpikeIn\")", CELL_MASK_MODE);
            end

            % Watershed-split the combined mask to separate any
            % touching/merged cells; ridgec records exactly which pixels
            % were the dividing lines, used below to flag which final
            % objects had to be split from a neighbor (WasSplit).
            Dc = bwdist(~cells);
            Dc = imhmax(Dc, 1);
            Lc = watershed(-Dc);
            ridgec = (Lc == 0);
            cells(Lc == 0) = 0;

            cc = bwconncomp(cells);
            L = labelmatrix(cc);
            cellprops = regionprops(cc, smgr, "MeanIntensity", 'PixelList','Centroid', ...
                'MajorAxisLength', 'MinorAxisLength','Circularity','Area');
            cellctrs = cat(1, cellprops.Centroid);
            cellar = cat(1, cellprops.MajorAxisLength)./cat(1, cellprops.MinorAxisLength);
            cella = cat(1, cellprops.Area);

            % Microcolonies: large (>=3000 px) DAPI-positive blobs that
            % aren't part of the single-cell (green-associated) mask --
            % i.e. multi-cell clumps too big to have survived binbl's
            % single-cell size filter in the first place.
            microcol = bin;
            microcol(greenOnlyMask) = false;
            microcol = bwareafilt(microcol, [3000 100000000]);
            ccm = bwconncomp(microcol);
            Lm = labelmatrix(ccm);
            microprops = regionprops(ccm, smgr, "MeanIntensity", 'Centroid','PixelList', 'Area');
            microctrs = cat(1, microprops.Centroid);
            microa = cat(1, microprops.Area);

            % Assign every fitted puncta to whichever single cell (L)
            % and/or microcolony (Lm) it physically falls inside; 0 in
            % either map means "not inside any object of that kind".
            objectID_red = [];
            objectIDm_red = [];
            objectID_fr = [];
            objectIDm_fr = [];

            for ii = 1:size(redfits_xy,1)
                xf = round(redfits_xy(ii,1));
                yf = round(redfits_xy(ii,2));
                objectID_red(ii) = L(yf, xf);
                objectIDm_red(ii) = Lm(yf, xf);
            end
            for ii = 1:size(frfits_xy,1)
                xf = round(frfits_xy(ii,1));
                yf = round(frfits_xy(ii,2));
                objectID_fr(ii) = L(yf, xf);
                objectIDm_fr(ii) = Lm(yf, xf);
            end

            if numel(cellprops) == 0
                data = data;
            else
                for k = 1:numel(cellprops)
                    mindist = [];
                    pix = cc.PixelIdxList{k};
                    thisMask = false(size(cells));
                    thisMask(pix) = true;

                    cellprops(k).Replicate = repnum;
                    cellprops(k).XYPos = string(xy(j));
                    cellprops(k).AspectRatio = cellar(k);
                    cellprops(k).FromBinbl = any(binbl(pix));
                    cellprops(k).FromBingr = any(bingr(pix));
                    cellprops(k).WasSplit = any(imdilate(thisMask, strel("disk",1)) & ridgec, 'all');
                    cellprops(k).SumRed = sum(redfits_mat(find(objectID_red==k),:),1);
                    cellprops(k).NRed = length(find(objectID_red==k));
                    cellprops(k).SumFarRed = sum(frfits_mat(find(objectID_fr==k),:),1);
                    cellprops(k).NFarRed = length(find(objectID_fr==k));
                    if numel(microprops) == 0
                        cellprops(k).MeanDistMicro = 0;
                        cellprops(k).MinDistMicro = 0;
                        cellprops(k).NearestMicroInt = 0;
                    else
                        for l = 1:numel(microprops)
                            micropix = cat(1, microprops(l).PixelList);
                            mindist(l) = min(pdist2(cellctrs(k,:), micropix));
                        end
                        [minVal, minIdx] = min(mindist);
                        cellprops(k).MeanDistMicro = mean(mindist);
                        cellprops(k).MinDistMicro = minVal;
                        cellprops(k).NearestMicroInt = microprops(minIdx).MeanIntensity;
                    end
                end
                data = [data; cellprops(:)];
            end

            % QC overlay: split every final cell object into green vs.
            % non-green (by overlap with the original bingr mask) and,
            % within each, split-from-a-neighbor vs. not (WasSplit).
            grIdx = [];
            ngrIdx = [];
            grIdx_nosplit = [];
            ngrIdx_nosplit = [];

            for k = 1:cc.NumObjects
                pix = cc.PixelIdxList{k};
                isGreen = any(bingr(pix));
                isSplit = cellprops(k).WasSplit;

                if isGreen
                    grIdx(end+1) = k;
                    if ~isSplit
                        grIdx_nosplit(end+1) = k;
                    end
                else
                    ngrIdx(end+1) = k;
                    if ~isSplit
                        ngrIdx_nosplit(end+1) = k;
                    end
                end
            end

            ccgr = cc;
            ccgr.PixelIdxList = ccgr.PixelIdxList(grIdx');
            ccgr.NumObjects = numel(grIdx);
            grmask = boundarymask(labelmatrix(ccgr));

            ccngr = cc;
            ccngr.PixelIdxList = ccngr.PixelIdxList(ngrIdx');
            ccngr.NumObjects = numel(ngrIdx);
            ngrmask = boundarymask(labelmatrix(ccngr));

            ccgr_ns = cc;
            ccgr_ns.PixelIdxList = ccgr_ns.PixelIdxList(grIdx_nosplit');
            ccgr_ns.NumObjects = numel(grIdx_nosplit);
            grmask_ns = boundarymask(labelmatrix(ccgr_ns));

            ccngr_ns = cc;
            ccngr_ns.PixelIdxList = ccngr_ns.PixelIdxList(ngrIdx_nosplit');
            ccngr_ns.NumObjects = numel(ngrIdx_nosplit);
            ngrmask_ns = boundarymask(labelmatrix(ccngr_ns));

            % Composite background: blue = DAPI, green = GFP, red = blank.
            Bch = adapthisteq(mat2gray(smbl));
            G = adapthisteq(mat2gray(smgr));
            R = zeros(size(Bch));
            C = cat(3, R, G, Bch);

            % Left panel = non-green cells: white outline for every
            % non-green cell, then red redrawn on top for just the
            % NOT-split ones -- so red = clean single cell, white (left
            % showing through) = was split from a neighbor.
            C_left = labeloverlay(C, ngrmask, 'Colormap', [1 1 1], 'Transparency', 0);
            C_left = labeloverlay(C_left, ngrmask_ns, 'Colormap', [1 0 0], 'Transparency', 0);
            % Right panel = same red/white split-vs-clean convention, for
            % green cells instead.
            C_right = labeloverlay(C, grmask, 'Colormap', [1 1 1], 'Transparency', 0);
            C_right = labeloverlay(C_right, grmask_ns, 'Colormap', [1 0 0], 'Transparency', 0);

            C_combined = cat(2, C_left, C_right);

            imwrite(im2uint8(C_combined), ...
                fullfile(outDir, strcat('Overlay_', repBaseName, '.tif')), ...
                'WriteMode','append')

            if numel(microprops) == 0
                datam = datam;
            else
                for g = 1:numel(microprops)
                    microprops(g).Replicate = repnum;
                    microprops(g).XYPos = string(xy(j));
                    microprops(g).SumRed = sum(redfits_mat(find(objectIDm_red==g),:),1);
                    microprops(g).NRed = length(find(objectIDm_red==g));
                    microprops(g).SumFarRed = sum(frfits_mat(find(objectIDm_fr==g),:),1);
                    microprops(g).NFarRed = length(find(objectIDm_fr==g));
                end
                datam = [datam; microprops(:)];
            end

            % microcolonies overlay
            maskm = boundarymask(Lm);
            microOverlay = labeloverlay(mat2gray(smbl),maskm,'Transparency',0);
            imwrite(im2uint8(microOverlay), ...
                fullfile(outDir, strcat('Microcolonies_', repBaseName, '.tif')), ...
                'WriteMode','append')
        end
    end

    redhistall = [redhistall; table2array(redhist(:,3:end))];
    frhistall = [frhistall; table2array(frhist(:,3:end))];
    areas = [areas;cella;microa];
    cellareas = [cellareas;cella];
    microareas = [microareas;microa];
    redsigmas = [redsigmas; redsig.Sigma];
    frsigmas = [frsigmas; frsig.Sigma];

    % save data
    save(fullfile(outDir, strcat('Struct_', thisGroup, '.mat')), 'data');
    save(fullfile(outDir, strcat('Struct_Micro_', thisGroup, '.mat')), 'datam');

    T1 = struct2table(rmfield(data, 'PixelList'));
    T2 = struct2table(rmfield(datam, 'PixelList'));

    writetable(T1, fullfile(outDir, strcat('Data_', thisGroup, '.xlsx')));
    writetable(T2, fullfile(outDir, strcat('DataMicro_', thisGroup, '.xlsx')));
    writetable(redhist, fullfile(outDir, strcat('RedHistograms_', thisGroup, '.xlsx')));
    writetable(frhist, fullfile(outDir, strcat('FarRedHistograms_', thisGroup, '.xlsx')));
    writetable(redsig, fullfile(outDir, strcat('RedSigmas_', thisGroup, '.xlsx')));
    writetable(frsig, fullfile(outDir, strcat('FarRedSigmas_', thisGroup, '.xlsx')));

end

function [fit_xy, fits, sigma] = extractpunctae(imname, color, im, outDir)
% Detects diffraction-limited smFISH puncta in one channel image and
% fits each with a 2D Gaussian for sub-pixel position and intensity.
%
%   1. A rough adaptive threshold (bin) marks candidate signal regions,
%      used only to estimate the local background level (the mean
%      intensity of everything OUTSIDE bin) and to reject noise-only LoG
%      peaks below -- it is not the final puncta segmentation itself.
%   2. The image is flat background-subtracted using that estimate, then
%      run through a Laplacian-of-Gaussian blob filter (sigma_seg = 2,
%      matched to the expected puncta size); every local maximum of the
%      LoG response that also falls inside the rough signal mask is kept
%      as a candidate puncta.
%   3. Each candidate gets its own Gaussian-fit ROI, sized to the
%      smaller of half the distance to its nearest neighboring puncta
%      (so fit windows don't overlap) and half the distance to the image
%      edge (so the window never runs off the image), capped at 20 px.
%   4. Two 2D Gaussian fits are made per ROI: a free-sigma fit (its own
%      best-fit width) and a second fit with sigma fixed at 1.4 px (a
%      more robust/comparable amplitude estimate assuming every puncta
%      shares the same point-spread-function width) -- returned as two
%      separate intensity metrics.
%
% Returns fit_xy (free-fit and fixed-fit sub-pixel [x,y] per puncta),
% fits (peak intensity, Gaussian-integrated area, background-subtracted
% raw pixel sum, and fixed-sigma peak intensity -- 4 columns), and sigma
% (the free fit's fitted width per puncta). Also saves two QC TIFFs per
% call: peak-call locations overlaid on the image, and each puncta's
% fit-ROI circle + raw detected position.

        bin = imbinarize(im, "adaptive", Sensitivity = 0.5);
        bin = bwareaopen(bin,5);
        im_inv = im2double(im).* double(~bin);
        im_inv_nan = im_inv;
        im_inv_nan(im_inv_nan == 0) = NaN;
        thresh = mean(im_inv_nan, "all", "omitnan");

        im_sub = im2double(im) - thresh;
        im_sub(im_sub < 0) = 0;

        % laplacian of Gaussian
        sigma_seg = 2;
        hsize = 2*ceil(3*sigma_seg)+1;
        h = fspecial('log', hsize, sigma_seg);
        response = -imfilter(im_sub, h, 'replicate');
        maxim = imregionalmax(response);
        peaks = maxim & bin;
        [y,x] = find(peaks==1);

        % initialize outputs (handles zero-peak case)
        fit_xy = zeros(0,4);
        fits = zeros(0,4);
        sigma = zeros(0,1);

        figure("visible", "off")
        peaksim = imdilate(peaks,strel("disk",1));
        imshowpair(adapthisteq(mat2gray(im_sub)),peaksim);
        F1 = getframe;
        imwrite(F1.cdata, fullfile(outDir, strcat('Peakcalls_', imname, '_', color, '.tif')), ...
            'WriteMode','append')

        % ROI sizes for Gaussians
        D = pdist2([x y], [x y]);
        D(D == 0) = inf;
        nearestDist = min(D,[],2);
        halfSizes = floor(nearestDist/2);
        halfSizes = min(halfSizes, 20);     % cap at 20 if needed
        [h, w] = size(im_sub);
        distToEdge = min([x, w - x, y, h - y], [], 2);
        halfdistToEdge = floor(distToEdge/2);
        theta = linspace(0, 2*pi, 200);
        for i = 1:length(x)
            r = min(halfSizes(i),halfdistToEdge(i));
            x0 = x(i);
            y0 = y(i);
            x1 = x0-r;
            x2 = x0+r;
            y1 = y0-r;
            y2 = y0+r;

            ROI = im_sub(y1:y2, x1:x2);
            gauss2D = @(p,xy) p(1)*exp(-((xy(:,1)-p(2)).^2 + ...
                              (xy(:,2)-p(3)).^2)/(2*p(4)^2)) + p(5);
            [X,Y] = meshgrid(1:(2*r+1), 1:(2*r+1));
            xy = [X(:), Y(:)];
            Z = ROI(:);
            p0 = [max(Z)-median(Z), r+1, r+1, 1, median(Z)];
            options = optimoptions('lsqcurvefit','Display','off');
            lb = [0, 1, 1, 0.5, 0];
            ub = [max(Z)*2, 2*r+1, 2*r+1, 4, max(Z)];
            p = lsqcurvefit(gauss2D, p0, xy, Z, lb, ub, options);

            sigma_fixed = 1.4;
            gauss2D_fixed = @(p,xy) p(1)*exp(-((xy(:,1)-p(2)).^2 + (xy(:,2)-p(3)).^2) / (2*sigma_fixed^2)) + p(4);
            p0_fixed = [max(Z)-median(Z), r+1, r+1, median(Z)];
            lb_fixed = [0, 1, 1, 0];
            ub_fixed = [max(Z)*2, 2*r+1, 2*r+1, max(Z)];
            p_fixed = lsqcurvefit(gauss2D_fixed, p0_fixed, xy, Z, lb_fixed, ub_fixed, options);

            x_fit = x1 - 1 + p(2);
            y_fit = y1 - 1 + p(3);
            x_fit_fixed = x1 - 1 + p_fixed(2);
            y_fit_fixed = y1 - 1 + p_fixed(3);
            fit_xy(i,:) = [x_fit, y_fit, x_fit_fixed, y_fit_fixed];
            fits(i,1) = p(1); % just peak intensity
            fits(i,2) = p(1) * 2*pi*p(4).^2; % area under gaussian fit
            fits(i,3) = sum(ROI(:) - p(5)); % background-subtracted sum
            fits(i,4) = p_fixed(1); % fixed sigma peak intensity
            sigma(i) = p(4);
        end

        base = adapthisteq(mat2gray(im_sub));
        rgb2 = repmat(base, [1 1 3]);
        circles = zeros(length(x), 3);  % [x_center, y_center, radius]
        for i = 1:length(x)
            r = min(halfSizes(i), halfdistToEdge(i));
            circles(i,:) = [x(i), y(i), r];
        end

        rgb2 = insertShape(rgb2, 'Circle', circles, ...
        'Color', 'red', 'LineWidth', 1);
        centers_xy = [x(:), y(:)];
        rgb2 = insertMarker(rgb2, centers_xy, '+', ...
        'Color', 'yellow', 'Size', 1);
        imwrite(rgb2, fullfile(outDir, strcat('ROIBoundaries_', imname, '_', color, '.tif')), ...
        'WriteMode','append');
end
