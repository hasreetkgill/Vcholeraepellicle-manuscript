%% Shift_Stack_ManualPlusAutoCorrection.m
% Author: Hasreet Gill, Bassler Lab, HHMI/Princeton University
% Developed/revised with Claude (Sonnet 5, Anthropic) - AI-assisted code review and comments.

% Builds a drift-corrected image stack from a folder of single-frame
% "Shift_1_<n>_Export...tif" snapshots (the bead images saved by the
% NIS-Elements XY-correction macro). Three versions of the stack are
% produced, at increasing levels of correction:
%   1. Uncorrected  -- frames as captured, no shift applied.
%   2. ManualOnly   -- each frame shifted by an X/Y value the user types
%                      in via a dialog (for drift too large/irregular for
%                      the automatic bead-centroid matching below to
%                      handle reliably).
%   3. Corrected    -- ManualOnly shift PLUS an automatic residual shift
%                      computed by matching bead centroids between frame 1
%                      (reference) and every other frame (same
%                      median-nearest-neighbor-shift idea as
%                      ComputeMedianShift in the acquisition macro).
% For each version, both the grayscale stack and a binarized/bead-mask
% stack are saved as multi-page TIFFs and as .mat files. A CSV table of
% every frame's manual shift, residual shift, and total (manual +
% residual) shift is also written, so the correction can be audited or
% reapplied elsewhere.

clear; clc; close all

%% ===================== Locate and sort input frames =====================

path = uigetdir;
cd(path);

% Only take files that look like "Shift_1_<number>_Export....tif" -- this
% excludes the "Shift_2_..." (post-correction) frames and any other stray
% tif in the folder.
bdsfiles = dir("Shift_1*.tif");
names = {bdsfiles.name};
keep = ~cellfun(@isempty, ...
    regexp(names, '^Shift_1_\d+_Export.*\.tif$', 'once'));
bdsfiles = bdsfiles(keep);
names = {bdsfiles.name};

% dir() does not sort numerically (Shift_1_10 would sort before
% Shift_1_2 by name), so extract the frame number and sort on that.
nums = cellfun(@(x) ...
    str2double(regexp(x, '(?<=Shift_1_)\d+', 'match', 'once')), ...
    names);
[~, idx] = sort(nums);
bdsfiles = bdsfiles(idx);

im1 = imread(bdsfiles(1).name);
dims = size(im1);
nF = length(bdsfiles);

%% ===================== Output filenames =====================

outTif_uncorr    = 'Shift_Stack_Uncorrected.tif';
outTifbin_uncorr = 'Shift_Stack_Binary_Uncorrected.tif';
outMat_uncorr    = 'Shift_Stack_Uncorrected.mat';
outMatbin_uncorr = 'Shift_Stack_Binary_Uncorrected.mat';

outTif_manual    = 'Shift_Stack_ManualOnly.tif';
outTifbin_manual = 'Shift_Stack_Binary_ManualOnly.tif';
outMat_manual    = 'Shift_Stack_ManualOnly.mat';
outMatbin_manual = 'Shift_Stack_Binary_ManualOnly.mat';

outTif_corr      = 'Shift_Stack_Corrected.tif';
outTifbin_corr   = 'Shift_Stack_Binary_Corrected.tif';
outMat_corr      = 'Shift_Stack_Corrected.mat';
outMatbin_corr   = 'Shift_Stack_Binary_Corrected.mat';

outCsvShifts     = 'Shift_Stack_TotalShifts.csv';

% Fresh run every time -- remove any output from a previous run in this
% folder so partially-overwritten stacks (e.g. fewer frames than before)
% can't linger and get mistaken for current results.
del = dir('Shift_Stack*');
if ~isempty(del)
    delete('Shift_Stack*');
end

% Max pixel distance for a bead in frame 1 to be paired with the "same"
% bead in another frame when computing the automatic residual shift
% (same role as MAX_DIST in the acquisition macro's ComputeMedianShift).
maxdist = 10;

allIms = cell(nF,1);
allBins = cell(nF,1);

bds_uncorr = zeros(dims(1), dims(2), nF, 'uint16');
bdsbin_uncorr = false(dims(1), dims(2), nF);

%% ===================== STAGE 1: load frames, save uncorrected stack =====================
% Also builds the binary bead mask for every frame here (adaptive
% threshold + drop objects smaller than 5 px) so it doesn't need to be
% recomputed later -- allIms/allBins are reused for the manual and
% corrected stages below.

for f = 1:nF

    imname = bdsfiles(f).name;
    im = imread(imname);
    allIms{f} = im;

    bin = imbinarize(im, "adaptive", Sensitivity=0.4);
    bin = bwareafilt(bin, [5 Inf]);
    allBins{f} = bin;

    bds_uncorr(:,:,f) = im;
    bdsbin_uncorr(:,:,f) = bin;

    binsave = uint8(bin * 255);

    if f == 1
        imwrite(binsave, outTifbin_uncorr, 'tif', 'Compression', 'none');
        imwrite(im, outTif_uncorr, 'tif', 'Compression', 'none');
    else
        imwrite(binsave, outTifbin_uncorr, 'tif', 'WriteMode', 'append', 'Compression', 'none');
        imwrite(im, outTif_uncorr, 'tif', 'WriteMode', 'append', 'Compression', 'none');
    end

    fprintf('Saved uncorrected frame %d\n', f);
end

save(outMat_uncorr, 'bds_uncorr');
save(outMatbin_uncorr, 'bdsbin_uncorr');

fprintf('Saved uncorrected stack.\n');

%% ===================== STAGE 2: manual correction, one dialog per frame =====================
% Drift tends to accumulate frame-to-frame, so each dialog is
% pre-populated with the PREVIOUS frame's shift (prevX/prevY) rather than
% 0, letting the user just confirm or nudge it instead of re-typing the
% cumulative value every time.

manualX = zeros(nF,1);
manualY = zeros(nF,1);

prevX = 0;
prevY = 0;

for f = 1:nF

    answer = inputdlg( ...
        {sprintf('Frame %d X shift', f), sprintf('Frame %d Y shift', f)}, ...
        sprintf('Manual correction: frame %d of %d', f, nF), ...
        [1 25; 1 25], ...
        {num2str(prevX), num2str(prevY)});

    % User hit Cancel: stop asking and leave every remaining frame at its
    % zero-initialized shift rather than erroring out or guessing.
    if isempty(answer)
        fprintf('Manual correction cancelled at frame %d. Remaining frames kept as 0 shift.\n', f);
        break;
    end

    newX = str2double(answer{1});
    newY = str2double(answer{2});

    if isnan(newX) || isnan(newY)
        errordlg('X and Y must be numeric. This frame will stay at 0 shift.', 'Invalid shift');
        newX = 0;
        newY = 0;
    end

    manualX(f) = newX;
    manualY(f) = newY;

    prevX = newX;
    prevY = newY;

    fprintf('Manual frame %d: X = %.3f, Y = %.3f\n', f, manualX(f), manualY(f));
end

%% ===================== Apply manual shifts, save manual-only stack =====================

bds_manual = zeros(dims(1), dims(2), nF, 'uint16');
bdsbin_manual = false(dims(1), dims(2), nF);

for f = 1:nF

    im = allIms{f};
    bin = allBins{f};

    imManual  = imtranslate(im,  [manualX(f), manualY(f)], 'OutputView', 'same', 'FillValues', 0);
    binManual = imtranslate(bin, [manualX(f), manualY(f)], 'OutputView', 'same', 'FillValues', 0);

    bds_manual(:,:,f) = imManual;
    bdsbin_manual(:,:,f) = binManual;

    binsave = uint8(binManual * 255);

    if f == 1
        imwrite(binsave, outTifbin_manual, 'tif', 'Compression', 'none');
        imwrite(imManual, outTif_manual, 'tif', 'Compression', 'none');
    else
        imwrite(binsave, outTifbin_manual, 'tif', 'WriteMode', 'append', 'Compression', 'none');
        imwrite(imManual, outTif_manual, 'tif', 'WriteMode', 'append', 'Compression', 'none');
    end

    fprintf('Saved manual-only frame %d\n', f);
end

save(outMat_manual, 'bds_manual', 'manualX', 'manualY');
save(outMatbin_manual, 'bdsbin_manual', 'manualX', 'manualY');

fprintf('Saved manual-only stack.\n');

%% ===================== STAGE 3: automatic residual correction =====================
% After the manual shift, small leftover drift may remain. This matches
% bead centroids in the manually-corrected binary mask of frame 1
% (reference) against every other frame and takes the median
% nearest-neighbor displacement (getMedianShift, below) as that frame's
% residual shift -- the same pairwise-distance/median approach as
% ComputeMedianShift in the acquisition macro, just run once per frame
% against a fixed reference instead of consecutive-frame-to-frame.

resX = zeros(nF,1);
resY = zeros(nF,1);

refProps = regionprops(bdsbin_manual(:,:,1), 'Centroid');

if isempty(refProps)
    refCtrs = [];
else
    refCtrs = cat(1, refProps.Centroid);
end

for f = 1:nF

    % Frame 1 is the reference itself, so its residual shift is 0 by
    % definition.
    if f == 1
        resX(f) = 0;
        resY(f) = 0;
        continue;
    end

    props = regionprops(bdsbin_manual(:,:,f), 'Centroid');

    if isempty(props)
        ctrs = [];
    else
        ctrs = cat(1, props.Centroid);
    end

    if isempty(refCtrs) || isempty(ctrs)
        medianX = 0;
        medianY = 0;
    else
        [~, ~, medianX, medianY] = getMedianShift(refCtrs, ctrs, maxdist);
    end

    resX(f) = medianX;
    resY(f) = medianY;

    fprintf('Auto residual frame %d: X = %.3f, Y = %.3f\n', f, resX(f), resY(f));
end

% Combined shift actually needed to bring each ORIGINAL frame in line
% with the reference: manual entry plus whatever residual drift the
% centroid matching found on top of it.
totalX = manualX + resX;
totalY = manualY + resY;

T = table((1:nF)', manualX, manualY, resX, resY, totalX, totalY, ...
    'VariableNames', {'Frame', 'ManualX', 'ManualY', 'ResidualX', 'ResidualY', 'TotalX', 'TotalY'});

writetable(T, outCsvShifts);

fprintf('Saved shift table:\n%s\n', outCsvShifts);

%% ===================== Apply total shift, save final corrected stack =====================
% Shifts the ORIGINAL (allIms/allBins), not the already-manually-shifted
% image, by the combined totalX/totalY in one step -- avoids compounding
% two separate imtranslate interpolation passes into the final result.

bds_corr = zeros(dims(1), dims(2), nF, 'uint16');
bdsbin_corr = false(dims(1), dims(2), nF);

for f = 1:nF

    im = allIms{f};
    bin = allBins{f};

    imShift  = imtranslate(im,  [totalX(f), totalY(f)], 'OutputView', 'same', 'FillValues', 0);
    binShift = imtranslate(bin, [totalX(f), totalY(f)], 'OutputView', 'same', 'FillValues', 0);

    bds_corr(:,:,f) = imShift;
    bdsbin_corr(:,:,f) = binShift;

    binsave = uint8(binShift * 255);

    if f == 1
        imwrite(binsave, outTifbin_corr, 'tif', 'Compression', 'none');
        imwrite(imShift, outTif_corr, 'tif', 'Compression', 'none');
    else
        imwrite(binsave, outTifbin_corr, 'tif', 'WriteMode', 'append', 'Compression', 'none');
        imwrite(imShift, outTif_corr, 'tif', 'WriteMode', 'append', 'Compression', 'none');
    end

    fprintf('Saved corrected frame %d with total X = %.3f, total Y = %.3f\n', ...
        f, totalX(f), totalY(f));
end

save(outMat_corr, 'bds_corr', 'manualX', 'manualY', 'resX', 'resY', 'totalX', 'totalY');
save(outMatbin_corr, 'bdsbin_corr', 'manualX', 'manualY', 'resX', 'resY', 'totalX', 'totalY');

fprintf('Saved corrected stack.\n');


function [xShifts, yShifts, medianX, medianY] = getMedianShift(prevCtrs, currCtrs, MAX_DIST)
% Pairs every centroid in prevCtrs with every centroid in currCtrs closer
% than MAX_DIST apart, collects their per-pair X/Y displacements, and
% returns the median of each as the frame's overall shift estimate (a
% MATLAB port of the ComputeMedianShift logic used in the acquisition
% macro). Returns 0/0 if no pair falls within MAX_DIST.

xShifts = [];
yShifts = [];

for i = 1:size(prevCtrs,1)
    for j = 1:size(currCtrs,1)

        dx = prevCtrs(i,1) - currCtrs(j,1);
        dy = prevCtrs(i,2) - currCtrs(j,2);

        d = sqrt(dx^2 + dy^2);

        if d < MAX_DIST
            xShifts(end+1) = dx; %#ok<AGROW>
            yShifts(end+1) = dy; %#ok<AGROW>
        end
    end
end

if isempty(xShifts)
    medianX = 0;
    medianY = 0;
else
    medianX = median(xShifts);
    medianY = median(yShifts);
end

end
