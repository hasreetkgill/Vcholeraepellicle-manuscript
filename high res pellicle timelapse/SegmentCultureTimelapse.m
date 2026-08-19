function SegmentCultureTimelapse(path)
% Author: Hasreet Gill, Bassler Lab, HHMI/Princeton University
% Revised with Claude (Sonnet 5, Anthropic)

% Cluster batch driver (for the "della" SLURM cluster) that segments and
% measures every Z-stack timelapse image under PATH.
%
% PATH is expected to contain one subfolder per imaging session/plate;
% each subfolder holds Z-stack TIFFs named "Surf_M<position>_<timepoint>_..."
% (optionally with a "488"/"561" channel tag if the acquisition was
% dual-channel). For every stage position (M) within a subfolder, this:
%   - Thresholds every timepoint's full Z-stack at two different
%     percentile-based cutoffs: a strict ("fine", 99th percentile) one
%     that also gets watershed-split to separate touching/merged objects,
%     and a looser ("coarse", 95th percentile) one without splitting --
%     giving both a detailed and a lenient view of "what's segmented".
%   - Extracts a "surface" projection (a small max-intensity window
%     around a fixed Z-slice, surfIdx) per timepoint at both thresholds,
%     and tracks that surface signal across time.
%   - Writes binary mask stacks (fine/watershed/coarse, full-stack and
%     surface-only), per-position 2D (per-slice) and 3D (per-object)
%     measurement tables as CSVs, and a per-position summary .mat file.
%
% Positions (M) are processed in parallel via a local parpool sized to
% the SLURM allocation (SLURM_CPUS_PER_TASK), since each position is an
% independent timelapse. A position is skipped entirely if its summary
% .mat already exists, so a killed/resubmitted cluster job resumes
% instead of redoing already-finished positions.

% Size the local parallel pool to exactly what SLURM gave this job.
% Reuses an existing pool if it already has the right number of workers;
% otherwise tears down any mismatched pool and starts a fresh one sized
% correctly. JobStorageLocation is pointed at /tmp/ so pool metadata
% doesn't collide with other concurrent jobs on shared cluster storage.
nCores = str2double(getenv('SLURM_CPUS_PER_TASK'));
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

% get folders
contents = dir(path);
is_subdir = [contents.isdir] & ~ismember({contents.name}, {'.', '..'});
subfolders = contents(is_subdir);
subfolderNames = {subfolders.name};

for d = 1:length(subfolderNames)

    folderPath = fullfile(path, subfolderNames{d});
    fprintf('Processing folder %s\n', subfolderNames{d});

    % ===== folders =====
    binDir          = fullfile(folderPath, 'Binary');
    binwDir         = fullfile(folderPath, 'Binary_Watershed');
    binCoarseDir    = fullfile(folderPath, 'Binary_Coarse');
    surfBinDir      = fullfile(folderPath, 'Surface_Binary');
    surfBinwDir     = fullfile(folderPath, 'Surface_Binary_Watershed');
    surfBinCoarseDir= fullfile(folderPath, 'Surface_Binary_Coarse');
    sumDir          = fullfile(folderPath, 'Summary');

    if ~exist(binDir, 'dir'),           mkdir(binDir);           end
    if ~exist(binwDir, 'dir'),          mkdir(binwDir);          end
    if ~exist(binCoarseDir, 'dir'),     mkdir(binCoarseDir);     end
    if ~exist(surfBinDir, 'dir'),       mkdir(surfBinDir);       end
    if ~exist(surfBinwDir, 'dir'),      mkdir(surfBinwDir);      end
    if ~exist(surfBinCoarseDir, 'dir'), mkdir(surfBinCoarseDir); end
    if ~exist(sumDir, 'dir'),           mkdir(sumDir);           end

    % Exclude any binary masks already written by a previous run of this
    % script from the input file list, so re-running doesn't try to
    % segment its own output.
    allfiles = dir(fullfile(folderPath, "*.tif"));
    allfiles(contains(string({allfiles.name}), "Binary")) = [];
    allnames = string({allfiles.name});

    % detect mode
    % Dual-channel acquisitions tag filenames with 488/561; single-channel
    % data has neither. In dual-channel mode, 488 is treated as the
    % "base" channel that drives the timepoint loop below, with its 561
    % partner processed alongside it (see isPaired branches further down).
    isPaired = any(contains(allnames, "488")) || any(contains(allnames, "561"));
    if isPaired
        baseMask = contains(allnames, "488");
        basefiles = allfiles(baseMask);
        basenames = string({basefiles.name});
    else
        basefiles = allfiles;
        basenames = string({basefiles.name});
    end

    % extract positions and time points
    % Filenames encode stage position (M) and timepoint (second number)
    % as "Surf_M<M>_<timepoint>_...". Every (M, timepoint) pair present
    % is recorded so each position's own file list can be reconstructed
    % and time-sorted below.
    tok = regexp(basenames, 'Surf_M(\d+)_(\d+)_', 'tokens', 'once');
    col1 = str2double(cellfun(@(x) x{1}, tok, 'UniformOutput', false));
    col2 = str2double(cellfun(@(x) x{2}, tok, 'UniformOutput', false));
    MY = table(col1', col2', 'VariableNames', {'M', 'Y'});
    M = unique(MY.M);

    % Each stage position is an independent timelapse, so positions are
    % farmed out across the parallel pool rather than run one at a time.
    parfor p = 1:length(M)

        thisM = M(p);

        % skip only if entire M done
        % The per-position Summary.mat is the very last file this loop
        % writes, so its presence means a prior run already finished this
        % position -- lets a resubmitted/resumed cluster job skip
        % everything already completed instead of redoing it.
        summaryCheck = fullfile(sumDir, "M" + thisM + "_Surface_Summary.mat");
        if isfile(summaryCheck)
            fprintf('Skipping completed outputs for M%d\n', thisM);
            continue
        end

        idx = find(MY.M == thisM);
        [~, ord] = sort(MY.Y(idx));
        files = basefiles(idx(ord));
        surfIdx = 26;  % fixed Z-slice treated as "the surface" for the surface-projection outputs

        % initialize vars
        T2s561 = table();
        T2sw561 = table();
        T2sCoarse561 = table();
        surfbinim561 = [];
        surfbinimw561 = [];
        surfbinimCoarse561 = [];
        T2Coarse561 = table();
        T3Coarse561 = table();

        % Per-timepoint results are collected into pre-sized cell arrays
        % (rather than growing tables/arrays in the loop) so each
        % iteration can run independently -- required for correctness
        % inside a parfor-nested loop, and also avoids slow repeated
        % concatenation. Everything is joined together (cat/vertcat) once
        % all timepoints for this position are done.
        surfbinimCell = cell(length(files),1);
        surfbinimwCell = cell(length(files),1);
        surfbinimCoarseCell = cell(length(files),1);
        T2cell = cell(length(files),1);
        T3cell = cell(length(files),1);
        T2wcell = cell(length(files),1);
        T3wcell = cell(length(files),1);
        T2CoarseCell = cell(length(files),1);
        T3CoarseCell = cell(length(files),1);
        surfbinimCell561 = {};
        surfbinimwCell561 = {};
        surfbinimCoarseCell561 = {};
        T2cell561 = {};
        T3cell561 = {};
        T2wcell561 = {};
        T3wcell561 = {};
        T2CoarseCell561 = {};
        T3CoarseCell561 = {};

        if isPaired
            surfbinimCell561 = cell(length(files),1);
            surfbinimwCell561 = cell(length(files),1);
            surfbinimCoarseCell561 = cell(length(files),1);
            T2cell561 = cell(length(files),1);
            T3cell561 = cell(length(files),1);
            T2wcell561 = cell(length(files),1);
            T3wcell561 = cell(length(files),1);
            T2CoarseCell561 = cell(length(files),1);
            T3CoarseCell561 = cell(length(files),1);
        end

        for f = 1:length(files)

            % -------- primary file --------
            imname = files(f).name;
            basename = extractBefore(imname, ".tif");
            tok2 = regexp(basename, 'Surf_M(\d+)_(\d+)_', 'tokens', 'once');
            Y = str2double(tok2{2});
            impath = fullfile(folderPath, imname);

            fprintf('Generating outputs for %s\n', imname);

            % Segments and measures this one timepoint's Z-stack at both
            % thresholds (plus watershed-split and surface-projection
            % variants) -- see process_one_stack below for the actual
            % segmentation logic.
            [T2new, T3new, T2wnew, T3wnew, T2CoarseNew, T3CoarseNew, ...
             binsurf, binsurfw, binsurfCoarse, binim, binimw, binimCoarse, dims] = ...
                process_one_stack(impath, Y, surfIdx);

            T2cell{f} = T2new;
            T3cell{f} = T3new;
            T2wcell{f} = T2wnew;
            T3wcell{f} = T3wnew;
            surfbinimCell{f} = binsurf;
            surfbinimwCell{f} = binsurfw;
            surfbinimCoarseCell{f} = binsurfCoarse;
            T2CoarseCell{f} = T2CoarseNew;
            T3CoarseCell{f} = T3CoarseNew;

            % ===== ACTIVE: fine-threshold stack outputs =====
            % Writes this timepoint's fine-threshold full 3D binary stack
            % and its watershed-split counterpart as TIFF + .mat.
            outBin  = fullfile(binDir, [basename, '_Binary.tif']);
            outBinw = fullfile(binwDir, [basename, '_Binary_Watershed.tif']);
            outMat  = fullfile(binDir, [basename, '_Binary.mat']);
            outMatw = fullfile(binwDir, [basename, '_Binary_Watershed.mat']);

            for j = 1:dims(3)
                binsave  = uint8(binim(:,:,j)) * 255;
                binwsave = uint8(binimw(:,:,j)) * 255;

                if j == 1
                    imwrite(binsave,  outBin,  'tif', 'Compression', 'none');
                    imwrite(binwsave, outBinw, 'tif', 'Compression', 'none');
                else
                    imwrite(binsave,  outBin,  'tif', 'WriteMode', 'append', 'Compression', 'none');
                    imwrite(binwsave, outBinw, 'tif', 'WriteMode', 'append', 'Compression', 'none');
                end
            end

            addImageJTimeMetadata(outBin, dims(3));
            addImageJTimeMetadata(outBinw, dims(3));

            Str = struct();
            Str.binim = binim;
            save(outMat, '-fromstruct', Str);

            Strw = struct();
            Strw.binimw = binimw;
            save(outMatw, '-fromstruct', Strw);

            % ===== ACTIVE: coarse stack only =====
            outBinCoarse = fullfile(binCoarseDir, [basename, '_Binary_Coarse.tif']);
            outMatCoarse = fullfile(binCoarseDir, [basename, '_Binary_Coarse.mat']);

            for j = 1:dims(3)
                binCoarseSave = uint8(binimCoarse(:,:,j)) * 255;

                if j == 1
                    imwrite(binCoarseSave, outBinCoarse, 'tif', 'Compression', 'none');
                else
                    imwrite(binCoarseSave, outBinCoarse, 'tif', 'WriteMode', 'append', 'Compression', 'none');
                end
            end

            addImageJTimeMetadata(outBinCoarse, dims(3));

            StrCoarse = struct();
            StrCoarse.binimCoarse = binimCoarse;
            save(outMatCoarse, '-fromstruct', StrCoarse);

            fprintf('Saved:\n%s\n%s\n', outBinCoarse, outMatCoarse);

            % -------- paired 561 --------
            % Same segmentation/output steps as the 488 (or single-
            % channel) file above, run again on this timepoint's matching
            % 561 file, reusing the timepoint (Y) parsed from the 488
            % filename since the pair shares one acquisition timestamp.
            if isPaired
                imname561 = replace(imname, "488", "561");
                basename561 = extractBefore(imname561, ".tif");
                impath561 = fullfile(folderPath, imname561);

                fprintf('Generating outputs for %s\n', imname561);

                [T2new561, T3new561, T2wnew561, T3wnew561, T2CoarseNew561, T3CoarseNew561, ...
                 binsurf561, binsurfw561, binsurfCoarse561, binim561, binimw561, binimCoarse561, dims561] = ...
                    process_one_stack(impath561, Y, surfIdx);

                % exact old results preserved in code, not used unless uncommented
                T2cell561{f} = T2new561;
                T3cell561{f} = T3new561;
                T2wcell561{f} = T2wnew561;
                T3wcell561{f} = T3wnew561;

                surfbinimCell561{f} = binsurf561;
                surfbinimwCell561{f} = binsurfw561;
                surfbinimCoarseCell561{f} = binsurfCoarse561;
                T2CoarseCell561{f} = T2CoarseNew561;
                T3CoarseCell561{f} = T3CoarseNew561;

                % ===== ACTIVE: fine-threshold stack outputs (561 partner) =====
                outBin561  = fullfile(binDir, [basename561, '_Binary.tif']);
                outBinw561 = fullfile(binwDir, [basename561, '_Binary_Watershed.tif']);
                outMat561  = fullfile(binDir, [basename561, '_Binary.mat']);
                outMatw561 = fullfile(binwDir, [basename561, '_Binary_Watershed.mat']);

                for j = 1:dims561(3)
                    binsave561  = uint8(binim561(:,:,j)) * 255;
                    binwsave561 = uint8(binimw561(:,:,j)) * 255;

                    if j == 1
                        imwrite(binsave561,  outBin561,  'tif', 'Compression', 'none');
                        imwrite(binwsave561, outBinw561, 'tif', 'Compression', 'none');
                    else
                        imwrite(binsave561,  outBin561,  'tif', 'WriteMode', 'append', 'Compression', 'none');
                        imwrite(binwsave561, outBinw561, 'tif', 'WriteMode', 'append', 'Compression', 'none');
                    end
                end

                addImageJTimeMetadata(outBin561, dims561(3));
                addImageJTimeMetadata(outBinw561, dims561(3));

                Str561 = struct();
                Str561.binim561 = binim561;
                save(outMat561, '-fromstruct', Str561);

                Strw561 = struct();
                Strw561.binimw561 = binimw561;
                save(outMatw561, '-fromstruct', Strw561);

                % ===== ACTIVE: coarse stack only =====
                outBinCoarse561 = fullfile(binCoarseDir, [basename561, '_Binary_Coarse.tif']);
                outMatCoarse561 = fullfile(binCoarseDir, [basename561, '_Binary_Coarse.mat']);

                for j = 1:dims561(3)
                    binCoarseSave561 = uint8(binimCoarse561(:,:,j)) * 255;

                    if j == 1
                        imwrite(binCoarseSave561, outBinCoarse561, 'tif', 'Compression', 'none');
                    else
                        imwrite(binCoarseSave561, outBinCoarse561, 'tif', 'WriteMode', 'append', 'Compression', 'none');
                    end
                end

                addImageJTimeMetadata(outBinCoarse561, dims561(3));

                StrCoarse561 = struct();
                StrCoarse561.binimCoarse561 = binimCoarse561;
                save(outMatCoarse561, '-fromstruct', StrCoarse561);

                fprintf('Saved:\n%s\n%s\n', outBinCoarse561, outMatCoarse561);
            end
        end

        % gather active results
        % Per-timepoint surface binaries are stacked along the 3rd
        % dimension into one time-series volume per channel/threshold;
        % per-timepoint coarse 3D tables are stacked into one long table.
        surfbinim = cat(3, surfbinimCell{:});
        surfbinimw = cat(3, surfbinimwCell{:});
        surfbinimCoarse = cat(3, surfbinimCoarseCell{:});
        T2Coarse = vertcat(T2CoarseCell{:});
        T3Coarse = vertcat(T3CoarseCell{:});

        if isPaired
            surfbinim561 = cat(3, surfbinimCell561{:});
            surfbinimw561 = cat(3, surfbinimwCell561{:});
            surfbinimCoarse561 = cat(3, surfbinimCoarseCell561{:});
            T2Coarse561 = vertcat(T2CoarseCell561{:});
            T3Coarse561 = vertcat(T3CoarseCell561{:});
        end

        % ===== ACTIVE: fine-threshold stack summaries =====
        % Concatenates every timepoint's fine-threshold 2D/3D and
        % watershed tables for this position and writes them out as
        % per-position CSVs.
        T2 = vertcat(T2cell{:});
        T3 = vertcat(T3cell{:});
        T2w = vertcat(T2wcell{:});
        T3w = vertcat(T3wcell{:});

        if isPaired
            T2_561 = vertcat(T2cell561{:});
            T3_561 = vertcat(T3cell561{:});
            T2w_561 = vertcat(T2wcell561{:});
            T3w_561 = vertcat(T3wcell561{:});

            writetable(T2,      fullfile(sumDir, "M" + thisM + "_2D_488.csv"));
            writetable(T3,      fullfile(sumDir, "M" + thisM + "_3D_488.csv"));
            writetable(T2w,     fullfile(sumDir, "M" + thisM + "_2D_Watershed_488.csv"));
            writetable(T3w,     fullfile(sumDir, "M" + thisM + "_3D_Watershed_488.csv"));
            writetable(T2_561,  fullfile(sumDir, "M" + thisM + "_2D_561.csv"));
            writetable(T3_561,  fullfile(sumDir, "M" + thisM + "_3D_561.csv"));
            writetable(T2w_561, fullfile(sumDir, "M" + thisM + "_2D_Watershed_561.csv"));
            writetable(T3w_561, fullfile(sumDir, "M" + thisM + "_3D_Watershed_561.csv"));

            datacell_old = {T2, T2w, T3, T3w, T2_561, T2w_561, T3_561, T3w_561};
        else
            T2 = vertcat(T2cell{:});
            T3 = vertcat(T3cell{:});
            T2w = vertcat(T2wcell{:});
            T3w = vertcat(T3wcell{:});

            writetable(T2,  fullfile(sumDir, "M" + thisM + "_2D.csv"));
            writetable(T3,  fullfile(sumDir, "M" + thisM + "_3D.csv"));
            writetable(T2w, fullfile(sumDir, "M" + thisM + "_2D_Watershed.csv"));
            writetable(T3w, fullfile(sumDir, "M" + thisM + "_3D_Watershed.csv"));

            datacell_old = {T2, T2w, T3, T3w};
        end

        % ===== coarse stack summaries + combined datacell =====
        if isPaired
            writetable(T2Coarse,    fullfile(sumDir, "M" + thisM + "_2D_Coarse_488.csv"));
            writetable(T3Coarse,    fullfile(sumDir, "M" + thisM + "_3D_Coarse_488.csv"));
            writetable(T2Coarse561, fullfile(sumDir, "M" + thisM + "_2D_Coarse_561.csv"));
            writetable(T3Coarse561, fullfile(sumDir, "M" + thisM + "_3D_Coarse_561.csv"));

            datacell = {T2, T2w, T3, T3w, ...
                        T2_561, T2w_561, T3_561, T3w_561, ...
                        T2Coarse, T3Coarse, T2Coarse561, T3Coarse561};
        else
            writetable(T2Coarse, fullfile(sumDir, "M" + thisM + "_2D_Coarse.csv"));
            writetable(T3Coarse, fullfile(sumDir, "M" + thisM + "_3D_Coarse.csv"));

            datacell = {T2, T2w, T3, T3w, T2Coarse, T3Coarse};
        end

        % ===== ACTIVE surface analyses =====
        % Runs the surface (near-coverslip max-projection) time-series
        % analysis across all of this position's timepoints at once --
        % see process_surface_stack / process_surface_stack_coarse below.
        [T2s, T2sw, ~] = process_surface_stack(surfbinim, length(files));
        T2sCoarse = process_surface_stack_coarse(surfbinimCoarse, length(files));

        if isPaired
            surfname488 = "Surf_M" + thisM + "_488";
            surfBin       = fullfile(surfBinDir,       surfname488 + "_Binary.tif");
            surfMat       = fullfile(surfBinDir,       surfname488 + "_Binary.mat");
            surfBinw      = fullfile(surfBinwDir,      surfname488 + "_Binary_Watershed.tif");
            surfMatw      = fullfile(surfBinwDir,      surfname488 + "_Binary_Watershed.mat");
            surfBinCoarse = fullfile(surfBinCoarseDir, surfname488 + "_Binary_Coarse.tif");
            surfMatCoarse = fullfile(surfBinCoarseDir, surfname488 + "_Binary_Coarse.mat");
        else
            surfname = "Surf_M" + thisM;
            surfBin       = fullfile(surfBinDir,       surfname + "_Binary.tif");
            surfMat       = fullfile(surfBinDir,       surfname + "_Binary.mat");
            surfBinw      = fullfile(surfBinwDir,      surfname + "_Binary_Watershed.tif");
            surfMatw      = fullfile(surfBinwDir,      surfname + "_Binary_Watershed.mat");
            surfBinCoarse = fullfile(surfBinCoarseDir, surfname + "_Binary_Coarse.tif");
            surfMatCoarse = fullfile(surfBinCoarseDir, surfname + "_Binary_Coarse.mat");
        end

        for k = 1:length(files)
            surfbinsave       = uint8(surfbinim(:,:,k)) * 255;
            surfbinwsave      = uint8(surfbinimw(:,:,k)) * 255;
            surfbinCoarseSave = uint8(surfbinimCoarse(:,:,k)) * 255;

            if k == 1
                imwrite(surfbinsave,       surfBin,       'tif', 'Compression', 'none');
                imwrite(surfbinwsave,      surfBinw,      'tif', 'Compression', 'none');
                imwrite(surfbinCoarseSave, surfBinCoarse, 'tif', 'Compression', 'none');
            else
                imwrite(surfbinsave,       surfBin,       'tif', 'WriteMode', 'append', 'Compression', 'none');
                imwrite(surfbinwsave,      surfBinw,      'tif', 'WriteMode', 'append', 'Compression', 'none');
                imwrite(surfbinCoarseSave, surfBinCoarse, 'tif', 'WriteMode', 'append', 'Compression', 'none');
            end
        end

        addImageJTimeMetadata(surfBin, length(files));
        addImageJTimeMetadata(surfBinw, length(files));
        addImageJTimeMetadata(surfBinCoarse, length(files));

        Ssurf = struct();
        Ssurf.surfbinim = surfbinim;
        save(surfMat, '-fromstruct', Ssurf);

        Ssurfw = struct();
        Ssurfw.surfbinimw = surfbinimw;
        save(surfMatw, '-fromstruct', Ssurfw);

        SsurfCoarse = struct();
        SsurfCoarse.surfbinimCoarse = surfbinimCoarse;
        save(surfMatCoarse, '-fromstruct', SsurfCoarse);

        if isPaired
            [T2s561, T2sw561, ~] = process_surface_stack(surfbinim561, length(files));
            T2sCoarse561 = process_surface_stack_coarse(surfbinimCoarse561, length(files));

            surfname561 = "Surf_M" + thisM + "_561";
            surfBin561       = fullfile(surfBinDir,       surfname561 + "_Binary.tif");
            surfMat561       = fullfile(surfBinDir,       surfname561 + "_Binary.mat");
            surfBinw561      = fullfile(surfBinwDir,      surfname561 + "_Binary_Watershed.tif");
            surfMatw561      = fullfile(surfBinwDir,      surfname561 + "_Binary_Watershed.mat");
            surfBinCoarse561 = fullfile(surfBinCoarseDir, surfname561 + "_Binary_Coarse.tif");
            surfMatCoarse561 = fullfile(surfBinCoarseDir, surfname561 + "_Binary_Coarse.mat");

            for k = 1:length(files)
                surfbinsave561       = uint8(surfbinim561(:,:,k)) * 255;
                surfbinwsave561      = uint8(surfbinimw561(:,:,k)) * 255;
                surfbinCoarseSave561 = uint8(surfbinimCoarse561(:,:,k)) * 255;

                if k == 1
                    imwrite(surfbinsave561,       surfBin561,       'tif', 'Compression', 'none');
                    imwrite(surfbinwsave561,      surfBinw561,      'tif', 'Compression', 'none');
                    imwrite(surfbinCoarseSave561, surfBinCoarse561, 'tif', 'Compression', 'none');
                else
                    imwrite(surfbinsave561,       surfBin561,       'tif', 'WriteMode', 'append', 'Compression', 'none');
                    imwrite(surfbinwsave561,      surfBinw561,      'tif', 'WriteMode', 'append', 'Compression', 'none');
                    imwrite(surfbinCoarseSave561, surfBinCoarse561, 'tif', 'WriteMode', 'append', 'Compression', 'none');
                end
            end

            addImageJTimeMetadata(surfBin561, length(files));
            addImageJTimeMetadata(surfBinw561, length(files));
            addImageJTimeMetadata(surfBinCoarse561, length(files));

            Ssurf561 = struct();
            Ssurf561.surfbinim561 = surfbinim561;
            save(surfMat561, '-fromstruct', Ssurf561);

            Ssurfw561 = struct();
            Ssurfw561.surfbinimw561 = surfbinimw561;
            save(surfMatw561, '-fromstruct', Ssurfw561);

            SsurfCoarse561 = struct();
            SsurfCoarse561.surfbinimCoarse561 = surfbinimCoarse561;
            save(surfMatCoarse561, '-fromstruct', SsurfCoarse561);

            writetable(T2s,          fullfile(sumDir, "M" + thisM + "_Surface_2D_488.csv"));
            writetable(T2sw,         fullfile(sumDir, "M" + thisM + "_Surface_2D_Watershed_488.csv"));
            writetable(T2s561,       fullfile(sumDir, "M" + thisM + "_Surface_2D_561.csv"));
            writetable(T2sw561,      fullfile(sumDir, "M" + thisM + "_Surface_2D_Watershed_561.csv"));
            writetable(T2sCoarse,    fullfile(sumDir, "M" + thisM + "_Surface_2D_Coarse_488.csv"));
            writetable(T2sCoarse561, fullfile(sumDir, "M" + thisM + "_Surface_2D_Coarse_561.csv"));

            surfcell = {T2s, T2sw, T2s561, T2sw561, T2sCoarse, T2sCoarse561};
        else
            writetable(T2s,       fullfile(sumDir, "M" + thisM + "_Surface_2D.csv"));
            writetable(T2sw,      fullfile(sumDir, "M" + thisM + "_Surface_2D_Watershed.csv"));
            writetable(T2sCoarse, fullfile(sumDir, "M" + thisM + "_Surface_2D_Coarse.csv"));

            surfcell = {T2s, T2sw, T2sCoarse};
        end

        % Last file written for this position -- its existence is what
        % the skip-check at the top of this parfor iteration looks for.
        surfdataMat = fullfile(sumDir, "M" + thisM + "_Surface_Summary.mat");
        Ssurfdata = struct();
        Ssurfdata.surfcell = surfcell;
        save(surfdataMat, '-fromstruct', Ssurfdata);

    end
end

end

%%
function [T2new, T3new, T2wnew, T3wnew, T2Coarse, T3Coarse, binsurf, binsurfw, binsurfCoarse, binim, binimw, binimCoarse, dims] = process_one_stack(imname, Y, surfIdx)
% Segments and measures one timepoint's full Z-stack at two independent,
% image-adaptive intensity thresholds (percentiles of that stack's own
% intensity distribution, rather than a fixed absolute cutoff):
%   - "fine" (99th percentile): strict threshold, further split by 3D
%     watershed to separate touching/merged objects, giving both an
%     un-split and a split measurement set.
%   - "coarse" (95th percentile): looser threshold, measured directly
%     with no splitting -- a lighter-weight, more permissive view of
%     occupied area/volume.
% Also extracts a small max-intensity "surface" projection around
% surfIdx (the fixed Z-slice treated as the coverslip/surface) at both
% thresholds, watershed-splitting the fine one.
%
% Returns per-slice 2D tables (Txnew, one row per detected cross-section
% per Z-slice) and per-object 3D tables (Txnew, one row per connected 3D
% object with real volume/shape) for both the fine and coarse thresholds,
% plus the raw binary volumes and the three single-slice surface masks.

im = tiffreadVolume(imname);
dims = size(im);

thr = prctile(im(:), 99);
thrCoarse = prctile(im(:), 95);

% ===== ACTIVE: fine-threshold (99th percentile) stack analysis =====
binim = false(dims(1), dims(2), dims(3));
cent2 = zeros(0,2);
ars = [];
mal  = [];
mil  = [];
nObj2 = 0;
z = [];
t = [];

for j = 1:dims(3)
    slice = im(:,:,j);

    bin = slice > thr;
    bin = bwareaopen(bin, 10);
    bin = imclose(bin, strel("disk",2));
    binim(:,:,j) = bin;

    CC2 = bwconncomp(bin);
    stats2 = regionprops(CC2, 'Area', 'Centroid', 'MajorAxisLength', 'MinorAxisLength');
    stats2 = struct2table(stats2);

    cent2 = [cent2; stats2.Centroid];
    ars = [ars; stats2.Area];
    mal  = [mal; stats2.MajorAxisLength];
    mil  = [mil; stats2.MinorAxisLength];
    z = [z; j * ones(height(stats2), 1)];
    t = [t; Y * ones(height(stats2), 1)];
    nObj2 = nObj2 + height(stats2);
end

T2new = table;
T2new.T = t;
T2new.ObjectID = (1:nObj2)';
T2new.Area = ars;
T2new.X = cent2(:,1);
T2new.Y = cent2(:,2);
T2new.Z = z;
T2new.MajorAxisLength = mal;
T2new.MinorAxisLength = mil;

CC3 = bwconncomp(binim);
stats3 = regionprops3(CC3, 'Volume', 'Centroid', 'PrincipalAxisLength');
if isempty(stats3)
    cent3 = zeros(0,3);
    vols = zeros(0,1);
    pal  = zeros(0,3);
else
    cent3 = stats3.Centroid;
    vols = stats3.Volume;
    pal  = stats3.PrincipalAxisLength;
end
nObj3 = height(stats3);
time = Y * ones(nObj3, 1);

T3new = table;
T3new.T = time;
T3new.ObjectID = (1:nObj3)';
T3new.Volume = vols;
T3new.X = cent3(:,1);
T3new.Y = cent3(:,2);
T3new.Z = round(cent3(:,3));
T3new.PrincipalAxisLength1 = pal(:,1);
T3new.PrincipalAxisLength2 = pal(:,2);
T3new.PrincipalAxisLength3 = pal(:,3);
T3new = sortrows(T3new, 'Z');

% Watershed-split every 3D object individually: crop a small padded
% sub-volume around just that object (fast -- avoids watershedding the
% whole stack at once), compute its distance transform, suppress shallow
% local maxima below a height of 0.75 (imhmax) so minor bumps in the
% distance map don't cause over-segmentation, then split via 3D watershed
% and paste the result back into the full-size binary.
binimw = false(size(binim));
for i = 1:CC3.NumObjects
    [r,c,z3] = ind2sub(size(binim), CC3.PixelIdxList{i});

    y1 = max(min(r)-1, 1);
    y2 = min(max(r)+1, size(binim,1));
    x1 = max(min(c)-1, 1);
    x2 = min(max(c)+1, size(binim,2));
    z1 = max(min(z3)-1, 1);
    z2 = min(max(z3)+1, size(binim,3));

    sub = false(y2-y1+1, x2-x1+1, z2-z1+1);
    subind = sub2ind(size(sub), r-y1+1, c-x1+1, z3-z1+1);
    sub(subind) = true;

    D = bwdist(~sub);
    D = imhmax(D, 0.75);
    L = watershed(-D);

    subw = sub;
    subw(L == 0) = 0;

    binimw(y1:y2, x1:x2, z1:z2) = binimw(y1:y2, x1:x2, z1:z2) | subw;
end

% Same per-slice 2D + whole-volume 3D measurements as above, but on the
% post-watershed (declumped) binary.
cent2 = zeros(0,2);
ars = [];
mal  = [];
mil  = [];
nObj2 = 0;
z = [];
t = [];

for j = 1:dims(3)
    binw = binimw(:,:,j);
    CC2 = bwconncomp(binw);
    stats2 = regionprops(CC2, 'Area', 'Centroid', 'MajorAxisLength', 'MinorAxisLength');
    stats2 = struct2table(stats2);

    cent2 = [cent2; stats2.Centroid];
    ars = [ars; stats2.Area];
    mal  = [mal; stats2.MajorAxisLength];
    mil  = [mil; stats2.MinorAxisLength];
    z = [z; j * ones(height(stats2), 1)];
    t = [t; Y * ones(height(stats2), 1)];
    nObj2 = nObj2 + height(stats2);
end

T2wnew = table;
T2wnew.T = t;
T2wnew.ObjectID = (1:nObj2)';
T2wnew.Area = ars;
T2wnew.X = cent2(:,1);
T2wnew.Y = cent2(:,2);
T2wnew.Z = z;
T2wnew.MajorAxisLength = mal;
T2wnew.MinorAxisLength = mil;

CC3 = bwconncomp(binimw);
stats3 = regionprops3(CC3, 'Volume', 'Centroid', 'PrincipalAxisLength');
if isempty(stats3)
    cent3 = zeros(0,3);
    vols = zeros(0,1);
    pal  = zeros(0,3);
else
    cent3 = stats3.Centroid;
    vols = stats3.Volume;
    pal  = stats3.PrincipalAxisLength;
end
nObj3 = height(stats3);
time = Y * ones(nObj3, 1);

T3wnew = table;
T3wnew.T = time;
T3wnew.ObjectID = (1:nObj3)';
T3wnew.Volume = vols;
T3wnew.X = cent3(:,1);
T3wnew.Y = cent3(:,2);
T3wnew.Z = round(cent3(:,3));
T3wnew.PrincipalAxisLength1 = pal(:,1);
T3wnew.PrincipalAxisLength2 = pal(:,2);
T3wnew.PrincipalAxisLength3 = pal(:,3);
T3wnew = sortrows(T3wnew, 'Z');

% comment out for full outputs
% Toggle: uncommenting this block blanks the fine-threshold outputs
% (2D/3D tables and both binary volumes) instead of returning them --
% useful to skip the more expensive per-slice/per-object fine analysis
% on large batch runs where only the coarse summary below is needed.
% T2new = table();
% T3new = table();
% T2wnew = table();
% T3wnew = table();
% binim = false(0,0,0);
% binimw = false(0,0,0);

%% coarse 95% stack binary
% Same per-slice 2D + whole-volume 3D measurement idea as the fine
% analysis above, but at the looser threshold and without any watershed
% splitting or shape (axis-length) columns -- a lighter, more permissive
% segmentation pass.
binimCoarse = false(dims(1), dims(2), dims(3));
cent2c = zeros(0,2);
arsc = [];
zc = [];
tc = [];
nObj2c = 0;

for j = 1:dims(3)
    slice = im(:,:,j);

    binc = slice > thrCoarse;
    binc = bwareaopen(binc, 10);
    binc = imclose(binc, strel("disk",2));
    binimCoarse(:,:,j) = binc;

    CC2c = bwconncomp(binc);
    stats2c = regionprops(CC2c, 'Area', 'Centroid');
    stats2c = struct2table(stats2c);

    cent2c = [cent2c; stats2c.Centroid];
    arsc = [arsc; stats2c.Area];
    zc = [zc; j * ones(height(stats2c), 1)];
    tc = [tc; Y * ones(height(stats2c), 1)];
    nObj2c = nObj2c + height(stats2c);
end

T2Coarse = table;
T2Coarse.T = tc;
T2Coarse.ObjectID = (1:nObj2c)';
T2Coarse.Area = arsc;
T2Coarse.X = cent2c(:,1);
T2Coarse.Y = cent2c(:,2);
T2Coarse.Z = zc;

CC3c = bwconncomp(binimCoarse);
stats3c = regionprops3(CC3c, 'Volume', 'Centroid');
if isempty(stats3c)
    cent3c = zeros(0,3);
    volc = zeros(0,1);
else
    cent3c = stats3c.Centroid;
    volc = stats3c.Volume;
end
nObj3c = height(stats3c);
timec = Y * ones(nObj3c, 1);

T3Coarse = table;
T3Coarse.T = timec;
T3Coarse.ObjectID = (1:nObj3c)';
T3Coarse.Volume = volc;
T3Coarse.X = cent3c(:,1);
T3Coarse.Y = cent3c(:,2);
T3Coarse.Z = round(cent3c(:,3));
T3Coarse = sortrows(T3Coarse, 'Z');

%% surface binaries
% Builds a single "at the surface" projection for this timepoint (max
% intensity over a small +/-3 slice window around surfIdx), then
% thresholds it at both cutoffs and watershed-splits the fine one --
% same distance-transform + h-maxima-suppression + watershed recipe as
% above, just applied directly to one 2D slice instead of per-3D-object.
imsurf = max(im(:,:,surfIdx-3:surfIdx+3), [], 3);

binsurf = imsurf > thr;
binsurf = bwareaopen(binsurf, 10);
binsurf = imclose(binsurf, strel("disk",2));

Ds = bwdist(~binsurf);
Ds = imhmax(Ds, 0.75);
Ls = watershed(-Ds);
binsurfw = binsurf;
binsurfw(Ls == 0) = 0;

binsurfCoarse = imsurf > thrCoarse;
binsurfCoarse = bwareaopen(binsurfCoarse, 10);
binsurfCoarse = imclose(binsurfCoarse, strel("disk",2));

end

%%
function [T2s, T2sw, surfbinimw] = process_surface_stack(surfbinim, nFiles)
% Given the fine-threshold surface projection stacked across every
% timepoint (one 2D slice per timepoint), measures 2D objects at each
% timepoint independently on both the as-is mask (T2s) and a freshly
% watershed-split version of it (T2sw) -- the split is redone here per
% timepoint (same distance-transform/h-maxima/watershed recipe used
% elsewhere) since this operates on the single-slice surface projection
% rather than the full 3D per-timepoint stack. Also returns the full
% watershed-split surface volume so the caller can save it as an output.

surfbinimw = false(size(surfbinim));

cent2s = zeros(0,2);
arss = [];
mals  = [];
mils  = [];
nObj2s = 0;
ts = [];
cent2sw = zeros(0,2);
arssw = [];
malsw  = [];
milsw  = [];
nObj2sw = 0;
tsw = [];

for s = 1:nFiles
    bin2d = surfbinim(:,:,s);
    Ds = bwdist(~bin2d);
    Ds = imhmax(Ds, 0.75);
    Ls = watershed(-Ds);

    bin2dw = bin2d;
    bin2dw(Ls == 0) = 0;
    surfbinimw(:,:,s) = bin2dw;

    CC2s = bwconncomp(bin2d);
    stats2s = regionprops(CC2s, 'Area', 'Centroid', 'MajorAxisLength', 'MinorAxisLength');
    stats2s = struct2table(stats2s);
    cent2s = [cent2s; stats2s.Centroid];
    arss = [arss; stats2s.Area];
    mals  = [mals; stats2s.MajorAxisLength];
    mils  = [mils; stats2s.MinorAxisLength];
    ts = [ts; s * ones(height(stats2s), 1)];
    nObj2s = nObj2s + height(stats2s);

    CC2sw = bwconncomp(bin2dw);
    stats2sw = regionprops(CC2sw, 'Area', 'Centroid', 'MajorAxisLength', 'MinorAxisLength');
    stats2sw = struct2table(stats2sw);
    cent2sw = [cent2sw; stats2sw.Centroid];
    arssw = [arssw; stats2sw.Area];
    malsw  = [malsw; stats2sw.MajorAxisLength];
    milsw  = [milsw; stats2sw.MinorAxisLength];
    tsw = [tsw; s * ones(height(stats2sw), 1)];
    nObj2sw = nObj2sw + height(stats2sw);
end

T2s = table;
T2s.T = ts;
T2s.ObjectID = (1:nObj2s)';
T2s.Area = arss;
T2s.X = cent2s(:,1);
T2s.Y = cent2s(:,2);
T2s.MajorAxisLength = mals;
T2s.MinorAxisLength = mils;

T2sw = table;
T2sw.T = tsw;
T2sw.ObjectID = (1:nObj2sw)';
T2sw.Area = arssw;
T2sw.X = cent2sw(:,1);
T2sw.Y = cent2sw(:,2);
T2sw.MajorAxisLength = malsw;
T2sw.MinorAxisLength = milsw;

end

%%
function T2sCoarse = process_surface_stack_coarse(surfbinimCoarse, nFiles)
% Same idea as process_surface_stack above but for the coarse-threshold
% surface projection: per-timepoint 2D Area/Centroid only, no watershed
% splitting.

cent2s = zeros(0,2);
arss = [];
nObj2s = 0;
ts = [];

for s = 1:nFiles
    bin2d = surfbinimCoarse(:,:,s);
    CC2s = bwconncomp(bin2d);
    stats2s = regionprops(CC2s, 'Area', 'Centroid');
    stats2s = struct2table(stats2s);

    cent2s = [cent2s; stats2s.Centroid];
    arss = [arss; stats2s.Area];
    ts = [ts; s * ones(height(stats2s), 1)];
    nObj2s = nObj2s + height(stats2s);
end

T2sCoarse = table;
T2sCoarse.T = ts;
T2sCoarse.ObjectID = (1:nObj2s)';
T2sCoarse.Area = arss;
T2sCoarse.X = cent2s(:,1);
T2sCoarse.Y = cent2s(:,2);

end

%%
function addImageJTimeMetadata(fname, nT)
% Stamps the given TIFF with an ImageJ hyperstack ImageDescription tag
% declaring it as a single-channel, single-slice, nT-frame time series,
% so ImageJ opens the appended stack as a proper T-series instead of an
% unlabeled stack.
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
