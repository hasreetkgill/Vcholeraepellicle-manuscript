%% SummarizeAndPlotCultureData.m
% Author: Hasreet Gill, Bassler Lab, HHMI/Princeton University
% Revised with Claude (Sonnet 5, Anthropic)

% Final analysis/plotting stage of the culture pipeline. Consumes the
% per-position Summary CSVs written by SegmentCultureTimelapse.m and
% CompileCultureTracking_SurfaceMiddleSummaries.m, plus two lookup
% workbooks (metadata.xlsx: per-experiment-date TimeStep and each stage
% position's genotype-code pairing; genotypes.xlsx: numeric code -> name),
% and produces:
%   1. Per-genotype, per-metric, per-data-type replicate-averaged time
%      courses (mean +/- SEM across dates/positions sharing a genotype),
%      for both the Surface_FromStacks and original Tracking/Binary/etc.
%      tables -- plus a "Surface normalized to Middle" version of each
%      Surface_FromStacks time course, dividing by the matching
%      Middle_FromStacks time course to correct for overall biofilm
%      growth/thickness trends (SURFACE PLOTS section).
%   2. Per-replicate (not averaged) "stack profile" plots of each metric
%      vs Z-slice, one line per timepoint, from the un-sliced full-stack
%      tables -- since these visualize 3D shape, not cross-replicate
%      reproducibility (STACK PROFILE PLOTS section).
%   3. For Coculture data specifically, the mean nearest-neighbor
%      distance from each red (561) object to the nearest qualifying
%      green (488) object per timepoint, averaged across replicates of
%      the same genotype pair (RED-TO-GREEN DISTANCE section).
%
% Every section re-derives the same culture-type/date-folder listing
% independently (rather than sharing state), so each can be run/re-run on
% its own.

clear
clc

parentPath = uigetdir();

metaFile = fullfile(parentPath,"metadata.xlsx");
genoFile = fullfile(parentPath,"genotypes.xlsx");

meta = readtable(metaFile);
geno = readtable(genoFile);

%% ---- User controls ----

% Approximate single-cell footprint/volume: used both to drop sub-single-
% cell noise objects and to convert a summed area/volume into an
% "estimated number of cells".
CELL_AREA = 120;
CELL_VOLUME = 1500;

% Size cutoff used to split objects into "bigger than" vs "smaller than"
% this threshold counts (a simple static size split, not a growth-based
% classification like the microcolony logic in the single-cell pipeline).
THRESH_AREA = 400;
THRESH_VOLUME = 4500;

% Track quality filters (max displacement / min duration) for the
% ORIGINAL compiled tracking tables (Tracking / Tracking_Watershed, built
% straight from TrackMate output before any surface/middle re-slicing).
ORIG_TRACK_MAX_DISP = 100;
ORIG_TRACK_MIN_DUR = 3;

% Same idea, but for the Surface_FromStacks_*_Tracking* /
% Middle_FromStacks_*_Tracking* tables (re-sliced from the full stack).
STACK_TRACK_MAX_DISP = 100;
STACK_TRACK_MIN_DUR = 3;

% Max acceptable XY distance for the earlier 2D-track-to-3D-object match
% (TrackingDistanceXY, from CompileCultureTracking...'s step 3) to be
% trusted; NaN (no match attempted) is kept rather than dropped.
TRACK_MATCH_MAX_DIST = 50;

% Separate, stricter size cutoffs used only in the red-to-green distance
% analysis, to decide which green objects count as real neighbor targets
% (distinct from the general CELL_AREA/CELL_VOLUME noise filter above).
NEIGHBOR_THRESH_AREA = 240;
NEIGHBOR_THRESH_VOLUME = 3000;

% The five per-frame summary metrics computed everywhere below.
metricNames = { ...
    "TotalMeasure", ...
    "EstimatedCells", ...
    "TotalObjects", ...
    "NumObjectsGTThresh", ...
    "NumObjectsLTThresh"};

%% ---- Convert metadata dates ----
% meta.Date is stored as MMDDYY; convert to the "20YYMMDD" form used by
% the date-folder names everywhere else in this pipeline.

metaDates_raw = string(meta.Date);
dates_meta = strings(size(metaDates_raw));

for i = 1:length(metaDates_raw)

    d  = metaDates_raw(i);
    mm = extractBetween(d,1,2);
    dd = extractBetween(d,3,4);
    yy = extractBetween(d,5,6);

    dates_meta(i) = "20" + yy + mm + dd;

end

% Which top-level dataset folder(s) to process -- toggle by
% commenting/uncommenting.
% cultureTypes = ["Monoculture","Coculture"];
cultureTypes = ["Protein addition"];

for c = 1:length(cultureTypes)

    cultureName = cultureTypes(c);
    disp(" ")
    disp("Processing " + cultureName)

    culturePath = fullfile(parentPath,cultureName);
    outDir = fullfile(culturePath,"SurfacePlots");

    if ~exist(outDir,'dir')
        mkdir(outDir);
    end

    genoData = struct();

    dateFolders = dir(culturePath);
    dateFolders = dateFolders([dateFolders.isdir]);
    dateFolders = dateFolders(~ismember({dateFolders.name},{'.','..','SurfacePlots'}));

    %% ===============================
    %% COLLECT DATA
    %% ===============================
    % Gathers every Surface_FromStacks_*/Middle_FromStacks_* summary CSV
    % (the base, un-sliced M*_2D*.csv/M*_3D*.csv tables are handled
    % separately in STACK PROFILE PLOTS below) across every date, and
    % accumulates each one as one more replicate under its genotype /
    % metric / data-type / channel combination for averaging.

    for d = 1:length(dateFolders)

        dateName = string(dateFolders(d).name);
        summaryPath = fullfile(culturePath,dateName,"Summary");

        metaIdx = find(dates_meta == dateName);

        % Dates with no matching metadata.xlsx row contribute nothing --
        % there's no TimeStep or genotype-pairing info to use for them.
        if isempty(metaIdx)
            continue;
        end

        timeStep = meta.TimeStep(metaIdx);

        Mcols = startsWith(meta.Properties.VariableNames,"M");
        XYcells = meta{metaIdx,Mcols};

        csvFiles = dir(fullfile(summaryPath,"M*.csv"));

        for f = 1:length(csvFiles)

            fname = csvFiles(f).name;

            if ~contains(fname,"Surface") && ~contains(fname,"Middle")
                continue;
            end

            if contains(fname,"SurfaceSliceSummary")
                continue;
            end

            tok = regexp(fname,"M(\d+)_",'tokens','once');
            M = str2double(tok{1});

            if M > length(XYcells)
                continue;
            end

            cellVal = XYcells{M};

            if isempty(cellVal)
                continue;
            end

            % Each M-column metadata cell is "redCode_greenCode" -- the
            % two genotype codes paired at that stage position.
            parts = split(string(cellVal),"_");
            redCode = str2double(parts(1));
            greenCode = str2double(parts(2));

            redIdx = find(geno.Var1 == redCode);
            greenIdx = find(geno.Var1 == greenCode);

            redName = string(geno.Var2(redIdx));
            greenName = string(geno.Var2(greenIdx));

            %% ---- dataset type ----
            % Order matters: e.g. a "3D_Tracking_Watershed" filename would
            % also satisfy contains(fname,"3D_Tracking") and
            % contains(fname,"3D"), so the most specific match must be
            % checked first in each if/elseif chain below.

            if contains(fname,"Middle_FromStacks")

                if contains(fname,"3D_Tracking_Watershed")
                    dtype = "Middle_FromStacks_3D_Tracking_Watershed";
                elseif contains(fname,"3D_Tracking")
                    dtype = "Middle_FromStacks_3D_Tracking";
                elseif contains(fname,"2D_Tracking_Watershed")
                    dtype = "Middle_FromStacks_2D_Tracking_Watershed";
                elseif contains(fname,"2D_Tracking")
                    dtype = "Middle_FromStacks_2D_Tracking";
                elseif contains(fname,"3D_Coarse")
                    dtype = "Middle_FromStacks_3D_Coarse";
                elseif contains(fname,"3D_Watershed")
                    dtype = "Middle_FromStacks_3D_Watershed";
                elseif contains(fname,"3D")
                    dtype = "Middle_FromStacks_3D";
                elseif contains(fname,"2D_Coarse")
                    dtype = "Middle_FromStacks_2D_Coarse";
                elseif contains(fname,"2D_Watershed")
                    dtype = "Middle_FromStacks_2D_Watershed";
                else
                    dtype = "Middle_FromStacks_2D";
                end

            elseif contains(fname,"Surface_FromStacks")

                if contains(fname,"3D_Tracking_Watershed")
                    dtype = "Surface_FromStacks_3D_Tracking_Watershed";
                elseif contains(fname,"3D_Tracking")
                    dtype = "Surface_FromStacks_3D_Tracking";
                elseif contains(fname,"2D_Tracking_Watershed")
                    dtype = "Surface_FromStacks_2D_Tracking_Watershed";
                elseif contains(fname,"2D_Tracking")
                    dtype = "Surface_FromStacks_2D_Tracking";
                elseif contains(fname,"3D_Coarse")
                    dtype = "Surface_FromStacks_3D_Coarse";
                elseif contains(fname,"3D_Watershed")
                    dtype = "Surface_FromStacks_3D_Watershed";
                elseif contains(fname,"3D")
                    dtype = "Surface_FromStacks_3D";
                elseif contains(fname,"2D_Coarse")
                    dtype = "Surface_FromStacks_2D_Coarse";
                elseif contains(fname,"2D_Watershed")
                    dtype = "Surface_FromStacks_2D_Watershed";
                else
                    dtype = "Surface_FromStacks_2D";
                end

            else

                if contains(fname,"Tracking_Watershed")
                    dtype = "Tracking_Watershed";
                elseif contains(fname,"Tracking")
                    dtype = "Tracking";
                elseif contains(fname,"Watershed")
                    dtype = "Watershed";
                elseif contains(fname,"Coarse")
                    dtype = "Coarse";
                else
                    dtype = "Binary";
                end

            end

            %% ---- genotype pairing ----
            % How replicates are grouped, and how the two channels are
            % treated, depends on the culture type:
            %   - Monoculture: one genotype, one "mono" series.
            %   - Protein addition: the red+green pair is treated as ONE
            %     combined replicate series (channel = "mono") rather
            %     than plotted as separate red/green traces -- here the
            %     two channels aren't competing genotypes.
            %   - Anything else (Coculture): red/green plotted as
            %     separate channels of the same genotype-pair group;
            %     files tagged neither 488 nor 561 are skipped.

            if cultureName == "Monoculture"

                genoKey = matlab.lang.makeValidName(redName);
                channel = "mono";

            elseif cultureName == "Protein addition"

                pairName = redName + "_" + greenName;
                genoKey = matlab.lang.makeValidName(pairName);
                channel = "mono";

            else

                pairName = redName + "_" + greenName;
                genoKey = matlab.lang.makeValidName(pairName);

                if contains(fname,"561")
                    channel = "red";
                elseif contains(fname,"488")
                    channel = "green";
                else
                    continue;
                end

            end

            %% ---- Read CSV ----

            Tall = readtable(fullfile(summaryPath,fname));

            if isempty(Tall)
                continue;
            end

            % Captured BEFORE any row filtering below, so a frame that
            % was imaged but ends up with zero qualifying objects still
            % gets an explicit (zero) value rather than being dropped.
            measuredFrames = unique(Tall.T);
            T = Tall;

            %% ---- Object size filter ----
            if contains(dtype,"3D")
                if ismember("Volume",string(T.Properties.VariableNames))
                    T = T(T.Volume >= CELL_VOLUME,:);
                end
            else
                if ismember("Area",string(T.Properties.VariableNames))
                    T = T(T.Area >= CELL_AREA,:);
                end
            end

            %% ---- Track filtering controls ----
            if contains(dtype,"FromStacks") && contains(dtype,"Tracking")
                T = T(T.TrackDisplacement < STACK_TRACK_MAX_DISP & ...
                      T.TrackDuration > STACK_TRACK_MIN_DUR,:);
            elseif contains(dtype,"Tracking")
                T = T(T.TrackDisplacement < ORIG_TRACK_MAX_DISP & ...
                      T.TrackDuration > ORIG_TRACK_MIN_DUR,:);
            end

            %% ---- 3D track match distance control ----
            if contains(dtype,"3D_Tracking")
                T = T(T.TrackingDistanceXY <= TRACK_MATCH_MAX_DIST | isnan(T.TrackingDistanceXY),:);
            end

            frames = measuredFrames;

            if isempty(frames)
                continue;
            end

            valsTotal = zeros(length(frames),1);
            valsCells = zeros(length(frames),1);
            valsNObj = zeros(length(frames),1);
            valsGT = zeros(length(frames),1);
            valsLT = zeros(length(frames),1);

            for t = 1:length(frames)

                idx2 = T.T == frames(t);

                % Middle_FromStacks rows span multiple Z-slices (the
                % +/-5-slice band around the estimated middle); divide by
                % how many of those slices actually contributed data at
                % this frame (counted from the UNFILTERED Tall, so the
                % slice count itself isn't affected by the size/track
                % filters above) to get a single-slice-equivalent value,
                % comparable to the single-slice Surface values.
                nSlicesHere = 1;
                if contains(dtype,"Middle_FromStacks") && ismember("Z",string(Tall.Properties.VariableNames))
                    idxAll = Tall.T == frames(t);
                    if contains(dtype,"3D")
                        nSlicesHere = numel(unique(round(Tall.Z(idxAll))));
                    else
                        nSlicesHere = numel(unique(Tall.Z(idxAll)));
                    end
                    if nSlicesHere < 1
                        nSlicesHere = 1;
                    end
                end

                if contains(dtype,"3D")
                    measureVals = T.Volume(idx2);
                    valsTotal(t) = sum(measureVals);
                    valsNObj(t) = sum(idx2);
                    valsGT(t) = sum(measureVals > THRESH_VOLUME);
                    valsLT(t) = sum(measureVals < THRESH_VOLUME);

                    if contains(dtype,"Middle_FromStacks")
                        valsTotal(t) = round(valsTotal(t) / nSlicesHere);
                        valsNObj(t) = round(valsNObj(t) / nSlicesHere);
                        valsGT(t) = round(valsGT(t) / nSlicesHere);
                        valsLT(t) = round(valsLT(t) / nSlicesHere);
                    end

                    valsCells(t) = valsTotal(t) / CELL_VOLUME;

                else
                    measureVals = T.Area(idx2);
                    valsTotal(t) = sum(measureVals);
                    valsNObj(t) = sum(idx2);
                    valsGT(t) = sum(measureVals > THRESH_AREA);
                    valsLT(t) = sum(measureVals < THRESH_AREA);

                    if contains(dtype,"Middle_FromStacks")
                        valsTotal(t) = round(valsTotal(t) / nSlicesHere);
                        valsNObj(t) = round(valsNObj(t) / nSlicesHere);
                        valsGT(t) = round(valsGT(t) / nSlicesHere);
                        valsLT(t) = round(valsLT(t) / nSlicesHere);
                    end

                    valsCells(t) = valsTotal(t) / CELL_AREA;
                end

            end

            timeVec = frames * timeStep;

            %% ---- Store replicate ----

            if ~isfield(genoData,genoKey)

                if cultureName == "Monoculture"
                    genoData.(genoKey).name = redName;
                else
                    genoData.(genoKey).name = redName + " + " + greenName;
                end

            end

            metricVals = {valsTotal, valsCells, valsNObj, valsGT, valsLT};

            for m = 1:length(metricNames)

                metricName = metricNames{m};

                if ~isfield(genoData.(genoKey),metricName)
                    genoData.(genoKey).(metricName) = struct();
                end

                if ~isfield(genoData.(genoKey).(metricName),dtype)
                    genoData.(genoKey).(metricName).(dtype) = struct();
                end

                if ~isfield(genoData.(genoKey).(metricName).(dtype),channel)
                    genoData.(genoKey).(metricName).(dtype).(channel).time = {};
                    genoData.(genoKey).(metricName).(dtype).(channel).vals = {};
                end

                % Each (date, M) contributing to this genoKey/dtype/
                % channel becomes one more replicate entry, averaged in
                % the next section.
                genoData.(genoKey).(metricName).(dtype).(channel).time{end+1} = timeVec;
                genoData.(genoKey).(metricName).(dtype).(channel).vals{end+1} = metricVals{m};

            end

        end
    end

    %% ===============================
    %% MEAN + SEM + SAVE
    %% ===============================
    % Averages every genoKey/metric/dtype/channel's accumulated
    % replicates across their union of timepoints, writes a per-replicate
    % + mean + SEM CSV, and plots a shaded-SEM mean curve overlaying all
    % channels present. Middle_FromStacks dtypes are skipped for direct
    % plotting -- they only exist to normalize the matching
    % Surface_FromStacks dtype below.

    genoNames = fieldnames(genoData);

    for g = 1:length(genoNames)

        key = genoNames{g};
        dataStruct = genoData.(key);
        metricFields = fieldnames(dataStruct);

        for mf = 1:length(metricFields)

            metricName = metricFields{mf};

            if metricName == "name"
                continue;
            end

            metricStruct = dataStruct.(metricName);
            dataTypes = fieldnames(metricStruct);

            for dt = 1:length(dataTypes)

                dtype = dataTypes{dt};

                if contains(dtype,"Middle_FromStacks")
                    continue;
                end

                channels = fieldnames(metricStruct.(dtype));

                outSubDir = fullfile(outDir,metricName,dtype);

                if ~exist(outSubDir,'dir')
                    mkdir(outSubDir);
                end

                fig = figure('Visible','off');
                hold on

                for ch = 1:length(channels)

                    chName = channels{ch};

                    times = metricStruct.(dtype).(chName).time;
                    vals = metricStruct.(dtype).(chName).vals;

                    % Align every replicate's own timepoints onto the
                    % union of all replicates' timepoints, leaving NaN
                    % wherever a replicate has no data at that time.
                    t = unique(vertcat(times{:}));
                    dataMat = nan(length(vals),length(t));

                    for r = 1:length(vals)
                        [tf,loc] = ismember(times{r},t);
                        dataMat(r,loc(tf)) = vals{r}(tf);
                    end

                    nRepHere = sum(~isnan(dataMat),1);
                    meanVals = mean(dataMat,1,'omitnan');
                    semVals = std(dataMat,0,1,'omitnan')./sqrt(nRepHere);
                    meanVals(nRepHere == 0) = nan;
                    semVals(nRepHere <= 1) = nan;   % can't estimate spread from a single replicate

                    outTable = table(t,'VariableNames',{'Time_min'});

                    for r = 1:size(dataMat,1)
                        outTable.("Rep"+r) = dataMat(r,:)';
                    end

                    outTable.Mean = meanVals';
                    outTable.SEM = semVals';

                    outFile = fullfile(outSubDir,key + "_" + metricName + "_" + dtype + "_" + chName + ".csv");
                    writeTableBlankNaN(outTable,outFile);

                    if chName == "red"
                        col = [1 0 0];
                    elseif chName == "green"
                        col = [0 0.6 0];
                    else
                        col = [0 0 0];
                    end

                    fill([t; flipud(t)], ...
                         [meanVals'-semVals'; flipud(meanVals'+semVals')], ...
                         col,'FaceAlpha',0.25,'EdgeColor','none');

                    plot(t,meanVals,'Color',col,'LineWidth',2);

                end

                xlabel("Time (minutes)");

                if metricName == "TotalMeasure"
                    if contains(dtype,"3D")
                        ylabel("Total Volume");
                    else
                        ylabel("Total Area");
                    end
                elseif metricName == "EstimatedCells"
                    ylabel("Estimated Number of Cells");
                elseif metricName == "TotalObjects"
                    ylabel("Total Number of Objects");
                elseif metricName == "NumObjectsGTThresh"
                    ylabel("Number of Objects > Threshold");
                elseif metricName == "NumObjectsLTThresh"
                    ylabel("Number of Objects < Threshold");
                end

                title(dataStruct.name + " " + metricName + " " + dtype);

                saveas(fig,fullfile(outSubDir,key + "_" + metricName + "_" + dtype + ".png"));
                close(fig)

                %% ===============================
                %% NORMALIZE SURFACE_FROMSTACKS TO MIDDLE_FROMSTACKS
                %% ===============================
                % For each Surface_FromStacks dtype, divide its mean by
                % the matching Middle_FromStacks dtype's mean at every
                % timepoint, so trends driven by overall biofilm growth
                % (which the middle band should also track) are factored
                % out, isolating surface-specific enrichment/depletion.
                % SEM on the ratio uses standard relative-error-in-
                % quadrature propagation.

                if contains(dtype,"Surface_FromStacks")

                    middleType = replace(dtype,"Surface_FromStacks","Middle_FromStacks");

                    if isfield(metricStruct,middleType)

                        middleStruct = metricStruct.(middleType);
                        middleChannels = fieldnames(middleStruct);

                        normType = replace(dtype,"Surface_FromStacks","SurfaceNormToMiddle");
                        normOutDir = fullfile(outDir,metricName,normType);

                        if ~exist(normOutDir,'dir')
                            mkdir(normOutDir);
                        end

                        figNorm = figure('Visible','off');
                        hold on

                        for ch = 1:length(channels)

                            chName = channels{ch};

                            if ~ismember(chName,middleChannels)
                                continue;
                            end

                            surfTimes = metricStruct.(dtype).(chName).time;
                            surfVals = metricStruct.(dtype).(chName).vals;

                            midTimes = metricStruct.(middleType).(chName).time;
                            midVals = metricStruct.(middleType).(chName).vals;

                            tNorm = unique([vertcat(surfTimes{:}); vertcat(midTimes{:})]);

                            surfMat = nan(length(surfVals),length(tNorm));
                            for r = 1:length(surfVals)
                                [tf,loc] = ismember(surfTimes{r},tNorm);
                                surfMat(r,loc(tf)) = surfVals{r}(tf);
                            end

                            midMat = nan(length(midVals),length(tNorm));
                            for r = 1:length(midVals)
                                [tf,loc] = ismember(midTimes{r},tNorm);
                                midMat(r,loc(tf)) = midVals{r}(tf);
                            end

                            nSurfHere = sum(~isnan(surfMat),1);
                            nMidHere = sum(~isnan(midMat),1);

                            surfMean = mean(surfMat,1,'omitnan');
                            surfSEM = std(surfMat,0,1,'omitnan')./sqrt(nSurfHere);
                            midMean = mean(midMat,1,'omitnan');
                            midSEM = std(midMat,0,1,'omitnan')./sqrt(nMidHere);

                            surfMean(nSurfHere == 0) = nan;
                            surfSEM(nSurfHere <= 1) = nan;
                            midMean(nMidHere == 0) = nan;
                            midSEM(nMidHere <= 1) = nan;

                            % Guard against exact-zero division; genuinely
                            % invalid (NaN) points are re-NaN-ed below
                            % regardless, but a true midMean == 0 with
                            % otherwise valid data would silently yield
                            % normMean == 0 rather than NaN/Inf.
                            midDenom = midMean;
                            midDenom(midDenom == 0) = 1;

                            normMean = surfMean ./ midDenom;
                            normSEM = normMean .* ...
                                sqrt((surfSEM./surfMean).^2 + ...
                                     (midSEM./midDenom).^2);
                            normMean(isnan(surfMean) | isnan(midMean)) = nan;
                            normSEM(isnan(normMean)) = nan;

                            outNormTable = table(tNorm,'VariableNames',{'Time_min'});
                            outNormTable.SurfaceMean = surfMean';
                            outNormTable.SurfaceSEM = surfSEM';
                            outNormTable.MiddleMean = midMean';
                            outNormTable.MiddleSEM = midSEM';
                            outNormTable.Mean = normMean';
                            outNormTable.SEM = normSEM';

                            outNormFile = fullfile(normOutDir,key + "_" + metricName + "_" + normType + "_" + chName + ".csv");
                            writeTableBlankNaN(outNormTable,outNormFile);

                            if chName == "red"
                                col = [1 0 0];
                            elseif chName == "green"
                                col = [0 0.6 0];
                            else
                                col = [0 0 0];
                            end

                            fill([tNorm; flipud(tNorm)], ...
                                 [normMean'-normSEM'; flipud(normMean'+normSEM')], ...
                                 col,'FaceAlpha',0.25,'EdgeColor','none');

                            plot(tNorm,normMean,'Color',col,'LineWidth',2);

                        end

                        xlabel("Time (minutes)");

                        if metricName == "TotalMeasure"
                            if contains(dtype,"3D")
                                ylabel("Normalized Volume");
                            else
                                ylabel("Normalized Area");
                            end
                        elseif metricName == "EstimatedCells"
                            ylabel("Normalized Estimated Cells");
                        elseif metricName == "TotalObjects"
                            ylabel("Normalized Total Objects");
                        elseif metricName == "NumObjectsGTThresh"
                            ylabel("Normalized Objects > Threshold");
                        elseif metricName == "NumObjectsLTThresh"
                            ylabel("Normalized Objects < Threshold");
                        end

                        title(dataStruct.name + " " + metricName + " " + normType);

                        saveas(figNorm,fullfile(normOutDir,key + "_" + metricName + "_" + normType + ".png"));
                        close(figNorm)

                    end

                end

            end
        end
    end
end

%% ===============================
%% STACK PROFILE PLOTS
%% ===============================
% Independent pass over the RAW, un-sliced per-Z-slice tables (base
% 2D/3D, Watershed, Coarse, Tracking) -- these show how each metric
% varies WITH Z at a given timepoint, so they're plotted per individual
% (date, M) replicate rather than averaged like the surface/middle
% summaries above.
%
% Note: the Tracking filter below always applies STACK_TRACK_MAX_DISP/
% STACK_TRACK_MIN_DUR here, whereas the COLLECT DATA section above
% applies ORIG_TRACK_* to these same "Tracking"/"Tracking_Watershed"
% dtype labels and reserves STACK_TRACK_* for the FromStacks-derived
% tracking tables -- worth confirming this difference is intentional.

for c = 1:length(cultureTypes)

    cultureName = cultureTypes(c);
    culturePath = fullfile(parentPath,cultureName);
    outDir = fullfile(culturePath,"StackProfiles");

    if ~exist(outDir,'dir')
        mkdir(outDir);
    end

    dateFolders = dir(culturePath);
    dateFolders = dateFolders([dateFolders.isdir]);
    dateFolders = dateFolders(~ismember({dateFolders.name},{'.','..','SurfacePlots','StackProfiles'}));

    for d = 1:length(dateFolders)

        dateName = string(dateFolders(d).name);
        summaryPath = fullfile(culturePath,dateName,"Summary");

        metaIdx = find(dates_meta == dateName);
        if isempty(metaIdx)
            continue;
        end

        Mcols = startsWith(meta.Properties.VariableNames,"M");
        XYcells = meta{metaIdx,Mcols};

        stackFiles = dir(fullfile(summaryPath,"M*.csv"));

        monoData = struct();
        cocoData = struct();

        for f = 1:length(stackFiles)

            fname = stackFiles(f).name;

            if contains(fname,"Surface") || contains(fname,"Middle") || contains(fname,"SurfaceSliceSummary")
                continue;
            end

            tok = regexp(fname,"M(\d+)_",'tokens','once');
            if isempty(tok)
                continue;
            end
            M = str2double(tok{1});

            if M > length(XYcells)
                continue;
            end

            cellVal = XYcells{M};
            if isempty(cellVal)
                continue;
            end

            parts = split(string(cellVal),"_");
            redCode = str2double(parts(1));
            greenCode = str2double(parts(2));

            redIdx = find(geno.Var1 == redCode);
            greenIdx = find(geno.Var1 == greenCode);

            redName = string(geno.Var2(redIdx));
            greenName = string(geno.Var2(greenIdx));

            %% ---- stack dtype ----
            % Anchored regexes (^M\d+_3D(_488|_561)?\.csv$ etc.) are
            % needed for the bare "3D"/"2D" cases specifically, since a
            % loose contains(fname,"3D") would also match "3D_Watershed"
            % etc. -- checked last, after the more specific substrings.

            if contains(fname,"3D_Tracking_Watershed")
                dtype = "3D_Tracking_Watershed";
            elseif contains(fname,"3D_Tracking")
                dtype = "3D_Tracking";
            elseif contains(fname,"2D_Tracking_Watershed")
                dtype = "2D_Tracking_Watershed";
            elseif contains(fname,"2D_Tracking")
                dtype = "2D_Tracking";
            elseif contains(fname,"3D_Coarse")
                dtype = "3D_Coarse";
            elseif contains(fname,"3D_Watershed")
                dtype = "3D_Watershed";
            elseif ~isempty(regexp(fname,"^M\d+_3D(_488|_561)?\.csv$",'once'))
                dtype = "3D";
            elseif contains(fname,"2D_Coarse")
                dtype = "2D_Coarse";
            elseif contains(fname,"2D_Watershed")
                dtype = "2D_Watershed";
            elseif ~isempty(regexp(fname,"^M\d+_2D(_488|_561)?\.csv$",'once'))
                dtype = "2D";
            else
                continue;
            end
            dtype = matlab.lang.makeValidName(dtype);

            %% ---- channel / genotype key ----
            % Unlike the COLLECT DATA section, the position number (M) is
            % baked directly into the storage key for Protein-addition
            % and Coculture data here (via pairKey below for Coculture),
            % since stack profiles are inherently per-individual-
            % replicate, not averaged across positions/dates.

            if cultureName == "Monoculture"
                genoKey = matlab.lang.makeValidName(redName);
                channel = "mono";

            elseif cultureName == "Protein addition"
                genoKey = matlab.lang.makeValidName(redName + "_" + greenName + "_M" + M);
                channel = "mono";

            else
                genoKey = matlab.lang.makeValidName(redName + "_" + greenName);

                if contains(fname,"561")
                    channel = "red";
                elseif contains(fname,"488")
                    channel = "green";
                else
                    continue;
                end
            end

            %% ---- Read CSV ----

            T = readtable(fullfile(summaryPath,fname));

            %% ---- Object size filter ----
            if contains(dtype,"3D")
                if ismember("Volume",string(T.Properties.VariableNames))
                    T = T(T.Volume >= CELL_VOLUME,:);
                end
            else
                if ismember("Area",string(T.Properties.VariableNames))
                    T = T(T.Area >= CELL_AREA,:);
                end
            end

            if isempty(T)
                continue;
            end

            %% ---- Track filtering controls ----
            if contains(dtype,"Tracking")
                T = T(T.TrackDisplacement < STACK_TRACK_MAX_DISP & ...
                      T.TrackDuration > STACK_TRACK_MIN_DUR,:);
            end

            %% ---- 3D track match distance control ----
            if contains(dtype,"3D_Tracking")
                T = T(T.TrackingDistanceXY <= TRACK_MATCH_MAX_DIST | isnan(T.TrackingDistanceXY),:);
            end

            if isempty(T)
                continue;
            end

            % Unlike COLLECT DATA's measuredFrames (captured before
            % filtering), timeVals here is taken AFTER filtering, so a
            % timepoint with no rows surviving the filters above simply
            % won't appear in the profile at all.
            timeVals = unique(T.T);
            if isempty(timeVals)
                continue;
            end

            metricTables = struct();

            for m = 1:length(metricNames)

                metricName = metricNames{m};
                rowsOut = table();

                for tt = 1:length(timeVals)

                    thisT = timeVals(tt);
                    idxT = T.T == thisT;

                    if contains(dtype,"3D")
                        zVals = unique(round(T.Z(idxT)));
                    else
                        zVals = unique(T.Z(idxT));
                    end

                    for zz = 1:length(zVals)

                        thisZ = zVals(zz);

                        if contains(dtype,"3D")
                            idxTZ = idxT & round(T.Z) == thisZ;
                            measureVals = T.Volume(idxTZ);
                        else
                            idxTZ = idxT & T.Z == thisZ;
                            measureVals = T.Area(idxTZ);
                        end

                        if metricName == "TotalMeasure"
                            thisVal = sum(measureVals);
                        elseif metricName == "EstimatedCells"
                            if contains(dtype,"3D")
                                thisVal = sum(measureVals) / CELL_VOLUME;
                            else
                                thisVal = sum(measureVals) / CELL_AREA;
                            end
                        elseif metricName == "TotalObjects"
                            thisVal = sum(idxTZ);
                        elseif metricName == "NumObjectsGTThresh"
                            if contains(dtype,"3D")
                                thisVal = sum(measureVals > THRESH_VOLUME);
                            else
                                thisVal = sum(measureVals > THRESH_AREA);
                            end
                        elseif metricName == "NumObjectsLTThresh"
                            if contains(dtype,"3D")
                                thisVal = sum(measureVals < THRESH_VOLUME);
                            else
                                thisVal = sum(measureVals < THRESH_AREA);
                            end
                        end

                        % Long-format: one row per (timepoint, Z-slice).
                        rowsOut = [rowsOut; table(thisT,thisZ,thisVal,'VariableNames',{'T','Z','Value'})];
                    end
                end

                metricTables.(metricName) = rowsOut;
            end

            %% ---- Store ----

            if cultureName == "Monoculture" || cultureName == "Protein addition"

                if ~isfield(monoData,genoKey)
                    monoData.(genoKey) = struct();
                    if cultureName == "Protein addition"
                        monoData.(genoKey).name = redName + " + " + greenName;
                    else
                        monoData.(genoKey).name = redName;
                    end
                end
                if ~isfield(monoData.(genoKey),dtype)
                    monoData.(genoKey).(dtype) = struct();
                end
                monoData.(genoKey).(dtype).M = M;

                for m = 1:length(metricNames)
                    metricName = metricNames{m};
                    monoData.(genoKey).(dtype).(metricName) = metricTables.(metricName);
                end

            else

                % pairKey bakes in date + M, so every Coculture replicate
                % also gets its own separate (un-averaged) profile entry.
                pairKey = matlab.lang.makeValidName("D" + dateName + "_M" + M + "_" + genoKey);

                if ~isfield(cocoData,pairKey)
                    cocoData.(pairKey) = struct();
                    cocoData.(pairKey).dateName = dateName;
                    cocoData.(pairKey).M = M;
                    cocoData.(pairKey).pairName = redName + " + " + greenName;
                    cocoData.(pairKey).redName = redName;
                    cocoData.(pairKey).greenName = greenName;
                end

                if ~isfield(cocoData.(pairKey),dtype)
                    cocoData.(pairKey).(dtype) = struct();
                end

                for m = 1:length(metricNames)
                    metricName = metricNames{m};
                    cocoData.(pairKey).(dtype).(channel).(metricName) = metricTables.(metricName);
                end

            end

        end

        %% ---- Save monoculture stack profiles ----
        % One line per timepoint (Value vs Z), shaded light-to-dark gray
        % from earliest to latest timepoint, so the Z-profile's evolution
        % over time is visible in a single plot.

        monoKeys = fieldnames(monoData);

        for g = 1:length(monoKeys)

            genoKey = monoKeys{g};
            dataStruct = monoData.(genoKey);
            dataFields = fieldnames(dataStruct);

            for df = 1:length(dataFields)

                dtype = dataFields{df};
                if dtype == "name"
                    continue;
                end

                M = dataStruct.(dtype).M;

                for m = 1:length(metricNames)

                    metricName = metricNames{m};
                    plotTable = dataStruct.(dtype).(metricName);

                    if isempty(plotTable)
                        continue;
                    end

                    outSubDir = fullfile(outDir,metricName,dtype);
                    if ~exist(outSubDir,'dir')
                        mkdir(outSubDir);
                    end

                    fig = figure('Visible','off');
                    hold on

                    tUnique = unique(plotTable.T);
                    nT = length(tUnique);

                    for tt = 1:nT

                        thisT = tUnique(tt);
                        idxT = plotTable.T == thisT;
                        Zplot = plotTable.Z(idxT);
                        Vplot = plotTable.Value(idxT);
                        [Zplot,ord] = sort(Zplot);
                        Vplot = Vplot(ord);

                        shade = 0.85 - 0.75*(tt-1)/max(nT-1,1);
                        col = [shade shade shade];

                        plot(Zplot,Vplot,'Color',col,'LineWidth',1.5);
                    end

                    xlabel("Z Position");

                    if metricName == "TotalMeasure"
                        if contains(dtype,"3D")
                            ylabel("Total Volume");
                        else
                            ylabel("Total Area");
                        end
                    elseif metricName == "EstimatedCells"
                        ylabel("Estimated Number of Cells");
                    elseif metricName == "TotalObjects"
                        ylabel("Total Number of Objects");
                    elseif metricName == "NumObjectsGTThresh"
                        ylabel("Number of Objects > Threshold");
                    elseif metricName == "NumObjectsLTThresh"
                        ylabel("Number of Objects < Threshold");
                    end

                    title(dataStruct.name + " " + dateName + " M" + M + " " + metricName + " " + dtype);

                    outBase = dateName + "_M" + M + "_" + metricName + "_" + dtype;
                    outCsv = fullfile(outSubDir,outBase + ".csv");
                    outPng = fullfile(outSubDir,outBase + ".png");

                    plotTable.Channel = repmat("mono",height(plotTable),1);
                    plotTable = movevars(plotTable,'Channel','Before','T');

                    writeTableBlankNaN(plotTable,outCsv)
                    saveas(fig,outPng)
                    close(fig)

                end
            end
        end

        %% ---- Save coculture stack profiles (true combined figures) ----
        % Side-by-side two-panel figure: green channel's Z-profile on the
        % left, red channel's on the right, same light-to-dark-by-time
        % gradient idea but tinted green/red respectively, for direct
        % visual comparison of the two genotypes at the same position.

        pairKeys = fieldnames(cocoData);

        for p = 1:length(pairKeys)

            pairKey = pairKeys{p};
            pairStruct = cocoData.(pairKey);
            pairFields = fieldnames(pairStruct);

            for pf = 1:length(pairFields)

                dtype = pairFields{pf};

                if ismember(dtype,["dateName","M","pairName","redName","greenName"])
                    continue;
                end

                for m = 1:length(metricNames)

                    metricName = metricNames{m};

                    hasGreen = isfield(pairStruct.(dtype),"green") && isfield(pairStruct.(dtype).green,metricName);
                    hasRed = isfield(pairStruct.(dtype),"red") && isfield(pairStruct.(dtype).red,metricName);

                    if ~hasGreen && ~hasRed
                        continue;
                    end

                    if hasGreen
                        greenTable = pairStruct.(dtype).green.(metricName);
                    else
                        greenTable = table();
                    end

                    if hasRed
                        redTable = pairStruct.(dtype).red.(metricName);
                    else
                        redTable = table();
                    end

                    if isempty(greenTable) && isempty(redTable)
                        continue;
                    end

                    outSubDir = fullfile(outDir,metricName,dtype);
                    if ~exist(outSubDir,'dir')
                        mkdir(outSubDir);
                    end

                    fig = figure('Visible','off');

                    subplot(1,2,1)
                    hold on
                    if ~isempty(greenTable)

                        tUnique = unique(greenTable.T);
                        nT = length(tUnique);

                        for tt = 1:nT

                            thisT = tUnique(tt);
                            idxT = greenTable.T == thisT;
                            Zplot = greenTable.Z(idxT);
                            Vplot = greenTable.Value(idxT);
                            [Zplot,ord] = sort(Zplot);
                            Vplot = Vplot(ord);

                            frac = 0.9 - 0.75*(tt-1)/max(nT-1,1);
                            col = [0.15 + 0.55*frac, 0.9 - 0.55*frac, 0.15 + 0.55*frac];

                            plot(Zplot,Vplot,'Color',col,'LineWidth',1.5);
                        end
                    end
                    xlabel("Z Position");
                    if metricName == "TotalMeasure"
                        if contains(dtype,"3D")
                            ylabel("Total Volume");
                        else
                            ylabel("Total Area");
                        end
                    elseif metricName == "EstimatedCells"
                        ylabel("Estimated Number of Cells");
                    elseif metricName == "TotalObjects"
                        ylabel("Total Number of Objects");
                    elseif metricName == "NumObjectsGTThresh"
                        ylabel("Number of Objects > Threshold");
                    elseif metricName == "NumObjectsLTThresh"
                        ylabel("Number of Objects < Threshold");
                    end
                    title(pairStruct.greenName)

                    subplot(1,2,2)
                    hold on
                    if ~isempty(redTable)

                        tUnique = unique(redTable.T);
                        nT = length(tUnique);

                        for tt = 1:nT

                            thisT = tUnique(tt);
                            idxT = redTable.T == thisT;
                            Zplot = redTable.Z(idxT);
                            Vplot = redTable.Value(idxT);
                            [Zplot,ord] = sort(Zplot);
                            Vplot = Vplot(ord);

                            frac = 0.9 - 0.75*(tt-1)/max(nT-1,1);
                            col = [0.9 - 0.55*frac, 0.15 + 0.55*frac, 0.15 + 0.55*frac];

                            plot(Zplot,Vplot,'Color',col,'LineWidth',1.5);
                        end
                    end
                    xlabel("Z Position");
                    if metricName == "TotalMeasure"
                        if contains(dtype,"3D")
                            ylabel("Total Volume");
                        else
                            ylabel("Total Area");
                        end
                    elseif metricName == "EstimatedCells"
                        ylabel("Estimated Number of Cells");
                    elseif metricName == "TotalObjects"
                        ylabel("Total Number of Objects");
                    elseif metricName == "NumObjectsGTThresh"
                        ylabel("Number of Objects > Threshold");
                    elseif metricName == "NumObjectsLTThresh"
                        ylabel("Number of Objects < Threshold");
                    end
                    title(pairStruct.redName)

                    sgtitle(pairStruct.pairName + " " + pairStruct.dateName + " M" + pairStruct.M + " " + metricName + " " + dtype)

                    outBase = pairStruct.dateName + "_M" + pairStruct.M + "_" + metricName + "_" + dtype;
                    outCsv = fullfile(outSubDir,outBase + ".csv");
                    outPng = fullfile(outSubDir,outBase + ".png");

                    comboTable = table();

                    if ~isempty(greenTable)
                        greenTable.Channel = repmat("green",height(greenTable),1);
                        greenTable = movevars(greenTable,'Channel','Before','T');
                        comboTable = [comboTable; greenTable];
                    end

                    if ~isempty(redTable)
                        redTable.Channel = repmat("red",height(redTable),1);
                        redTable = movevars(redTable,'Channel','Before','T');
                        comboTable = [comboTable; redTable];
                    end

                    writeTableBlankNaN(comboTable,outCsv)
                    saveas(fig,outPng)
                    close(fig)

                end
            end
        end

    end
end

%% ===============================
%% COCULTURE RED-TO-GREEN DISTANCE PLOTS
%% ===============================
% Only ever produces output when "Coculture" is included in cultureTypes
% above (skipped entirely otherwise). For each Coculture replicate
% (date, M), computes the mean nearest-neighbor XY distance from every
% qualifying red (561) object to the nearest qualifying green (488)
% object at each timepoint, then averages that distance across
% replicates sharing the same genotype pair.

distData = struct();

for c = 1:length(cultureTypes)

    cultureName = cultureTypes(c);

    if cultureName ~= "Coculture"
        continue;
    end

    culturePath = fullfile(parentPath,cultureName);
    outDir = fullfile(culturePath,"RedToGreenDistance");

    if ~exist(outDir,'dir')
        mkdir(outDir);
    end

    dateFolders = dir(culturePath);
    dateFolders = dateFolders([dateFolders.isdir]);
    dateFolders = dateFolders(~ismember({dateFolders.name},{'.','..','SurfacePlots','StackProfiles','RedToGreenDistance'}));

    for d = 1:length(dateFolders)

        dateName = string(dateFolders(d).name);
        summaryPath = fullfile(culturePath,dateName,"Summary");

        metaIdx = find(dates_meta == dateName);

        if isempty(metaIdx)
            continue;
        end

        timeStep = meta.TimeStep(metaIdx);

        Mcols = startsWith(meta.Properties.VariableNames,"M");
        XYcells = meta{metaIdx,Mcols};

        csvFiles = dir(fullfile(summaryPath,"M*.csv"));

        % Reset per date: repData only ever pairs up red/green files
        % from THIS date's positions before the distance calc below.
        repData = struct();

        for f = 1:length(csvFiles)

            fname = csvFiles(f).name;

            if ~contains(fname,"Surface")
                continue;
            end

            if contains(fname,"Middle") || contains(fname,"SurfaceSliceSummary")
                continue;
            end

            tok = regexp(fname,"M(\d+)_",'tokens','once');
            if isempty(tok)
                continue;
            end
            M = str2double(tok{1});

            if M > length(XYcells)
                continue;
            end

            cellVal = XYcells{M};

            if isempty(cellVal)
                continue;
            end

            parts = split(string(cellVal),"_");
            redCode = str2double(parts(1));
            greenCode = str2double(parts(2));

            redIdx = find(geno.Var1 == redCode);
            greenIdx = find(geno.Var1 == greenCode);

            redName = string(geno.Var2(redIdx));
            greenName = string(geno.Var2(greenIdx));

            %% ---- dtype ----

            if contains(fname,"Surface_FromStacks")

                if contains(fname,"3D_Tracking_Watershed")
                    dtype = "Surface_FromStacks_3D_Tracking_Watershed";
                elseif contains(fname,"3D_Tracking")
                    dtype = "Surface_FromStacks_3D_Tracking";
                elseif contains(fname,"2D_Tracking_Watershed")
                    dtype = "Surface_FromStacks_2D_Tracking_Watershed";
                elseif contains(fname,"2D_Tracking")
                    dtype = "Surface_FromStacks_2D_Tracking";
                elseif contains(fname,"3D_Coarse")
                    dtype = "Surface_FromStacks_3D_Coarse";
                elseif contains(fname,"3D_Watershed")
                    dtype = "Surface_FromStacks_3D_Watershed";
                elseif contains(fname,"3D")
                    dtype = "Surface_FromStacks_3D";
                elseif contains(fname,"2D_Coarse")
                    dtype = "Surface_FromStacks_2D_Coarse";
                elseif contains(fname,"2D_Watershed")
                    dtype = "Surface_FromStacks_2D_Watershed";
                else
                    dtype = "Surface_FromStacks_2D";
                end

            else

                if contains(fname,"Tracking_Watershed")
                    dtype = "Tracking_Watershed";
                elseif contains(fname,"Tracking")
                    dtype = "Tracking";
                elseif contains(fname,"Watershed")
                    dtype = "Watershed";
                elseif contains(fname,"Coarse")
                    dtype = "Coarse";
                else
                    dtype = "Binary";
                end

            end

            if contains(fname,"561")
                channel = "red";
            elseif contains(fname,"488")
                channel = "green";
            else
                continue;
            end

            pairKey = matlab.lang.makeValidName(redName + "_" + greenName);
            repKey = matlab.lang.makeValidName("D" + dateName + "_M" + M);

            if ~isfield(repData,repKey)
                repData.(repKey) = struct();
                repData.(repKey).pairName = redName + " + " + greenName;
                repData.(repKey).pairKey = pairKey;
                repData.(repKey).dateName = dateName;
                repData.(repKey).M = M;
            end

            if ~isfield(repData.(repKey),dtype)
                repData.(repKey).(dtype) = struct();
            end

            Tall = readtable(fullfile(summaryPath,fname));

            if isempty(Tall)
                continue;
            end

            T = Tall;

            %% ---- Object size filter ----
            if contains(dtype,"3D")
                if ismember("Volume",string(T.Properties.VariableNames))
                    T = T(T.Volume >= CELL_VOLUME,:);
                end
            else
                if ismember("Area",string(T.Properties.VariableNames))
                    T = T(T.Area >= CELL_AREA,:);
                end
            end

            %% ---- existing tracking filters ----
            if contains(dtype,"FromStacks") && contains(dtype,"Tracking")
                T = T(T.TrackDisplacement < STACK_TRACK_MAX_DISP & ...
                      T.TrackDuration > STACK_TRACK_MIN_DUR,:);
            elseif contains(dtype,"Tracking")
                T = T(T.TrackDisplacement < ORIG_TRACK_MAX_DISP & ...
                      T.TrackDuration > ORIG_TRACK_MIN_DUR,:);
            end

            %% ---- 3D track match distance control ----
            if contains(dtype,"3D_Tracking")
                T = T(T.TrackingDistanceXY <= TRACK_MATCH_MAX_DIST | isnan(T.TrackingDistanceXY),:);
            end

            % Filtered per-channel table, plus the UNFILTERED set of
            % frames actually imaged (same "was it imaged at all" role as
            % measuredFrames earlier).
            repData.(repKey).(dtype).(channel) = T;
            repData.(repKey).(dtype).(channel + "Frames") = unique(Tall.T);

        end

        repKeys = fieldnames(repData);

        for r = 1:length(repKeys)

            repKey = repKeys{r};
            pairKey = repData.(repKey).pairKey;
            pairName = repData.(repKey).pairName;

            dtypeFields = fieldnames(repData.(repKey));

            for df = 1:length(dtypeFields)

                dtype = dtypeFields{df};

                if ismember(dtype,["pairName","pairKey","dateName","M"])
                    continue;
                end

                if ~isfield(repData.(repKey).(dtype),"red") || ~isfield(repData.(repKey).(dtype),"green")
                    continue;
                end

                Tred = repData.(repKey).(dtype).red;
                Tgreen = repData.(repKey).(dtype).green;
                redFrames = repData.(repKey).(dtype).redFrames;
                greenFrames = repData.(repKey).(dtype).greenFrames;

                allT = unique([redFrames; greenFrames]);
                if isempty(allT)
                    continue;
                end

                distVals = nan(length(allT),1);
                nRedVals = nan(length(allT),1);
                nGreenVals = nan(length(allT),1);

                for tt = 1:length(allT)

                    thisT = allT(tt);

                    redMeasured = ismember(thisT,redFrames);
                    greenMeasured = ismember(thisT,greenFrames);

                    % A frame not actually imaged in one of the two
                    % channels leaves a true NaN gap here, distinct from
                    % the "imaged but nothing qualified" zero cases below.
                    if ~redMeasured || ~greenMeasured
                        distVals(tt) = nan;
                        nRedVals(tt) = nan;
                        nGreenVals(tt) = nan;
                        continue;
                    end

                    redNow = Tred(Tred.T == thisT,:);
                    greenNow = Tgreen(Tgreen.T == thisT,:);

                    % Stricter size cutoff than the general noise filter,
                    % specific to what counts as a real neighbor target.
                    if contains(dtype,"3D")
                        greenMeasure = greenNow.Volume;
                        greenKeep = greenMeasure > NEIGHBOR_THRESH_VOLUME;
                    else
                        greenMeasure = greenNow.Area;
                        greenKeep = greenMeasure > NEIGHBOR_THRESH_AREA;
                    end

                    greenNow = greenNow(greenKeep,:);

                    nGreenVals(tt) = height(greenNow);

                    if isempty(greenNow)
                        distVals(tt) = 0;
                        nRedVals(tt) = 0;
                        continue;
                    end

                    % For Tracking dtypes, only measure from red objects
                    % whose track started at this frame. Note:
                    % TRACK_START in Tred was already converted to
                    % 1-based when it was written by
                    % CompileCultureTracking...m, so comparing it against
                    % (thisT - 1) here is worth double-checking against
                    % that convention.
                    if contains(dtype,"Tracking")

                        if ismember("TRACK_START", string(Tred.Properties.VariableNames))
                            redUse = redNow(redNow.TRACK_START == (thisT - 1),:);
                        else
                            redUse = redNow([],:);
                        end

                    else
                        redUse = redNow;
                    end

                    nRedVals(tt) = height(redUse);

                    if isempty(redUse)
                        distVals(tt) = 0;
                        continue;
                    end

                    redX = redUse.X;
                    redY = redUse.Y;
                    greenX = greenNow.X;
                    greenY = greenNow.Y;

                    D = (redX - greenX').^2 + (redY - greenY').^2;
                    minD = sqrt(min(D,[],2));

                    distVals(tt) = mean(minD,'omitnan');

                end

                timeVec = allT * timeStep;

                if ~isfield(distData,pairKey)
                    distData.(pairKey) = struct();
                    distData.(pairKey).name = pairName;
                end

                if ~isfield(distData.(pairKey),dtype)
                    distData.(pairKey).(dtype) = struct();
                    distData.(pairKey).(dtype).time = {};
                    distData.(pairKey).(dtype).dist = {};
                    distData.(pairKey).(dtype).nRed = {};
                    distData.(pairKey).(dtype).nGreen = {};
                end

                distData.(pairKey).(dtype).time{end+1} = timeVec;
                distData.(pairKey).(dtype).dist{end+1} = distVals;
                distData.(pairKey).(dtype).nRed{end+1} = nRedVals;
                distData.(pairKey).(dtype).nGreen{end+1} = nGreenVals;

            end
        end
    end
end

%% ===============================
%% SAVE RED-TO-GREEN DISTANCE PLOTS
%% ===============================
% Same align-to-union-of-timepoints + mean/SEM pattern as the SurfacePlots
% section, plus mean red/green object counts as extra diagnostic columns.
% Note outSubDir below hardcodes a "Coculture" folder name rather than
% using a cultureName variable (this section runs after the cultureTypes
% loop above, so none is in scope) -- harmless while distData is only
% ever populated when "Coculture" is in cultureTypes, but worth knowing
% it doesn't dynamically follow the culture type that produced the data.

pairKeys = fieldnames(distData);

for p = 1:length(pairKeys)

    pairKey = pairKeys{p};
    dataStruct = distData.(pairKey);
    dtypeFields = fieldnames(dataStruct);

    for df = 1:length(dtypeFields)

        dtype = dtypeFields{df};

        if dtype == "name"
            continue;
        end

        times = dataStruct.(dtype).time;
        distVals = dataStruct.(dtype).dist;
        nRedVals = dataStruct.(dtype).nRed;
        nGreenVals = dataStruct.(dtype).nGreen;

        t = unique(vertcat(times{:}));

        distMat = nan(length(distVals),length(t));
        nRedMat = nan(length(distVals),length(t));
        nGreenMat = nan(length(distVals),length(t));

        for r = 1:length(distVals)
            [tf,loc] = ismember(times{r},t);
            distMat(r,loc(tf)) = distVals{r}(tf);
            nRedMat(r,loc(tf)) = nRedVals{r}(tf);
            nGreenMat(r,loc(tf)) = nGreenVals{r}(tf);
        end

        nDistHere = sum(~isnan(distMat),1);
        meanVals = mean(distMat,1,'omitnan');
        semVals = std(distMat,0,1,'omitnan') ./ sqrt(nDistHere);
        meanVals(nDistHere == 0) = nan;
        semVals(nDistHere <= 1) = nan;

        meanNRed = mean(nRedMat,1,'omitnan');
        meanNGreen = mean(nGreenMat,1,'omitnan');

        outSubDir = fullfile(parentPath,"Coculture","RedToGreenDistance",dtype);

        if ~exist(outSubDir,'dir')
            mkdir(outSubDir);
        end

        outTable = table(t,'VariableNames',{'Time_min'});

        for r = 1:size(distMat,1)
            outTable.("Rep"+r+"_Distance") = distMat(r,:)';
        end
        for r = 1:size(nRedMat,1)
            outTable.("Rep"+r+"_NRed") = nRedMat(r,:)';
        end
        for r = 1:size(nGreenMat,1)
            outTable.("Rep"+r+"_NGreen") = nGreenMat(r,:)';
        end

        outTable.MeanDistance = meanVals';
        outTable.SEMDistance = semVals';
        outTable.MeanNRed = meanNRed';
        outTable.MeanNGreen = meanNGreen';

        outCsv = fullfile(outSubDir,pairKey + "_" + dtype + ".csv");
        writeTableBlankNaN(outTable,outCsv)

        fig = figure('Visible','off');
        hold on

        fill([t; flipud(t)], ...
             [meanVals'-semVals'; flipud(meanVals'+semVals')], ...
             [1 0 0],'FaceAlpha',0.25,'EdgeColor','none');

        plot(t,meanVals,'Color',[1 0 0],'LineWidth',2);

        xlabel("Time (minutes)");
        ylabel("Mean distance to nearest green object");
        title(dataStruct.name + " " + dtype + " Red-to-Green Distance");

        saveas(fig,fullfile(outSubDir,pairKey + "_" + dtype + ".png"));
        close(fig)

    end
end

function writeTableBlankNaN(T,outFile)
% Writes a table to CSV with NaN scalar numeric cells rendered as blank
% instead of the literal text "NaN", for cleaner reading downstream
% (e.g. in Excel or a plotting tool that would otherwise treat "NaN" as
% a string).
C = table2cell(T);
for ii = 1:numel(C)
    if isnumeric(C{ii}) && isscalar(C{ii}) && isnan(C{ii})
        C{ii} = '';
    end
end
writecell([T.Properties.VariableNames; C],outFile);
end
