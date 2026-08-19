%% SummarizeAndPlotSpikeInData.m
% Author: Hasreet Gill, Bassler Lab, HHMI/Princeton University
% Revised with Claude (Sonnet 5, Anthropic)

% Plotting/analysis stage for the spike-in attachment assay, run after
% SegmentSpikeInAttachment.m has written ImageSummary_All.xlsx. For a
% filtered subset of genotypes/replicates/dates, produces:
%   PART 1: one line plot (mean +/- SEM across replicates vs. time) per
%           genotype, per metric -- split into separate colored series
%           by induction Condition ("+ara -nspd" etc.) for genotypes
%           that have one, or by culture age (Hours) for genotypes that
%           don't.
%   PART 2: the transpose of Part 1 -- one line plot per Condition, per
%           metric, with all genotypes sharing that condition as
%           separate series.
%   PART 3: bar plots per genotype at fixed target times (barTimes):
%           an "hours-grouped" bar+line comparison across culture ages
%           for genotypes with no condition, and a per-condition bar
%           comparison for genotypes that have one.
%   PART 4: the transpose of Part 3's condition bars -- one bar plot per
%           condition, comparing all genotypes under that condition.
% Every PNG gets a matching CSV of the exact plotted values (mean, SEM,
% and every individual replicate value) via writeLineCSV/writeBarCSV.
%
% Two things worth knowing about the SETTINGS below:
%   - infile and the "041526" output subfolder tag are both hardcoded
%     literals that happen to match by convention (the exported summary
%     file is manually dated/copied, and the output folder tag is kept
%     in sync by hand) -- both need updating together when running this
%     against a newer export.
%   - barReplicateGroups lets the same genotype/condition/time cell be
%     satisfied by either replicates [1 2 3] or [4 5 6]; getBarReplicate-
%     ValuesAndError below tries both and keeps whichever triplet has
%     more actual (non-NaN) data, rather than assuming replicates always
%     come in one fixed set of three.

clear; clc; close all;

%% ================= SETTINGS =================

path = uigetdir(pwd, 'Choose folder with ImageSummary_All.xlsx');
if isequal(path,0)
    error('No folder selected.');
end
cd(path);

infile = "ImageSummary_All_041526.xlsx";

% Genotype allowlist -- leave empty to include every genotype present.
useGenotypes = [ ...
    % "luxOD61E"
    % "luxOD61E ∆mshA"
    % "luxOD61A ∆mshA"
    % "luxOD61A ∆mshA pbad-vpsT"
    % "luxOD61E ∆mshA pbad-vpsT"
];

% Genotype denylist, applied after the allowlist above.
excludeGenotypes = [ ...
    "luxOD61E"
];

useReplicates = [1 2 3 4 5 6];

% Whole dates to drop entirely (e.g. a failed/discarded run).
excludeDates = [ ...
    "040126"
    "040926"
];

metricsToPlot = [ ...
    "ProportionWithTracks"
    "ProportionGT1Track"
    "MeanFractionAreaFilledWithSpots"
    "MeanSpotsPerObject"
];

% Declared separately from metricsToPlot so the two lists COULD diverge,
% even though they're identical right now.
barMetrics = [ ...
    "ProportionWithTracks"
    "ProportionGT1Track"
    "MeanFractionAreaFilledWithSpots"
    "MeanSpotsPerObject"
];

barTimes = [30 60 120];

barReplicateGroups = { ...
    [1 2 3]
    [4 5 6]
};

%% ================= LOAD =================

T = readtable(infile, 'VariableNamingRule','preserve');

% Coerce to string explicitly -- readtable may otherwise infer these as
% categorical/cell, which behaves inconsistently with the == / ismember
% comparisons used throughout below.
vars = string(T.Properties.VariableNames);
if ismember("Genotype", vars),  T.Genotype  = string(T.Genotype); end
if ismember("Date", vars),      T.Date      = string(T.Date); end
if ismember("XY", vars),        T.XY        = string(T.XY); end
if ismember("ImageName", vars), T.ImageName = string(T.ImageName); end

%% ================= LABELS =================
% Derives a human-readable induction Condition label per row. Plain
% strains (no "pbad-vpsT" construct and not an "empty" vector control)
% get no condition label at all (blank "") -- only the inducible-vpsT
% construct and its empty-vector control have a meaningful arabinose/
% spermidine induction state to report. This blank-vs-labeled split is
% exactly what later sections use to decide whether a genotype's plots
% should be split by Condition or by Hours instead.

T.BaseGenotype = strtrim(T.Genotype);
T.Condition = strings(height(T),1);

for i = 1:height(T)
    if ~contains(T.Genotype(i), "pbad-vpsT") && ~contains(T.Genotype(i), "empty")
        T.Condition(i) = "";
    else
        if ~isnan(T.Arabinose(i)) && T.Arabinose(i) == 1
            araStr = "+ara";
        else
            araStr = "-ara";
        end

        if ~isnan(T.Nspd(i)) && T.Nspd(i) == 1
            nspdStr = "+nspd";
        else
            nspdStr = "-nspd";
        end

        T.Condition(i) = araStr + " " + nspdStr;
    end
end

%% ================= FILTER =================
% Each filter below is an independent AND-mask on T; order doesn't
% change the result, just narrows the row set further at each step.

if ~isempty(useGenotypes)
    T = T(ismember(T.BaseGenotype, useGenotypes), :);
end

if ~isempty(excludeGenotypes)
    T = T(~ismember(T.BaseGenotype, excludeGenotypes), :);
end

if ~isempty(useReplicates)
    T = T(ismember(T.Replicate, useReplicates), :);
end

if ~isempty(excludeDates)
    T = T(~ismember(T.Date, excludeDates), :);
end

if isempty(T)
    error('No rows remain after filtering.');
end

%% ================= TIME + CLEAN LABELS =================

% 1-based timepoint index -> real elapsed minutes, same convention used
% throughout this pipeline family.
T.Time_min = (T.T - 1) .* T.TimeStep_m;
T.Condition = strtrim(T.Condition);

%% ================= OUTPUT FOLDERS =================

rootOut = fullfile(path, "041526");
if ~exist(rootOut, 'dir'), mkdir(rootOut); end

metricDirs = strings(size(metricsToPlot));
comboDirs  = strings(size(metricsToPlot));

for m = 1:length(metricsToPlot)
    metric = metricsToPlot(m);

    metricDir = fullfile(rootOut, metric);
    if ~exist(metricDir, 'dir'), mkdir(metricDir); end
    metricDirs(m) = metricDir;

    comboDir = fullfile(rootOut, metric + "_AllGenotypes_ByCondition");
    if ~exist(comboDir, 'dir'), mkdir(comboDir); end
    comboDirs(m) = comboDir;
end

barGenoDirs = containers.Map;
barCondDirs = containers.Map;

for bm = 1:length(barMetrics)
    barMetric = barMetrics(bm);

    barRoot = fullfile(rootOut, "Bars_" + barMetric);
    if ~exist(barRoot, 'dir'), mkdir(barRoot); end

    genoDir = fullfile(barRoot, "PerGenotype");
    condDir = fullfile(barRoot, "PerCondition");

    if ~exist(genoDir, 'dir'), mkdir(genoDir); end
    if ~exist(condDir, 'dir'), mkdir(condDir); end

    for tt = 1:length(barTimes)
        d1 = fullfile(genoDir, sprintf('%dmin', barTimes(tt)));
        d2 = fullfile(condDir, sprintf('%dmin', barTimes(tt)));
        if ~exist(d1, 'dir'), mkdir(d1); end
        if ~exist(d2, 'dir'), mkdir(d2); end
    end

    barGenoDirs(char(barMetric)) = genoDir;
    barCondDirs(char(barMetric)) = condDir;
end

%% ================= GLOBAL Y LIMITS =================
% Precomputes one shared y-axis range per (non-proportion) metric across
% the WHOLE filtered dataset, so every individual genotype/condition line
% plot for that metric uses the same scale and stays visually comparable
% -- proportion metrics always just use a fixed [0,1] axis instead (set
% inline wherever plots are drawn below).

globalY = struct;
for m = 1:length(metricsToPlot)
    metric = metricsToPlot(m);
    vals = T.(metric);
    vals = vals(~isnan(vals));

    if isempty(vals)
        globalY(m).min = 0;
        globalY(m).max = 1;
    else
        globalY(m).min = min(vals);
        globalY(m).max = max(vals);
    end
end

%% ================= STYLE =================

genotypes = sort(unique(T.BaseGenotype));
condListAll = sort(unique(T.Condition));
condListAll = condListAll(condListAll ~= "");

%% ============================================================
%% PART 1: LINE PLOTS — ONE PNG/CSV PER GENOTYPE
%% ============================================================

for m = 1:length(metricsToPlot)

    metric = metricsToPlot(m);

    if contains(metric, "Proportion")
        yMin = 0;
        yMax = 1;
    else
        yMin = globalY(m).min;
        yMax = globalY(m).max;
        yPad = 0.05 * (yMax - yMin);
        if yPad == 0, yPad = 0.02; end   % constant-valued metric: still give it a visible margin
        yMin = max(0, yMin - yPad);
        yMax = yMax + yPad;
    end

    outDir = metricDirs(m);

    for g = 1:length(genotypes)

        genoName = genotypes(g);
        Tg = T(T.BaseGenotype == genoName, :);
        if isempty(Tg), continue; end

        Tgm = Tg(~isnan(Tg.(metric)), :);
        if isempty(Tgm), continue; end

        % Split by Condition if this genotype has any (inducible-vpsT/
        % empty-vector strains); otherwise fall back to splitting by
        % culture age (Hours) for plain strains.
        condObserved = sort(unique(Tgm.Condition));
        condObserved = condObserved(condObserved ~= "");
        useConditionPlot = ~isempty(condObserved);

        fig = figure('Visible','off', 'Position',[100 100 900 600]);
        hold on;

        legHandles = gobjects(0);
        legLabels  = strings(0,1);
        xMaxThisPlot = 0;

        csvSeries = strings(0,1);
        csvX = {};
        csvY = {};
        csvE = {};

        if useConditionPlot

            colorsThis = lines(length(condObserved));

            for c = 1:length(condObserved)
                condName = condObserved(c);
                Tc = Tgm(Tgm.Condition == condName, :);
                if isempty(Tc), continue; end

                % Mean/SEM across whatever replicates share this exact
                % elapsed-minutes value.
                times = sort(unique(Tc.Time_min));
                meanVals = nan(size(times));
                semVals  = nan(size(times));

                for t = 1:length(times)
                    idx = Tc.Time_min == times(t);
                    vals = Tc.(metric)(idx);
                    meanVals(t) = mean(vals, 'omitnan');

                    n = sum(~isnan(vals));
                    if n > 1
                        semVals(t) = std(vals, 'omitnan') / sqrt(n);
                    end
                end

                good = ~isnan(meanVals);
                x = times(good);
                y = meanVals(good);
                e = semVals(good);

                if isempty(x), continue; end
                xMaxThisPlot = max(xMaxThisPlot, max(x));

                % Shaded SEM ribbon only where SEM is defined (>=2
                % replicates); excluded from the legend so each series
                % shows just one line entry, not a line + a ribbon swatch.
                goodErr = ~isnan(e);
                if any(goodErr)
                    xf = x(goodErr);
                    yf = y(goodErr);
                    ef = e(goodErr);

                    fh = fill([xf; flipud(xf)], [yf-ef; flipud(yf+ef)], ...
                        colorsThis(c,:), 'FaceAlpha', 0.2, 'EdgeColor', 'none');
                    fh.Annotation.LegendInformation.IconDisplayStyle = 'off';
                end

                ph = plot(x, y, 'Color', colorsThis(c,:), 'LineWidth', 2);
                legHandles(end+1) = ph;
                legLabels(end+1) = condName;

                csvSeries(end+1,1) = condName;
                csvX{end+1} = x;
                csvY{end+1} = y;
                csvE{end+1} = e;
            end

            plotLabel = "conditions";

        else

            hoursList = sort(unique(Tgm.Hours));
            hourColors = lines(length(hoursList));

            for h = 1:length(hoursList)
                hourVal = hoursList(h);
                Th = Tgm(Tgm.Hours == hourVal, :);
                if isempty(Th), continue; end

                times = sort(unique(Th.Time_min));
                meanVals = nan(size(times));
                semVals  = nan(size(times));

                for t = 1:length(times)
                    idx = Th.Time_min == times(t);
                    vals = Th.(metric)(idx);
                    meanVals(t) = mean(vals, 'omitnan');

                    n = sum(~isnan(vals));
                    if n > 1
                        semVals(t) = std(vals, 'omitnan') / sqrt(n);
                    end
                end

                good = ~isnan(meanVals);
                x = times(good);
                y = meanVals(good);
                e = semVals(good);

                if isempty(x), continue; end
                xMaxThisPlot = max(xMaxThisPlot, max(x));

                goodErr = ~isnan(e);
                if any(goodErr)
                    xf = x(goodErr);
                    yf = y(goodErr);
                    ef = e(goodErr);

                    fh = fill([xf; flipud(xf)], [yf-ef; flipud(yf+ef)], ...
                        hourColors(h,:), 'FaceAlpha', 0.2, 'EdgeColor', 'none');
                    fh.Annotation.LegendInformation.IconDisplayStyle = 'off';
                end

                ph = plot(x, y, 'Color', hourColors(h,:), 'LineWidth', 2);
                legHandles(end+1) = ph;
                legLabels(end+1) = hourVal + " h";

                csvSeries(end+1,1) = hourVal + " h";
                csvX{end+1} = x;
                csvY{end+1} = y;
                csvE{end+1} = e;
            end

            plotLabel = "hours";
        end

        if xMaxThisPlot == 0
            close(fig);
            continue
        end

        xlim([0 xMaxThisPlot]);
        if contains(metric, "Proportion")
            ylim([0 1]);
        else
            ylim([yMin yMax]);
        end

        xlabel("Time (min)");
        ylabel(metric, 'Interpreter', 'none');
        title(genoName + " - " + metric + " (" + plotLabel + ")", 'Interpreter', 'none');

        if ~isempty(legHandles)
            legend(legHandles, legLabels, 'Interpreter', 'none', 'Location', 'eastoutside');
        end

        grid on;
        box on;

        safeName = regexprep(char(genoName), '[^\w\d-]', '_');   % genotype names can contain "∆", spaces, etc.
        outPng = fullfile(outDir, safeName + ".png");
        exportgraphics(fig, outPng, 'Resolution', 300);
        close(fig);

        writeLineCSV(outPng, csvSeries, csvX, csvY, csvE, "Time_min", metric);
    end
end

%% ============================================================
%% PART 2: LINE PLOTS — ONE PNG/CSV PER CONDITION
%% ============================================================
% Transpose of Part 1: one plot per condition, with every genotype that
% has data under that condition as its own colored series. Only
% genotypes with a non-blank Condition ever appear here.

for m = 1:length(metricsToPlot)

    metric = metricsToPlot(m);

    if contains(metric, "Proportion")
        yMin = 0;
        yMax = 1;
    else
        yMin = globalY(m).min;
        yMax = globalY(m).max;
        yPad = 0.05 * (yMax - yMin);
        if yPad == 0, yPad = 0.02; end
        yMin = max(0, yMin - yPad);
        yMax = yMax + yPad;
    end

    outDir = comboDirs(m);

    for c = 1:length(condListAll)

        condName = condListAll(c);
        TcAll = T(T.Condition == condName, :);
        if isempty(TcAll), continue; end

        genosThis = sort(unique(TcAll.BaseGenotype));
        genoColors = lines(length(genosThis));

        fig = figure('Visible','off', 'Position',[100 100 950 650]);
        hold on;

        legHandles = gobjects(0);
        legLabels  = strings(0,1);
        xMaxThisPlot = 0;

        csvSeries = strings(0,1);
        csvX = {};
        csvY = {};
        csvE = {};

        for g = 1:length(genosThis)
            genoName = genosThis(g);
            Tg = TcAll(TcAll.BaseGenotype == genoName, :);
            Tg = Tg(~isnan(Tg.(metric)), :);
            if isempty(Tg), continue; end

            times = sort(unique(Tg.Time_min));
            meanVals = nan(size(times));
            semVals  = nan(size(times));

            for t = 1:length(times)
                idx = Tg.Time_min == times(t);
                vals = Tg.(metric)(idx);
                meanVals(t) = mean(vals, 'omitnan');

                n = sum(~isnan(vals));
                if n > 1
                    semVals(t) = std(vals, 'omitnan') / sqrt(n);
                end
            end

            good = ~isnan(meanVals);
            x = times(good);
            y = meanVals(good);
            e = semVals(good);

            if isempty(x), continue; end
            xMaxThisPlot = max(xMaxThisPlot, max(x));

            goodErr = ~isnan(e);
            if any(goodErr)
                xf = x(goodErr);
                yf = y(goodErr);
                ef = e(goodErr);

                fh = fill([xf; flipud(xf)], [yf-ef; flipud(yf+ef)], ...
                    genoColors(g,:), 'FaceAlpha', 0.2, 'EdgeColor', 'none');
                fh.Annotation.LegendInformation.IconDisplayStyle = 'off';
            end

            ph = plot(x, y, 'Color', genoColors(g,:), 'LineWidth', 2);
            legHandles(end+1) = ph;
            legLabels(end+1) = genoName;

            csvSeries(end+1,1) = genoName;
            csvX{end+1} = x;
            csvY{end+1} = y;
            csvE{end+1} = e;
        end

        if xMaxThisPlot == 0
            close(fig);
            continue
        end

        xlim([0 xMaxThisPlot]);
        if contains(metric, "Proportion")
            ylim([0 1]);
        else
            ylim([yMin yMax]);
        end

        xlabel("Time (min)");
        ylabel(metric, 'Interpreter', 'none');
        title(condName + " - all genotypes - " + metric, 'Interpreter', 'none');

        if ~isempty(legHandles)
            legend(legHandles, legLabels, 'Interpreter', 'none', 'Location', 'eastoutside');
        end

        grid on;
        box on;

        safeCond = regexprep(char(condName), '[^\w\d-]', '_');
        outPng = fullfile(outDir, safeCond + ".png");
        exportgraphics(fig, outPng, 'Resolution', 300);
        close(fig);

        writeLineCSV(outPng, csvSeries, csvX, csvY, csvE, "Time_min", metric);
    end
end

%% ============================================================
%% PART 3: BAR PLOTS
%% ============================================================
% Per-genotype bar comparisons at each fixed target time in barTimes:
% "hours-grouped" bars (+ a companion line plot) for genotypes with no
% condition, split by culture age; and per-condition bars for genotypes
% that do have a condition.

allGenotypes = sort(unique(T.BaseGenotype));

for bm = 1:length(barMetrics)

    barMetric = barMetrics(bm);
    genoBarDir = barGenoDirs(char(barMetric));

    if ~ismember(barMetric, string(T.Properties.VariableNames))
        warning("Column %s not found in table, skipping.", barMetric);
        continue
    end

    hourGroupedDir = fullfile(genoBarDir, "HoursGrouped");
    if ~exist(hourGroupedDir, 'dir'), mkdir(hourGroupedDir); end

    for tt = 1:length(barTimes)
        d1 = fullfile(hourGroupedDir, sprintf('%dmin', barTimes(tt)));
        if ~exist(d1, 'dir'), mkdir(d1); end
    end

    %% ---------------- HOURS GROUPED BAR + LINE PLOTS ----------------

    genoHoursOnly = strings(0,1);
    hourListAll = [];

    for g = 1:length(allGenotypes)
        genoName = allGenotypes(g);
        Tg = T(T.BaseGenotype == genoName & ~isnan(T.(barMetric)), :);
        if isempty(Tg), continue; end

        condObserved = sort(unique(Tg.Condition));
        condObserved = condObserved(condObserved ~= "");
        useConditionPlot = ~isempty(condObserved);

        if ~useConditionPlot
            genoHoursOnly(end+1) = genoName; %#ok<AGROW>
            hourListAll = [hourListAll; unique(Tg.Hours)]; %#ok<AGROW>
        end
    end

    genoHoursOnly = sort(unique(genoHoursOnly));
    hourListAll = sort(unique(hourListAll));

    for tt = 1:length(barTimes)

        targetTime = barTimes(tt);

        if ~isempty(genoHoursOnly) && ~isempty(hourListAll)

            meanMat = nan(length(hourListAll), length(genoHoursOnly));
            errMat  = nan(length(hourListAll), length(genoHoursOnly));
            repCell = cell(length(hourListAll), length(genoHoursOnly));

            for h = 1:length(hourListAll)
                hourVal = hourListAll(h);

                for g = 1:length(genoHoursOnly)
                    genoName = genoHoursOnly(g);

                    Tsel = T(T.BaseGenotype == genoName & ...
                             T.Hours == hourVal & ...
                             ~isnan(T.(barMetric)), :);

                    [valsRep, errVal, ~] = getBarReplicateValuesAndError( ...
                        Tsel, barMetric, targetTime, barReplicateGroups);

                    repCell{h,g} = valsRep;
                    meanMat(h,g) = mean(valsRep, 'omitnan');
                    errMat(h,g)  = errVal;
                end
            end

            if ~all(all(isnan(meanMat)))

                fig = figure('Visible','off', 'Position',[100 100 1200 650]);
                bh = bar(meanMat, 'grouped');
                hold on;

                genoColors = lines(length(genoHoursOnly));

                for g = 1:length(genoHoursOnly)
                    bh(g).FaceColor = genoColors(g,:);
                end

                for g = 1:length(genoHoursOnly)
                    xg = bh(g).XEndPoints;
                    errorbar(xg, meanMat(:,g), errMat(:,g), ...
                        'k', 'LineStyle', 'none', 'LineWidth', 1);
                end

                % Overlay each individual replicate value as a small
                % jittered dot (darker shade of the bar's color) on top
                % of its bar, so the spread behind each mean/SEM is
                % visible directly rather than only implied by the
                % error bar.
                for g = 1:length(genoHoursOnly)
                    xg = bh(g).XEndPoints;
                    dotColor = max(0, genoColors(g,:) * 0.55);

                    for h = 1:length(hourListAll)
                        valsRep = repCell{h,g};
                        valsRep = valsRep(~isnan(valsRep));

                        if isempty(valsRep), continue; end

                        jitter = linspace(-0.015, 0.015, max(length(valsRep),2));
                        jitter = jitter(1:length(valsRep));

                        scatter(xg(h) + jitter, valsRep, 28, ...
                            'MarkerFaceColor', dotColor, ...
                            'MarkerEdgeColor', dotColor, ...
                            'MarkerFaceAlpha', 0.95, ...
                            'MarkerEdgeAlpha', 0.95);
                    end
                end

                set(gca, 'XTick', 1:length(hourListAll), ...
                    'XTickLabel', string(hourListAll) + "h");

                if contains(barMetric, "Proportion")
                    ylim([0 1]);
                end

                xlabel("Culture time");
                ylabel(barMetric, 'Interpreter', 'none');
                title("Hours-only genotypes - " + barMetric + " - " + targetTime + " min", ...
                    'Interpreter', 'none');
                legend(genoHoursOnly, 'Interpreter', 'none', 'Location', 'eastoutside');

                grid on;
                box on;

                outFileBar = fullfile(hourGroupedDir, sprintf('%dmin', targetTime), ...
                    "HoursGrouped_AllGenotypes.png");
                exportgraphics(fig, outFileBar, 'Resolution', 300);
                close(fig);

                writeBarCSV(outFileBar, string(hourListAll) + "h", genoHoursOnly, ...
                    meanMat, errMat, repCell, "Hours", barMetric);

                % Same hours-grouped data, redrawn as a connected line
                % plot (one line per genotype across hours) instead of
                % grouped bars, for an alternative view of the same
                % values -- saved alongside the bar chart above.
                fig = figure('Visible','off', 'Position',[100 100 1200 650]);
                hold on;

                legHandles = gobjects(0);
                legLabels = strings(0,1);

                lineSeries = strings(0,1);
                lineX = {};
                lineY = {};
                lineE = {};

                for g = 1:length(genoHoursOnly)

                    genoName = genoHoursOnly(g);

                    y = meanMat(:,g);
                    e = errMat(:,g);
                    x = hourListAll;

                    good = ~isnan(y);
                    xg = x(good);
                    yg = y(good);
                    eg = e(good);

                    if isempty(xg), continue; end

                    goodErr = ~isnan(eg);
                    if any(goodErr)
                        xf = xg(goodErr);
                        yf = yg(goodErr);
                        ef = eg(goodErr);

                        fh = fill([xf; flipud(xf)], ...
                                  [yf-ef; flipud(yf+ef)], ...
                                  genoColors(g,:), ...
                                  'FaceAlpha', 0.2, ...
                                  'EdgeColor', 'none');
                        fh.Annotation.LegendInformation.IconDisplayStyle = 'off';
                    end

                    ph = plot(xg, yg, 'Color', genoColors(g,:), 'LineWidth', 2);

                    legHandles(end+1) = ph;
                    legLabels(end+1) = genoName;

                    lineSeries(end+1,1) = genoName;
                    lineX{end+1} = xg;
                    lineY{end+1} = yg;
                    lineE{end+1} = eg;
                end

                if contains(barMetric, "Proportion")
                    ylim([0 1]);
                end

                xlim([min(hourListAll) max(hourListAll)]);
                xlabel("Culture time");
                ylabel(barMetric, 'Interpreter', 'none');
                title("Hours-only genotypes - " + barMetric + " - " + targetTime + " min - line plot", ...
                    'Interpreter', 'none');
                legend(legHandles, legLabels, 'Interpreter', 'none', 'Location', 'eastoutside');

                grid on;
                box on;

                outFileLine = fullfile(hourGroupedDir, sprintf('%dmin', targetTime), ...
                    "HoursGrouped_AllGenotypes_Line.png");
                exportgraphics(fig, outFileLine, 'Resolution', 300);
                close(fig);

                writeLineCSV(outFileLine, lineSeries, lineX, lineY, lineE, "Hours", barMetric);
            end
        end
    end

    %% ---------------- CONDITION PER-GENOTYPE BAR PLOTS ----------------

    for tt = 1:length(barTimes)

        targetTime = barTimes(tt);

        for g = 1:length(allGenotypes)

            genoName = allGenotypes(g);
            Tg = T(T.BaseGenotype == genoName & ~isnan(T.(barMetric)), :);
            if isempty(Tg), continue; end

            condObserved = sort(unique(Tg.Condition));
            condObserved = condObserved(condObserved ~= "");
            useConditionPlot = ~isempty(condObserved);

            if ~useConditionPlot
                continue
            end

            xCats_bar   = condObserved;
            xLabels_bar = xCats_bar;
            colors_bar  = lines(length(xCats_bar));
            xLabelText  = "Condition";

            meanVals_bar = nan(1, length(xCats_bar));
            errVals_bar  = nan(1, length(xCats_bar));
            repVals_bar  = cell(1, length(xCats_bar));

            for ci = 1:length(xCats_bar)
                Tsel = Tg(Tg.Condition == xCats_bar(ci), :);

                [valsRep, errVal, ~] = getBarReplicateValuesAndError( ...
                    Tsel, barMetric, targetTime, barReplicateGroups);

                repVals_bar{ci}  = valsRep;
                meanVals_bar(ci) = mean(valsRep, 'omitnan');
                errVals_bar(ci)  = errVal;
            end

            if all(isnan(meanVals_bar)), continue; end

            titleText = genoName + " - " + barMetric + " - conditions - " + targetTime + " min";

            safeName   = regexprep(char(genoName), '[^\w\d-]', '_');
            outFileBar = fullfile(genoBarDir, sprintf('%dmin', targetTime), safeName + ".png");

            saveBarFigure(meanVals_bar, errVals_bar, repVals_bar, colors_bar, ...
                xLabels_bar, xLabelText, barMetric, titleText, outFileBar, false);

            writeBarCSV(outFileBar, xLabels_bar, "", ...
                meanVals_bar, errVals_bar, repVals_bar, xLabelText, barMetric);
        end
    end
end

%% ============================================================
%% PART 4: BAR PLOTS — PER CONDITION
%% ============================================================
% Transpose of Part 3's condition bars: one bar plot per condition,
% comparing every genotype that has data under that condition.

for bm = 1:length(barMetrics)

    barMetric = barMetrics(bm);
    condBarDir = barCondDirs(char(barMetric));

    if ~ismember(barMetric, string(T.Properties.VariableNames))
        warning("Column %s not found in table, skipping.", barMetric);
        continue
    end

    for tt = 1:length(barTimes)

        targetTime = barTimes(tt);

        for ci = 1:length(condListAll)

            condName = condListAll(ci);
            TcAll = T(T.Condition == condName & ~isnan(T.(barMetric)), :);
            if isempty(TcAll), continue; end

            genosThis = sort(unique(TcAll.BaseGenotype));
            colors_bar = lines(length(genosThis));

            meanVals_bar = nan(1, length(genosThis));
            errVals_bar  = nan(1, length(genosThis));
            repVals_bar  = cell(1, length(genosThis));

            for g = 1:length(genosThis)
                genoName = genosThis(g);
                Tsel = TcAll(TcAll.BaseGenotype == genoName, :);

                [valsRep, errVal, ~] = getBarReplicateValuesAndError( ...
                    Tsel, barMetric, targetTime, barReplicateGroups);

                repVals_bar{g}  = valsRep;
                meanVals_bar(g) = mean(valsRep, 'omitnan');
                errVals_bar(g)  = errVal;
            end

            if all(isnan(meanVals_bar)), continue; end

            safeCond  = regexprep(char(condName), '[^\w\d-]', '_');
            titleText = condName + " - " + barMetric + " - " + targetTime + " min";
            outFileBar = fullfile(condBarDir, sprintf('%dmin', targetTime), safeCond + ".png");

            saveBarFigure(meanVals_bar, errVals_bar, repVals_bar, colors_bar, ...
                genosThis, "Genotype", barMetric, titleText, outFileBar, true);

            writeBarCSV(outFileBar, genosThis, "", ...
                meanVals_bar, errVals_bar, repVals_bar, "Genotype", barMetric);
        end
    end
end

disp('Done: all PNG plots and matching CSV files saved.');

%% ================= HELPER FUNCTIONS =================

function [valsRep, errVal, usedReplicates] = getBarReplicateValuesAndError(Tsub, metricName, targetTime, replicateGroups)
% For one already-filtered (genotype x hour, or genotype x condition)
% subset, finds each candidate replicate's value nearest to targetTime
% (not an exact match -- real elapsed times don't always land on round
% numbers), tries every group in replicateGroups (e.g. reps [1 2 3] vs.
% [4 5 6]) and keeps whichever group has the most actual (non-NaN)
% replicate values, then returns those per-replicate values plus one
% combined error estimate.
%
% Error estimate depends on the metric: MeanFractionAreaFilledWithSpots
% and MeanSpotsPerObject are themselves already per-image MEANS across
% objects, so their spread has to account for that nested within-image
% variability -- each kept replicate's own within-image
% Std*/StdSpotsPerObject value is pooled (sqrt(sum(sigma^2))/n) rather
% than computing the across-replicate SEM of the replicate-level means.
% Every other metric uses the plain across-replicate SEM instead.

if isempty(Tsub)
    valsRep = nan(1,3);
    errVal = NaN;
    usedReplicates = [];
    return
end

bestVals = [];
bestSigma = [];
bestReps = [];
bestN = -1;

for gg = 1:length(replicateGroups)

    replicateList = replicateGroups{gg};
    valsThis = nan(1, length(replicateList));
    sigThis  = nan(1, length(replicateList));

    for r = 1:length(replicateList)

        repVal = replicateList(r);
        Tr = Tsub(Tsub.Replicate == repVal, :);

        if isempty(Tr), continue; end

        d = abs(Tr.Time_min - targetTime);
        minD = min(d);
        Tnear = Tr(d == minD, :);

        vals = Tnear.(metricName);
        vals = vals(~isnan(vals));

        if isempty(vals), continue; end

        valsThis(r) = mean(vals, 'omitnan');

        switch metricName
            case "MeanFractionAreaFilledWithSpots"
                if ismember("StdFractionAreaFilledWithSpots", string(Tnear.Properties.VariableNames))
                    s = Tnear.StdFractionAreaFilledWithSpots;
                    s = s(~isnan(s));
                    if ~isempty(s)
                        sigThis(r) = mean(s, 'omitnan');
                    end
                end

            case "MeanSpotsPerObject"
                if ismember("StdSpotsPerObject", string(Tnear.Properties.VariableNames))
                    s = Tnear.StdSpotsPerObject;
                    s = s(~isnan(s));
                    if ~isempty(s)
                        sigThis(r) = mean(s, 'omitnan');
                    end
                end
        end
    end

    nGood = sum(~isnan(valsThis));

    if nGood > bestN
        bestN = nGood;
        bestVals = valsThis;
        bestSigma = sigThis;
        bestReps = replicateList;
    end
end

if isempty(bestVals) || bestN <= 0
    valsRep = nan(1,3);
    errVal = NaN;
    usedReplicates = [];
    return
end

valsRep = bestVals;
usedReplicates = bestReps;

good = ~isnan(valsRep);
n = sum(good);

if n == 0
    errVal = NaN;
    return
end

switch metricName
    case {"MeanFractionAreaFilledWithSpots", "MeanSpotsPerObject"}
        sig = bestSigma(good);
        sig = sig(~isnan(sig));

        if isempty(sig)
            errVal = NaN;
        else
            errVal = sqrt(sum(sig.^2)) / n;
        end

    otherwise
        if n > 1
            errVal = std(valsRep(good), 'omitnan') / sqrt(n);
        else
            errVal = NaN;
        end
end

end

function saveBarFigure(meanVals, errVals, repVals, colorsThis, xLabels, xLabelText, yLabelText, titleText, outFile, doXTickAngle)
% Generic single-row-of-bars plot: one colored bar per category, a black
% error bar, and each individual replicate value overlaid as a small
% jittered dot in a darkened shade of that bar's color.

fig = figure('Visible','off', 'Position',[100 100 1000 650]);
bh = bar(meanVals, 'FaceColor', 'flat');
hold on;

for i = 1:length(meanVals)
    bh.CData(i,:) = colorsThis(i,:);
end

errorbar(1:length(meanVals), meanVals, errVals, ...
    'k', 'LineStyle', 'none', 'LineWidth', 1);

for i = 1:length(meanVals)
    vR = repVals{i};
    vR = vR(~isnan(vR));
    if isempty(vR), continue; end

    dotColor = max(0, colorsThis(i,:) * 0.55);
    jitter = linspace(-0.015, 0.015, max(length(vR),2));
    jitter = jitter(1:length(vR));

    scatter(i + jitter, vR, 28, ...
        'MarkerFaceColor', dotColor, ...
        'MarkerEdgeColor', dotColor, ...
        'MarkerFaceAlpha', 0.95, ...
        'MarkerEdgeAlpha', 0.95);
end

set(gca, 'XTick', 1:length(meanVals), 'XTickLabel', xLabels);
if doXTickAngle
    xtickangle(30);
end

if contains(yLabelText, "Proportion")
    ylim([0 1]);
end

xlabel(xLabelText);
ylabel(yLabelText, 'Interpreter', 'none');
title(titleText, 'Interpreter', 'none');

grid on;
box on;

exportgraphics(fig, outFile, 'Resolution', 300);
close(fig);

end

function writeLineCSV(outPng, seriesCell, xCell, yCell, errCell, xName, yName)
% Flattens the per-series cell-array-of-vectors used for line plots into
% one long-format table (Series, x, y, Error -- one row per series per
% x-value) and writes it next to the PNG with a matching filename.

PlotSeries = strings(0,1);
X = [];
Y = [];
Error = [];

for s = 1:length(seriesCell)
    x = xCell{s};
    y = yCell{s};
    e = errCell{s};

    for i = 1:length(x)
        PlotSeries(end+1,1) = seriesCell(s);
        X(end+1,1) = x(i);
        Y(end+1,1) = y(i);

        if isempty(e) || length(e) < i
            Error(end+1,1) = NaN;
        else
            Error(end+1,1) = e(i);
        end
    end
end

Tout = table(PlotSeries, X, Y, Error);
Tout.Properties.VariableNames = ["Series", xName, yName, "Error"];

outCsv = replace(outPng, ".png", ".csv");
writetable(Tout, outCsv);

end

function writeBarCSV(outPng, groupLabels, barLabels, meanVals, errVals, repVals, xName, yName)
% Same flattening idea as writeLineCSV, but for bar-chart data -- handles
% both a single row of bars (meanVals a vector: one group category, no
% further sub-grouping) and a full grouped-bar matrix (meanVals a matrix:
% rows and columns are two different category axes, as in the
% HoursGrouped bar chart). Writes one row per (group, bar, replicate)
% combination, or one row with ReplicateIndex/ReplicateValue left NaN
% if that group/bar cell had no valid replicate values at all.

Group = strings(0,1);
Bar = strings(0,1);
Mean = [];
Error = [];
ReplicateIndex = [];
ReplicateValue = [];

if isvector(meanVals)

    for i = 1:length(meanVals)
        vals = repVals{i};
        vals = vals(~isnan(vals));

        if isempty(vals)
            Group(end+1,1) = groupLabels(i);
            Bar(end+1,1) = "";
            Mean(end+1,1) = meanVals(i);
            Error(end+1,1) = errVals(i);
            ReplicateIndex(end+1,1) = NaN;
            ReplicateValue(end+1,1) = NaN;
        else
            for r = 1:length(vals)
                Group(end+1,1) = groupLabels(i);
                Bar(end+1,1) = "";
                Mean(end+1,1) = meanVals(i);
                Error(end+1,1) = errVals(i);
                ReplicateIndex(end+1,1) = r;
                ReplicateValue(end+1,1) = vals(r);
            end
        end
    end

else

    for i = 1:size(meanVals,1)
        for j = 1:size(meanVals,2)
            vals = repVals{i,j};
            vals = vals(~isnan(vals));

            if isempty(vals)
                Group(end+1,1) = groupLabels(i);
                Bar(end+1,1) = barLabels(j);
                Mean(end+1,1) = meanVals(i,j);
                Error(end+1,1) = errVals(i,j);
                ReplicateIndex(end+1,1) = NaN;
                ReplicateValue(end+1,1) = NaN;
            else
                for r = 1:length(vals)
                    Group(end+1,1) = groupLabels(i);
                    Bar(end+1,1) = barLabels(j);
                    Mean(end+1,1) = meanVals(i,j);
                    Error(end+1,1) = errVals(i,j);
                    ReplicateIndex(end+1,1) = r;
                    ReplicateValue(end+1,1) = vals(r);
                end
            end
        end
    end
end

Tout = table(Group, Bar, Mean, Error, ReplicateIndex, ReplicateValue);
Tout.Properties.VariableNames = [xName, "Bar", yName, "Error", "ReplicateIndex", "ReplicateValue"];

outCsv = replace(outPng, ".png", ".csv");
writetable(Tout, outCsv);

end
