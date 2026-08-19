%% CompileCultureTracking_SurfaceMiddleSummaries.m
% Author: Hasreet Gill, Bassler Lab, HHMI/Princeton University
% Revised with Claude (Sonnet 5, Anthropic)

% Post-processing script that runs AFTER SegmentCultureTimelapse.m has segmented a
% dataset (Binary/Binary_Watershed/Surface_Binary/Surface_Binary_Watershed
% folders per date) AND after TrackMate (Fiji) has been run on those
% binary stacks to produce "<name>-tracks.csv" / "<name>-spots.csv" pairs
% inside those same folders. For every culture type / date folder under
% parentPath, this:
%   1. Compiles TrackMate's surface-projection tracking output into one
%      per-position (M) CSV (SURFACE TRACKING COMPILE).
%   2. Compiles TrackMate's full-Z-stack tracking output the same way,
%      combining every timepoint file for a given position into one
%      table (FULL STACK TRACKING COMPILE) -- see note below on what
%      TrackMate's "frame" axis means in each case.
%   3. Maps each fine-threshold 3D object (from culture_della.m's own
%      M*_3D*.csv) onto its nearest 2D tracked spot at the same
%      timepoint/Z-slice, attaching that spot's track ID and stats
%      (MAP 2D TRACKS ONTO 3D OBJECTS).
%   4. Estimates which Z-slice is "the surface" at every timepoint, by
%      four independent methods, and writes a per-position comparison
%      table (SURFACE SLICE SUMMARY PER M).
%   5. Smooths the track-based surface-slice estimate over time and uses
%      it (plus a simple "half the surface index" middle estimate) to
%      re-slice every 2D/3D/tracking/coarse/watershed table down to just
%      the surface slice or a +/-5-slice middle band, writing a parallel
%      "Surface_FromStacks_*" / "Middle_FromStacks_*" CSV for each
%      (FROM STACKS SURFACE + MIDDLE CSVs).
%
% TrackMate's own "frame" column means two different things here: for the
% surface-projection files it's the real imaging timepoint (one 2D slice
% per timepoint was tracked over time); for the full-stack files it's
% really the Z-slice WITHIN one timepoint's volume (TrackMate tracked
% each single-timepoint stack's Z-slices as if they were a movie, since
% object identity across the Z-stack of one instant is what's needed
% there, not tracking across real time). Sections 1 and 2 below convert
% each case's TrackMate output into a common column layout accordingly.

clear
clc

parentPath = uigetdir();

% Which top-level dataset folder(s) under parentPath to process. Toggle
% by commenting/uncommenting -- currently limited to "Protein addition".
% cultureTypes = ["Monoculture","Coculture"];
cultureTypes = ["Protein addition"];

for c = 1:length(cultureTypes)

    cultureName = cultureTypes(c);
    culturePath = fullfile(parentPath,cultureName);

    % Date-stamped experiment folders under this culture type, excluding
    % "SurfacePlots" (an output folder from a separate plotting step, not
    % an experiment date) alongside "." and "..".
    dateFolders = dir(culturePath);
    dateFolders = dateFolders([dateFolders.isdir]);
    dateFolders = dateFolders(~ismember({dateFolders.name},{'.','..','SurfacePlots'}));

    for d = 1:length(dateFolders)

        dateName = string(dateFolders(d).name);
        datePath = fullfile(culturePath,dateName);
        sumFolder = fullfile(datePath,"Summary");

        if ~isfolder(sumFolder)
            mkdir(sumFolder)
        end

        %% ===============================
        %% SURFACE TRACKING COMPILE
        %% ===============================
        % TrackMate was run per position/channel on the surface-projection
        % binary (and its watershed variant), one -tracks/-spots CSV pair
        % per position. Here TrackMate's frame axis is real time.

        folderList = { ...
            "Surface_Binary", "Tracking"; ...
            "Surface_Binary_Watershed", "Tracking_Watershed"};

        for k = 1:size(folderList,1)

            srcFolder = fullfile(datePath,folderList{k,1});
            outLabel = folderList{k,2};

            if ~isfolder(srcFolder)
                continue
            end

            trackFiles = dir(fullfile(srcFolder,"*-tracks.csv"));

            for f = 1:length(trackFiles)

                trackName = trackFiles(f).name;
                spotName = erase(trackName,"-tracks.csv") + "-spots.csv";

                trackPath = fullfile(srcFolder,trackName);
                spotPath = fullfile(srcFolder,spotName);

                if ~isfile(spotPath)
                    continue
                end

                tok = regexp(trackName,"Surf_M(\d+)",'tokens','once');
                M = str2double(tok{1});

                if contains(trackName,"488")
                    outName = "M" + M + "_Surface_" + outLabel + "_488.csv";
                elseif contains(trackName,"561")
                    outName = "M" + M + "_Surface_" + outLabel + "_561.csv";
                else
                    outName = "M" + M + "_Surface_" + outLabel + ".csv";
                end

                outPath = fullfile(sumFolder,outName);

                % TrackMate CSV exports have 3 extra header/units rows
                % below the real column-name row, so real data starts at
                % row 5 -- DataLines skips past them.
                optsT = detectImportOptions(trackPath);
                optsT.VariableNamesLine = 1;
                optsT.DataLines = [5 Inf];
                Ttracks = readtable(trackPath,optsT);

                optsS = detectImportOptions(spotPath);
                optsS.VariableNamesLine = 1;
                optsS.DataLines = [5 Inf];
                Tspots = readtable(spotPath,optsS);

                % Spots with no TRACK_ID were detected but never linked
                % into a track by TrackMate -- drop them, only tracked
                % spots are useful here.
                keepSpots = ~isnan(Tspots.TRACK_ID);
                Tkeep = Tspots(keepSpots,:);

                % Attach each spot's parent track's summary stats.
                TtrackKeep = Ttracks(:,["TRACK_ID","TRACK_START","TRACK_STOP","TRACK_DISPLACEMENT","TRACK_DURATION","TOTAL_DISTANCE_TRAVELED"]);
                Tkeep = innerjoin(Tkeep,TtrackKeep,"Keys","TRACK_ID");

                % Rename/reformat into this pipeline's column convention;
                % TrackMate's FRAME/TRACK_START/TRACK_STOP are 0-based, so
                % +1 brings them in line with this pipeline's 1-based
                % frame numbering.
                Tout = table();
                Tout.T = Tkeep.FRAME + 1;
                Tout.ObjectID = Tkeep.ID;
                Tout.Area = Tkeep.AREA;
                Tout.X = Tkeep.POSITION_X;
                Tout.Y = Tkeep.POSITION_Y;
                Tout.TRACK_ID = Tkeep.TRACK_ID;
                Tout.TRACK_START = Tkeep.TRACK_START + 1;
                Tout.TRACK_STOP = Tkeep.TRACK_STOP + 1;
                Tout.TrackDisplacement = Tkeep.TRACK_DISPLACEMENT;
                Tout.TrackDuration = Tkeep.TRACK_DURATION;
                Tout.TotalDistanceTraveled = Tkeep.TOTAL_DISTANCE_TRAVELED;

                Tout = sortrows(Tout,"T");

                writetable(Tout,outPath)

            end
        end

        %% ===============================
        %% FULL STACK TRACKING COMPILE
        %% ===============================
        % Same TrackMate compile as above, but for the full Z-stack
        % binaries -- there, TrackMate tracked Z-slices WITHIN one
        % timepoint's volume, so its "frame" is really Z, and the actual
        % timepoint (Tval) instead comes from the filename. Every
        % timepoint file for the same position/channel is accumulated
        % into one combined table per position before being written out.

        folderList = { ...
            "Binary", "Tracking"; ...
            "Binary_Watershed", "Tracking_Watershed"};

        for k = 1:size(folderList,1)

            srcFolder = fullfile(datePath,folderList{k,1});
            outLabel = folderList{k,2};

            if ~isfolder(srcFolder)
                continue
            end

            trackFiles = dir(fullfile(srcFolder,"*-tracks.csv"));

            allTout = struct();

            for f = 1:length(trackFiles)

                trackName = trackFiles(f).name;
                spotName = erase(trackName,"-tracks.csv") + "-spots.csv";

                trackPath = fullfile(srcFolder,trackName);
                spotPath = fullfile(srcFolder,spotName);

                if ~isfile(spotPath)
                    continue
                end

                tok = regexp(trackName,"Surf_M(\d+)_(\d+)_",'tokens','once');
                M = str2double(tok{1});
                Tval = str2double(tok{2});

                if contains(trackName,"488")
                    outKey = "M" + M + "_488";
                elseif contains(trackName,"561")
                    outKey = "M" + M + "_561";
                else
                    outKey = "M" + M;
                end

                optsT = detectImportOptions(trackPath);
                optsT.VariableNamesLine = 1;
                optsT.DataLines = [5 Inf];
                Ttracks = readtable(trackPath,optsT);

                optsS = detectImportOptions(spotPath);
                optsS.VariableNamesLine = 1;
                optsS.DataLines = [5 Inf];
                Tspots = readtable(spotPath,optsS);

                % Some timepoints may have had no detected/tracked
                % objects at all -- skip rather than error.
                if isempty(Tspots)
                    continue
                end

                keepSpots = ~isnan(Tspots.TRACK_ID);
                Tkeep = Tspots(keepSpots,:);

                if isempty(Tkeep)
                    continue
                end

                TtrackKeep = Ttracks(:,["TRACK_ID","TRACK_START","TRACK_STOP","TRACK_DISPLACEMENT","TRACK_DURATION","TOTAL_DISTANCE_TRAVELED"]);
                Tkeep = innerjoin(Tkeep,TtrackKeep,"Keys","TRACK_ID");

                % T is the real timepoint (constant for this whole file);
                % Z repurposes TrackMate's POSITION_T (its native time
                % field), which here actually holds the Z-slice index.
                Tout = table();
                Tout.T = repmat(Tval,height(Tkeep),1);
                Tout.ObjectID = Tkeep.ID;
                Tout.Area = Tkeep.AREA;
                Tout.X = Tkeep.POSITION_X;
                Tout.Y = Tkeep.POSITION_Y;
                Tout.Z = Tkeep.POSITION_T + 1;
                Tout.TRACK_ID = Tkeep.TRACK_ID;
                Tout.TRACK_START = Tkeep.TRACK_START + 1;
                Tout.TRACK_STOP = Tkeep.TRACK_STOP + 1;
                Tout.TrackDisplacement = Tkeep.TRACK_DISPLACEMENT;
                Tout.TrackDuration = Tkeep.TRACK_DURATION;
                Tout.TotalDistanceTraveled = Tkeep.TOTAL_DISTANCE_TRAVELED;

                if ~isfield(allTout,outKey)
                    allTout.(outKey) = Tout;
                else
                    allTout.(outKey) = [allTout.(outKey); Tout];
                end

            end

            outKeys = fieldnames(allTout);

            for q = 1:length(outKeys)

                thisKey = outKeys{q};
                Tout = allTout.(thisKey);

                Tout = sortrows(Tout,["T","Z"]);

                if endsWith(thisKey,"_488")
                    Mstr = erase(thisKey,"_488");
                    outName = Mstr + "_2D_" + outLabel + "_488.csv";
                elseif endsWith(thisKey,"_561")
                    Mstr = erase(thisKey,"_561");
                    outName = Mstr + "_2D_" + outLabel + "_561.csv";
                else
                    Mstr = thisKey;
                    outName = Mstr + "_2D_" + outLabel + ".csv";
                end

                outPath = fullfile(sumFolder,outName);
                writetable(Tout,outPath)

            end
        end

        %% ===============================
        %% MAP 2D TRACKS ONTO 3D OBJECTS
        %% ===============================
        % For every fine-threshold full-3D summary table culture_della.m
        % already wrote (M*_3D*.csv, excluding Surface/Middle/Tracking/
        % Coarse variants -- those either aren't 3D-object tables or have
        % no matching 2D-tracking counterpart here), attach each 3D
        % object's nearest same-(timepoint,Z) 2D tracked spot and its
        % track ID/stats.

        threeDFiles = dir(fullfile(sumFolder,"M*_3D*.csv"));

        for f = 1:length(threeDFiles)

            file3D = threeDFiles(f).name;

            if contains(file3D,"Surface") || contains(file3D,"Middle") || contains(file3D,"Tracking")
                continue
            end

            if contains(file3D,"Coarse")
                continue
            end

            tok = regexp(file3D,"M(\d+)_",'tokens','once');
            M = str2double(tok{1});

            if contains(file3D,"_488")
                suffix = "_488";
            elseif contains(file3D,"_561")
                suffix = "_561";
            else
                suffix = "";
            end

            if contains(file3D,"Watershed")
                file2DTr = "M" + M + "_2D_Tracking_Watershed" + suffix + ".csv";
                outName  = "M" + M + "_3D_Tracking_Watershed" + suffix + ".csv";
            else
                file2DTr = "M" + M + "_2D_Tracking" + suffix + ".csv";
                outName  = "M" + M + "_3D_Tracking" + suffix + ".csv";
            end

            path3D   = fullfile(sumFolder,file3D);
            path2DTr = fullfile(sumFolder,file2DTr);
            outPath  = fullfile(sumFolder,outName);

            if ~isfile(path3D) || ~isfile(path2DTr)
                continue
            end

            T3  = readtable(path3D);
            T2t = readtable(path2DTr);

            Tout = T3;
            Tout.TrackingObjectID = nan(height(T3),1);
            Tout.TRACK_ID = nan(height(T3),1);
            Tout.TRACK_START = nan(height(T3),1);
            Tout.TRACK_STOP = nan(height(T3),1);
            Tout.TrackDisplacement = nan(height(T3),1);
            Tout.TrackDuration = nan(height(T3),1);
            Tout.TotalDistanceTraveled = nan(height(T3),1);
            Tout.TrackingDistanceXY = nan(height(T3),1);

            if isempty(T3) || isempty(T2t)
                writetable(Tout,outPath)
                continue
            end

            allT = unique(T3.T);

            for tt = 1:length(allT)

                thisT = allT(tt);

                idx3T = T3.T == thisT;
                idx2T = T2t.T == thisT;

                if ~any(idx3T) || ~any(idx2T)
                    continue
                end

                zVals = unique(round(T3.Z(idx3T)));

                for z = 1:length(zVals)

                    thisZ = zVals(z);

                    objIdx = find(idx3T & round(T3.Z) == thisZ);
                    spotIdx = find(idx2T & T2t.Z == thisZ);

                    if isempty(objIdx) || isempty(spotIdx)
                        continue
                    end

                    objX = T3.X(objIdx);
                    objY = T3.Y(objIdx);

                    spotX = T2t.X(spotIdx);
                    spotY = T2t.Y(spotIdx);

                    % Nearest 2D tracked spot per 3D object (greedy,
                    % not a unique one-to-one assignment -- more than
                    % one 3D object can claim the same spot).
                    D = (objX - spotX').^2 + (objY - spotY').^2;
                    [minD, minI] = min(D,[],2);

                    matchIdx = spotIdx(minI);

                    Tout.TrackingObjectID(objIdx) = T2t.ObjectID(matchIdx);
                    Tout.TRACK_ID(objIdx) = T2t.TRACK_ID(matchIdx);
                    Tout.TRACK_START(objIdx) = T2t.TRACK_START(matchIdx);
                    Tout.TRACK_STOP(objIdx) = T2t.TRACK_STOP(matchIdx);
                    Tout.TrackDisplacement(objIdx) = T2t.TrackDisplacement(matchIdx);
                    Tout.TrackDuration(objIdx) = T2t.TrackDuration(matchIdx);
                    Tout.TotalDistanceTraveled(objIdx) = T2t.TotalDistanceTraveled(matchIdx);
                    Tout.TrackingDistanceXY(objIdx) = sqrt(minD);  % match-quality check

                end
            end

            writetable(Tout,outPath)

        end
    end
end

%% ===============================
%% SURFACE SLICE SUMMARY PER M
%% ===============================
% Independent second pass over every culture/date: for each position,
% estimates which Z-slice is "the surface" at every timepoint by four
% different signals, so downstream steps can pick a consistent, smoothed
% surface Z instead of relying on the fixed surfIdx used at acquisition
% time in culture_della.m.

surfaceSliceSummary = struct();

for c = 1:length(cultureTypes)

    cultureName = cultureTypes(c);
    culturePath = fullfile(parentPath,cultureName);

    dateFolders = dir(culturePath);
    dateFolders = dateFolders([dateFolders.isdir]);
    dateFolders = dateFolders(~ismember({dateFolders.name},{'.','..','SurfacePlots'}));

    for d = 1:length(dateFolders)

        dateName = string(dateFolders(d).name);
        datePath = fullfile(culturePath,dateName);
        sumFolder = fullfile(datePath,"Summary");

        if ~isfolder(sumFolder)
            continue
        end

        allFiles = dir(fullfile(sumFolder,"M*_2D*.csv"));

        keepIdx = ~contains({allFiles.name},"Surface");
        base2DFiles = allFiles(keepIdx);

        for f = 1:length(base2DFiles)

            fname = base2DFiles(f).name;

            % Only the plain fine-threshold per-slice 2D table (from
            % culture_della.m) drives which M/channel combos get a
            % surface-slice estimate here.
            if contains(fname,"Tracking") || contains(fname,"Watershed") || contains(fname,"Coarse")
                continue
            end

            if contains(fname,"_488")
                tok = regexp(fname,"M(\d+)_2D_488\.csv",'tokens','once');
                M = str2double(tok{1});
                suffix = "_488";
                file2D = "M" + M + "_2D_488.csv";
                file3D = "M" + M + "_3D_488.csv";
                fileTr = "M" + M + "_2D_Tracking_488.csv";
                outName = "M" + M + "_SurfaceSliceSummary_488.csv";

            elseif contains(fname,"_561")
                tok = regexp(fname,"M(\d+)_2D_561\.csv",'tokens','once');
                M = str2double(tok{1});
                suffix = "_561";
                file2D = "M" + M + "_2D_561.csv";
                file3D = "M" + M + "_3D_561.csv";
                fileTr = "M" + M + "_2D_Tracking_561.csv";
                outName = "M" + M + "_SurfaceSliceSummary_561.csv";

            else
                tok = regexp(fname,"M(\d+)_2D\.csv",'tokens','once');
                M = str2double(tok{1});
                suffix = "";
                file2D = "M" + M + "_2D.csv";
                file3D = "M" + M + "_3D.csv";
                fileTr = "M" + M + "_2D_Tracking.csv";
                outName = "M" + M + "_SurfaceSliceSummary.csv";
            end

            path2D = fullfile(sumFolder,file2D);
            path3D = fullfile(sumFolder,file3D);
            pathTr = fullfile(sumFolder,fileTr);

            if ~isfile(path2D) || ~isfile(path3D) || ~isfile(pathTr)
                continue
            end

            T2 = readtable(path2D);
            T3 = readtable(path3D);
            Ttr = readtable(pathTr);

            % Restrict to long-lived, barely-moving tracks as a proxy for
            % "real, stationary surface-attached objects" -- these are
            % what the track-based surface estimates below are built on.
            Ttr = Ttr(Ttr.TrackDuration > 10 & Ttr.TrackDisplacement < 2,:);

            allT = unique([T2.T; T3.T; Ttr.T]);

            Tout = table();
            Tout.T = allT;
            Tout.SurfaceSliceby2D = nan(length(allT),1);
            Tout.SurfaceSliceby3D = nan(length(allT),1);
            Tout.SurfaceSlicebyTracks = nan(length(allT),1);
            Tout.SurfaceSlicebyTrackingArea = nan(length(allT),1);

            for t = 1:length(allT)

                thisT = allT(t);

                % Method 1: Z-slice with the largest summed 2D area.
                idx2 = T2.T == thisT;
                if any(idx2)
                    z2 = T2.Z(idx2);
                    a2 = T2.Area(idx2);
                    u2 = unique(z2);
                    area2 = zeros(length(u2),1);
                    for j = 1:length(u2)
                        area2(j) = sum(a2(z2 == u2(j)));
                    end
                    [~,mx2] = max(area2);
                    Tout.SurfaceSliceby2D(t) = u2(mx2);
                end

                % Method 2: Z-slice with the largest summed 3D volume.
                idx3 = T3.T == thisT;
                if any(idx3)
                    z3 = T3.Z(idx3);
                    a3 = T3.Volume(idx3);
                    u3 = unique(z3);
                    area3 = zeros(length(u3),1);
                    for j = 1:length(u3)
                        area3(j) = sum(a3(z3 == u3(j)));
                    end
                    [~,mx3] = max(area3);
                    Tout.SurfaceSliceby3D(t) = u3(mx3);
                end

                % Method 3: Z-slice with the most qualifying stable
                % tracked spots (simple count).
                idxTr = Ttr.T == thisT;
                if any(idxTr)
                    zTr = Ttr.Z(idxTr);
                    uTr = unique(zTr);
                    nTr = zeros(length(uTr),1);
                    for j = 1:length(uTr)
                        nTr(j) = sum(zTr == uTr(j));
                    end
                    [~,mxTr] = max(nTr);
                    Tout.SurfaceSlicebyTracks(t) = uTr(mxTr);
                end

                % Method 4: same as method 3, but weighted by the total
                % area of those stable tracked spots rather than a count.
                idxTr = Ttr.T == thisT;
                if any(idxTr)
                    zTr = Ttr.Z(idxTr);
                    aTr = Ttr.Area(idxTr);
                    uTr = unique(zTr);
                    areaTr = zeros(length(uTr),1);
                    for j = 1:length(uTr)
                        areaTr(j) = sum(aTr(zTr == uTr(j)));
                    end
                    [~,mxTrA] = max(areaTr);
                    Tout.SurfaceSlicebyTrackingArea(t) = uTr(mxTrA);
                end

            end

            outPath = fullfile(sumFolder,outName);
            writetable(Tout,outPath)

            % Also keep a copy in memory, nested by culture/date/M so the
            % whole run's surface-slice estimates are available together
            % (names sanitized since e.g. "Protein addition" isn't a
            % valid struct field name as-is).
            cultureKey = matlab.lang.makeValidName(cultureName);
            dateKey = matlab.lang.makeValidName(dateName);
            saveKey = "M" + M + suffix;

            if ~isfield(surfaceSliceSummary,cultureKey)
                surfaceSliceSummary.(cultureKey) = struct();
            end

            if ~isfield(surfaceSliceSummary.(cultureKey),dateKey)
                surfaceSliceSummary.(cultureKey).(dateKey) = struct();
            end

            surfaceSliceSummary.(cultureKey).(dateKey).(matlab.lang.makeValidName(saveKey)) = Tout;

        end
    end
end

%% ===============================
%% FROM STACKS SURFACE + MIDDLE CSVs
%% ===============================
% Third independent pass: using the track-count-based surface-slice
% estimate (SurfaceSlicebyTracks) from the previous section -- smoothed
% over time -- re-slices every 2D/3D/tracking/coarse/watershed table down
% to just that surface Z, and separately to a +/-5-slice band around an
% estimated "middle of the biofilm" Z, writing a parallel CSV for each.

for c = 1:length(cultureTypes)

    cultureName = cultureTypes(c);
    culturePath = fullfile(parentPath,cultureName);

    dateFolders = dir(culturePath);
    dateFolders = dateFolders([dateFolders.isdir]);
    dateFolders = dateFolders(~ismember({dateFolders.name},{'.','..','SurfacePlots'}));

    for d = 1:length(dateFolders)

        dateName = string(dateFolders(d).name);
        datePath = fullfile(culturePath,dateName);
        sumFolder = fullfile(datePath,"Summary");

        if ~isfolder(sumFolder)
            continue
        end

        sliceFiles = dir(fullfile(sumFolder,"M*_SurfaceSliceSummary*.csv"));

        for f = 1:length(sliceFiles)

            sliceName = sliceFiles(f).name;

            % Builds, per channel suffix, the full set of input table
            % names to re-slice (every combination of 2D/3D x plain/
            % Watershed/Coarse, plus the Tracking-derived 2D/3D tables)
            % and the matching Surface_FromStacks_*/Middle_FromStacks_*
            % output names for each.
            if contains(sliceName,"_488")
                tok = regexp(sliceName,"M(\d+)_SurfaceSliceSummary_488\.csv",'tokens','once');
                M = str2double(tok{1});
                suffix = "_488";
                sliceFile = "M" + M + "_SurfaceSliceSummary_488.csv";

                in2D   = "M" + M + "_2D_488.csv";
                in2Dw  = "M" + M + "_2D_Watershed_488.csv";
                in2Dc  = "M" + M + "_2D_Coarse_488.csv";
                in3D   = "M" + M + "_3D_488.csv";
                in3Dw  = "M" + M + "_3D_Watershed_488.csv";
                in3Dc  = "M" + M + "_3D_Coarse_488.csv";
                inTr   = "M" + M + "_2D_Tracking_488.csv";
                inTrw  = "M" + M + "_2D_Tracking_Watershed_488.csv";
                in3DTr = "M" + M + "_3D_Tracking_488.csv";
                in3DTrw = "M" + M + "_3D_Tracking_Watershed_488.csv";

                outSurf2D   = "M" + M + "_Surface_FromStacks_2D_488.csv";
                outSurf2Dw  = "M" + M + "_Surface_FromStacks_2D_Watershed_488.csv";
                outSurf2Dc  = "M" + M + "_Surface_FromStacks_2D_Coarse_488.csv";
                outSurf3D   = "M" + M + "_Surface_FromStacks_3D_488.csv";
                outSurf3Dw  = "M" + M + "_Surface_FromStacks_3D_Watershed_488.csv";
                outSurf3Dc  = "M" + M + "_Surface_FromStacks_3D_Coarse_488.csv";
                outSurfTr   = "M" + M + "_Surface_FromStacks_2D_Tracking_488.csv";
                outSurfTrw  = "M" + M + "_Surface_FromStacks_2D_Tracking_Watershed_488.csv";
                outSurf3DTr = "M" + M + "_Surface_FromStacks_3D_Tracking_488.csv";
                outSurf3DTrw = "M" + M + "_Surface_FromStacks_3D_Tracking_Watershed_488.csv";

                outMid2D   = "M" + M + "_Middle_FromStacks_2D_488.csv";
                outMid2Dw  = "M" + M + "_Middle_FromStacks_2D_Watershed_488.csv";
                outMid2Dc  = "M" + M + "_Middle_FromStacks_2D_Coarse_488.csv";
                outMid3D   = "M" + M + "_Middle_FromStacks_3D_488.csv";
                outMid3Dw  = "M" + M + "_Middle_FromStacks_3D_Watershed_488.csv";
                outMid3Dc  = "M" + M + "_Middle_FromStacks_3D_Coarse_488.csv";
                outMidTr   = "M" + M + "_Middle_FromStacks_2D_Tracking_488.csv";
                outMidTrw  = "M" + M + "_Middle_FromStacks_2D_Tracking_Watershed_488.csv";
                outMid3DTr = "M" + M + "_Middle_FromStacks_3D_Tracking_488.csv";
                outMid3DTrw = "M" + M + "_Middle_FromStacks_3D_Tracking_Watershed_488.csv";

            elseif contains(sliceName,"_561")
                tok = regexp(sliceName,"M(\d+)_SurfaceSliceSummary_561\.csv",'tokens','once');
                M = str2double(tok{1});
                suffix = "_561";
                sliceFile = "M" + M + "_SurfaceSliceSummary_561.csv";

                in2D   = "M" + M + "_2D_561.csv";
                in2Dw  = "M" + M + "_2D_Watershed_561.csv";
                in2Dc  = "M" + M + "_2D_Coarse_561.csv";
                in3D   = "M" + M + "_3D_561.csv";
                in3Dw  = "M" + M + "_3D_Watershed_561.csv";
                in3Dc  = "M" + M + "_3D_Coarse_561.csv";
                inTr   = "M" + M + "_2D_Tracking_561.csv";
                inTrw  = "M" + M + "_2D_Tracking_Watershed_561.csv";
                in3DTr = "M" + M + "_3D_Tracking_561.csv";
                in3DTrw = "M" + M + "_3D_Tracking_Watershed_561.csv";

                outSurf2D   = "M" + M + "_Surface_FromStacks_2D_561.csv";
                outSurf2Dw  = "M" + M + "_Surface_FromStacks_2D_Watershed_561.csv";
                outSurf2Dc  = "M" + M + "_Surface_FromStacks_2D_Coarse_561.csv";
                outSurf3D   = "M" + M + "_Surface_FromStacks_3D_561.csv";
                outSurf3Dw  = "M" + M + "_Surface_FromStacks_3D_Watershed_561.csv";
                outSurf3Dc  = "M" + M + "_Surface_FromStacks_3D_Coarse_561.csv";
                outSurfTr   = "M" + M + "_Surface_FromStacks_2D_Tracking_561.csv";
                outSurfTrw  = "M" + M + "_Surface_FromStacks_2D_Tracking_Watershed_561.csv";
                outSurf3DTr = "M" + M + "_Surface_FromStacks_3D_Tracking_561.csv";
                outSurf3DTrw = "M" + M + "_Surface_FromStacks_3D_Tracking_Watershed_561.csv";

                outMid2D   = "M" + M + "_Middle_FromStacks_2D_561.csv";
                outMid2Dw  = "M" + M + "_Middle_FromStacks_2D_Watershed_561.csv";
                outMid2Dc  = "M" + M + "_Middle_FromStacks_2D_Coarse_561.csv";
                outMid3D   = "M" + M + "_Middle_FromStacks_3D_561.csv";
                outMid3Dw  = "M" + M + "_Middle_FromStacks_3D_Watershed_561.csv";
                outMid3Dc  = "M" + M + "_Middle_FromStacks_3D_Coarse_561.csv";
                outMidTr   = "M" + M + "_Middle_FromStacks_2D_Tracking_561.csv";
                outMidTrw  = "M" + M + "_Middle_FromStacks_2D_Tracking_Watershed_561.csv";
                outMid3DTr = "M" + M + "_Middle_FromStacks_3D_Tracking_561.csv";
                outMid3DTrw = "M" + M + "_Middle_FromStacks_3D_Tracking_Watershed_561.csv";

            else
                tok = regexp(sliceName,"M(\d+)_SurfaceSliceSummary\.csv",'tokens','once');
                M = str2double(tok{1});
                suffix = "";
                sliceFile = "M" + M + "_SurfaceSliceSummary.csv";

                in2D   = "M" + M + "_2D.csv";
                in2Dw  = "M" + M + "_2D_Watershed.csv";
                in2Dc  = "M" + M + "_2D_Coarse.csv";
                in3D   = "M" + M + "_3D.csv";
                in3Dw  = "M" + M + "_3D_Watershed.csv";
                in3Dc  = "M" + M + "_3D_Coarse.csv";
                inTr   = "M" + M + "_2D_Tracking.csv";
                inTrw  = "M" + M + "_2D_Tracking_Watershed.csv";
                in3DTr = "M" + M + "_3D_Tracking.csv";
                in3DTrw = "M" + M + "_3D_Tracking_Watershed.csv";

                outSurf2D   = "M" + M + "_Surface_FromStacks_2D.csv";
                outSurf2Dw  = "M" + M + "_Surface_FromStacks_2D_Watershed.csv";
                outSurf2Dc  = "M" + M + "_Surface_FromStacks_2D_Coarse.csv";
                outSurf3D   = "M" + M + "_Surface_FromStacks_3D.csv";
                outSurf3Dw  = "M" + M + "_Surface_FromStacks_3D_Watershed.csv";
                outSurf3Dc  = "M" + M + "_Surface_FromStacks_3D_Coarse.csv";
                outSurfTr   = "M" + M + "_Surface_FromStacks_2D_Tracking.csv";
                outSurfTrw  = "M" + M + "_Surface_FromStacks_2D_Tracking_Watershed.csv";
                outSurf3DTr = "M" + M + "_Surface_FromStacks_3D_Tracking.csv";
                outSurf3DTrw = "M" + M + "_Surface_FromStacks_3D_Tracking_Watershed.csv";

                outMid2D   = "M" + M + "_Middle_FromStacks_2D.csv";
                outMid2Dw  = "M" + M + "_Middle_FromStacks_2D_Watershed.csv";
                outMid2Dc  = "M" + M + "_Middle_FromStacks_2D_Coarse.csv";
                outMid3D   = "M" + M + "_Middle_FromStacks_3D.csv";
                outMid3Dw  = "M" + M + "_Middle_FromStacks_3D_Watershed.csv";
                outMid3Dc  = "M" + M + "_Middle_FromStacks_3D_Coarse.csv";
                outMidTr   = "M" + M + "_Middle_FromStacks_2D_Tracking.csv";
                outMidTrw  = "M" + M + "_Middle_FromStacks_2D_Tracking_Watershed.csv";
                outMid3DTr = "M" + M + "_Middle_FromStacks_3D_Tracking.csv";
                outMid3DTrw = "M" + M + "_Middle_FromStacks_3D_Tracking_Watershed.csv";
            end

            slicePath = fullfile(sumFolder,sliceFile);

            if ~isfile(slicePath)
                continue
            end

            Tslice = readtable(slicePath);

            % Commit to the track-count-based estimate as the surface Z
            % to use from here on, then smooth it over time: back-fill a
            % missing first value from the first valid one, and for every
            % later timepoint replace a missing value OR an implausible
            % >5-slice jump with the previous timepoint's value.
            surfaceSlice = Tslice.SurfaceSlicebyTracks;

            if ~isempty(surfaceSlice)

                if isnan(surfaceSlice(1))
                    firstGood = find(~isnan(surfaceSlice),1,'first');
                    if ~isempty(firstGood)
                        surfaceSlice(1:firstGood-1) = surfaceSlice(firstGood);
                    end
                end

                for r = 2:length(surfaceSlice)
                    if isnan(surfaceSlice(r))
                        surfaceSlice(r) = surfaceSlice(r-1);
                    elseif abs(surfaceSlice(r) - surfaceSlice(r-1)) > 5
                        surfaceSlice(r) = surfaceSlice(r-1);
                    end
                end
            end

            % "Middle of the biofilm" is approximated as simply half the
            % surface Z index (treating the coverslip as Z=0), not from
            % any independent mid-biofilm signal.
            middleSlice = round(surfaceSlice/2);

            fileList = { ...
                in2D,    outSurf2D,    outMid2D,    "2D"; ...
                in2Dw,   outSurf2Dw,   outMid2Dw,   "2D"; ...
                in2Dc,   outSurf2Dc,   outMid2Dc,   "2D"; ...
                in3D,    outSurf3D,    outMid3D,    "3D"; ...
                in3Dw,   outSurf3Dw,   outMid3Dw,   "3D"; ...
                in3Dc,   outSurf3Dc,   outMid3Dc,   "3D"; ...
                inTr,    outSurfTr,    outMidTr,    "2D"; ...
                inTrw,   outSurfTrw,   outMidTrw,   "2D"; ...
                in3DTr,  outSurf3DTr,  outMid3DTr,  "3D"; ...
                in3DTrw, outSurf3DTrw, outMid3DTrw, "3D"};

            for q = 1:size(fileList,1)

                inName      = fileList{q,1};
                outSurfName = fileList{q,2};
                outMidName  = fileList{q,3};
                modeType    = fileList{q,4};

                inPath = fullfile(sumFolder,inName);

                if ~isfile(inPath)
                    continue
                end

                Tin = readtable(inPath);

                ToutSurf = Tin([],:);
                ToutMid  = Tin([],:);

                allT = unique(Tin.T);

                for tt = 1:length(allT)

                    thisT = allT(tt);

                    idxSlice = Tslice.T == thisT;
                    if ~any(idxSlice)
                        continue
                    end

                    sSurf = surfaceSlice(idxSlice);
                    sMid  = middleSlice(idxSlice);

                    if isnan(sSurf)
                        continue
                    end

                    idxT = Tin.T == thisT;

                    % 2D tables already have an integer slice Z; 3D
                    % tables have a (near-)continuous centroid Z, so it's
                    % rounded before comparing to the estimated surface/
                    % middle slice.
                    if modeType == "2D"
                        idxSurf = idxT & Tin.Z == sSurf;

                        zMin = max(1, sMid - 5);
                        zMax = sMid + 5;
                        idxMid = idxT & Tin.Z >= zMin & Tin.Z <= zMax;
                    else
                        idxSurf = idxT & round(Tin.Z) == sSurf;

                        zMin = max(1, sMid - 5);
                        zMax = sMid + 5;
                        idxMid = idxT & round(Tin.Z) >= zMin & round(Tin.Z) <= zMax;
                    end

                    ToutSurf = [ToutSurf; Tin(idxSurf,:)];
                    ToutMid  = [ToutMid; Tin(idxMid,:)];
                end

                outSurfPath = fullfile(sumFolder,outSurfName);
                outMidPath  = fullfile(sumFolder,outMidName);

                writetable(ToutSurf,outSurfPath)
                writetable(ToutMid,outMidPath)

            end

        end
    end
end
