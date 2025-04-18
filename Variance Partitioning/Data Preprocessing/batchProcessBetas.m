%% ========= CONFIGURATION SECTION =========
% Base directories for the betas processing pipeline:
BASE_ROI_DIR       = 'D:\ML_project\ROIs';           % Where ROI mask folders (subj*) are located.
SUBJECT_TRIALS_DIR = 'D:\ML_project\subject_trials';   % Where input Excel files are (named like subjectX_data.xlsx).
HDF5_BASE_DIR      = 'D:\ML_project\Subject_betas';    % Where HDF5 files (for betas) are stored (subjX folders).
OUTPUT_BASE_DIR    = 'D:\ML_project\RDM_results\final';  % Where aggregated results will be saved (organized by subjX\region).

% To process only specific subjects, set subjectIndices (e.g. [3]). Otherwise, leave empty.
subjectIndices = [8];  % <-- Change as needed.
%% ===========================================

%% ========= BATCH PROCESSING CODE =========
% This code:
% 1. Iterates over subject folders in BASE_ROI_DIR (expected to be named "subj*").
% 2. For each subject, it searches the ROI mask subfolders (those ending with _segmented)
%    for .nii files. Each .nii file corresponds to an ROI; the ROI (or region) is taken as
%    the name of the ROI mask file (or its parent folder, as desired).
% 3. For each ROI mask, it uses the subject number to build paths for:
%    - The input Excel file (subjectX_data.xlsx)
%    - The HDF5 folder (for that subject)
%    - The output aggregated Excel file (saved under OUTPUT_BASE_DIR\subjX\region\aggregated_results_unformatted.xlsx)
% 4. The processing function (processSubject) performs all steps:
%    a. It loads the Excel data, gathers unique trials,
%    b. Searches for session HDF5 files and builds a trial–to–session map,
%    c. Loads betas from HDF5 files, computes row-wise averages,
%    d. Appends the results to the output Excel.
%
% A summary is recorded for each subject-region (ROI) and then printed at the end.

% Initialize summary storage.
betasSummary = struct('subjectNum', {}, 'region', {}, 'numRowsProcessed', {}, 'numUniqueTrials', {}, 'numMappedTrials', {});

% List all subject folders (named like "subj*") in BASE_ROI_DIR.
subjectFolders = dir(fullfile(BASE_ROI_DIR, 'subj*'));
fprintf('[DEBUG] Found %d subject folders in %s\n', length(subjectFolders), BASE_ROI_DIR);

for k = 1:length(subjectFolders)
    if subjectFolders(k).isdir
        subjectName = subjectFolders(k).name;
        subjectNum = sscanf(subjectName, 'subj%d');
        if ~isempty(subjectIndices) && ~ismember(subjectNum, subjectIndices)
            fprintf('[DEBUG] Skipping subject %d (not in subjectIndices).\n', subjectNum);
            continue;
        end
        
        % ROI mask files are assumed to be within subfolders (e.g., "subj3\*_segmented\*.nii")
        roiRoot = fullfile(subjectFolders(k).folder, subjectName);
        segmentedFolders = dir(fullfile(roiRoot, '*_segmented'));
        for j = 1:length(segmentedFolders)
            if segmentedFolders(j).isdir
                niiFiles = dir(fullfile(segmentedFolders(j).folder, segmentedFolders(j).name, '*.nii'));
                for i = 1:length(niiFiles)
                    roiMaskFile = fullfile(niiFiles(i).folder, niiFiles(i).name);
                    % You can choose to derive the region name from the mask file itself or its folder.
                    % Here we take the region as the file name (without extension):
                    [~, regionName, ~] = fileparts(roiMaskFile);
                    
                    % Build output Excel path.
                    outputExcelPath = fullfile(OUTPUT_BASE_DIR, sprintf('subj%d', subjectNum), regionName, 'aggregated_results_unformatted.xlsx');
                    fprintf('[DEBUG] Considering ROI mask: %s (region: %s)\n', roiMaskFile, regionName);
                    
                    % If output already exists, skip processing.
                    if isfile(outputExcelPath)
                        fprintf('[DEBUG] Skipping ROI %s because output file exists: %s\n', roiMaskFile, outputExcelPath);
                        continue;
                    end
                    
                    fprintf('[DEBUG] Processing ROI mask: %s\n', roiMaskFile);
                    % Process this ROI mask.
                    [nRows, nUnique, nMapped] = processSubject(subjectNum, roiMaskFile, SUBJECT_TRIALS_DIR, HDF5_BASE_DIR, OUTPUT_BASE_DIR);
                    
                    % Record summary info.
                    betasSummary(end+1) = struct('subjectNum', subjectNum, 'region', regionName, 'numRowsProcessed', nRows, 'numUniqueTrials', nUnique, 'numMappedTrials', nMapped);
                end
            end
        end
    end
end

%% ========= SUMMARY REPORT =========
fprintf('\n\n========== BETAS PROCESSING SUMMARY ==========\n');
uniqueKeys = unique(arrayfun(@(s) sprintf('subj%d_%s', s.subjectNum, s.region), betasSummary, 'UniformOutput', false));
for i = 1:length(uniqueKeys)
    key = uniqueKeys{i};
    idx = find(arrayfun(@(s) strcmp(sprintf('subj%d_%s', s.subjectNum, s.region), key), betasSummary));
    subj = betasSummary(idx(1)).subjectNum;
    reg = betasSummary(idx(1)).region;
    totalRows = sum([betasSummary(idx).numRowsProcessed]);
    totalUnique = sum([betasSummary(idx).numUniqueTrials]);
    totalMapped = sum([betasSummary(idx).numMappedTrials]);
    fprintf('  subj%d: %s -> Rows processed: %d, Unique trials: %d, Mapped trials: %d\n', subj, reg, totalRows, totalUnique, totalMapped);
end
fprintf('===============================================\n');

%% ========= FUNCTIONS =========

% ----- Function: processSubject -----
function [numRowsProcessed, numUniqueTrials, numMappedTrials] = processSubject(subjectNum, roiMaskFile, SUBJECT_TRIALS_DIR, HDF5_BASE_DIR, OUTPUT_BASE_DIR)
    % Build paths.
    inputExcelPath  = sprintf('%s\\subject%d_data.xlsx', SUBJECT_TRIALS_DIR, subjectNum);
    [~, regionName, ~] = fileparts(roiMaskFile);
    outputExcelPath = sprintf('%s\\subj%d\\%s\\aggregated_results_unformatted.xlsx', OUTPUT_BASE_DIR, subjectNum, regionName);
    roiMaskPath     = roiMaskFile;
    hdf5FolderPath  = sprintf('%s\\subj%d', HDF5_BASE_DIR, subjectNum);
    
    fprintf('[DEBUG] Subject %d, region %s: Using ROI mask: %s\n', subjectNum, regionName, roiMaskFile);
    
    % --- Load Excel Table ---
    fprintf('[DEBUG] Reading input Excel file: %s\n', inputExcelPath);
    data = readtable(inputExcelPath);
    numRowsProcessed = height(data);
    if isempty(data)
        error('[ERROR] Excel file is empty: %s', inputExcelPath);
    end
    columnNames = data.Properties.VariableNames;
    cocoIdColumn = columnNames{1};  % Assume first column is CocoID.
    trialColumns = columnNames(2:4);  % Assume next three columns are trials.
    
    rowsToProcess = 1:min(99999, height(data));
    data = data(rowsToProcess, :);
    fprintf('[DEBUG] Processing %d rows from Excel file.\n', height(data));
    
    % --- Gather Unique Trials ---
    allTrials = [];
    for iRow = 1:height(data)
        rowTrials = table2array(data(iRow, trialColumns));
        rowTrials = double(rowTrials);
        allTrials = [allTrials; rowTrials];
    end
    uniqueTrials = unique(allTrials);
    numUniqueTrials = numel(uniqueTrials);
    fprintf('[DEBUG] Found %d unique trial IDs in Excel file.\n', numUniqueTrials);
    
    % --- Search for HDF5 Files ---
    sessionFiles_h5   = dir(fullfile(hdf5FolderPath, '*.h5'));
    sessionFiles_hdf5 = dir(fullfile(hdf5FolderPath, '*.hdf5'));
    sessionFiles = [sessionFiles_h5; sessionFiles_hdf5];
    if isempty(sessionFiles)
        error('[ERROR] No HDF5 files found in: %s', hdf5FolderPath);
    end
    maxAvailableTrials = 750 * numel(sessionFiles);
    fprintf('[DEBUG] Total available trials (from HDF5): %d\n', maxAvailableTrials);
    
    % --- Build Trial-to-Session Map ---
    trialToSessionMap = containers.Map('KeyType','double','ValueType','char');
    numMappedTrials = 0;
    for s = 1:length(sessionFiles)
        sessionFilePath = fullfile(sessionFiles(s).folder, sessionFiles(s).name);
        [~, fname, ~] = fileparts(sessionFilePath);
        tokens = regexp(fname, 'betas_session(\d+)$', 'tokens', 'once');
        if isempty(tokens)
            continue;
        end
        sessNum = str2double(tokens{1});
        if isnan(sessNum)
            continue;
        end
        startID = (sessNum - 1)*750 + 1;
        endID   = sessNum*750;
        sessionTrialIDs = startID:endID;
        sessionTrialIDs = sessionTrialIDs(sessionTrialIDs <= maxAvailableTrials);
        for tID = sessionTrialIDs
            if ismember(tID, uniqueTrials)
                trialToSessionMap(tID) = sessionFilePath;
                numMappedTrials = numMappedTrials + 1;
            end
        end
    end
    fprintf('[DEBUG] Mapped %d trial IDs to session files.\n', numMappedTrials);
    
    % --- Load Betas from HDF5 Files and Build (trial -> betas) Map ---
    trialBetasMap = containers.Map('KeyType','double','ValueType','any');
    for s = 1:length(sessionFiles)
        sessionFilePath = fullfile(sessionFiles(s).folder, sessionFiles(s).name);
        [~, fname, ~] = fileparts(sessionFilePath);
        tokens = regexp(fname, 'betas_session(\d+)$', 'tokens', 'once');
        if isempty(tokens)
            continue;
        end
        sessNum = str2double(tokens{1});
        if isnan(sessNum)
            continue;
        end
        startID = (sessNum - 1)*750 + 1;
        endID   = sessNum*750;
        sessionTrialIDs = intersect(startID:endID, uniqueTrials);
        if isempty(sessionTrialIDs)
            continue;
        end
        fprintf('[DEBUG] Loading session file: %s (session %d) for %d trials.\n', sessionFiles(s).name, sessNum, numel(sessionTrialIDs));
        localBetasMap = loadSessionBetas_real(sessionFilePath, roiMaskPath, sessionTrialIDs);
        for iT = 1:numel(sessionTrialIDs)
            tID = sessionTrialIDs(iT);
            trialBetasMap(tID) = localBetasMap(tID);
        end
    end
    fprintf('[DEBUG] Loaded betas; trialBetasMap now has %d keys.\n', length(trialBetasMap.keys));
    
    % --- Compute Row-wise Averages ---
    resultsTable = table();
    fprintf('[DEBUG] Computing row-wise averages...\n');
    for iRow = 1:height(data)
        cocoId = data.(cocoIdColumn)(iRow);
        trials = table2array(data(iRow, trialColumns));
        trials = double(trials);
        trials = trials(trials <= maxAvailableTrials);
        validTrialIDs = trials(ismember(trials, cell2mat(trialBetasMap.keys)));
        if isempty(validTrialIDs)
            fprintf('[WARNING] Row %d (CocoID=%d) has no valid trial IDs; skipping.\n', iRow, cocoId);
            continue;
        end
        avgBetas = zeros(size(trialBetasMap(validTrialIDs(1))));
        for j = 1:length(validTrialIDs)
            avgBetas = avgBetas + trialBetasMap(validTrialIDs(j));
        end
        avgBetas = avgBetas / length(validTrialIDs);
        colNames = arrayfun(@(x) sprintf('Beta%d', x), 1:length(avgBetas), 'UniformOutput', false);
        outRow = array2table(avgBetas(:)', 'VariableNames', colNames);
        outRow = addvars(outRow, cocoId, 'Before', 'Beta1', 'NewVariableNames', 'CocoID');
        resultsTable = [resultsTable; outRow];
    end
    fprintf('[DEBUG] Computed averages; results table has %d rows.\n', height(resultsTable));
    
    % --- Append (or write) results to output Excel file ---
    if isfile(outputExcelPath)
        aggregatedResults = readtable(outputExcelPath, 'PreserveVariableNames', true);
    else
        aggregatedResults = table();
    end
    if isempty(aggregatedResults)
        aggregatedResults = resultsTable;
    else
        aggregatedResults = [aggregatedResults; resultsTable];
    end
    fprintf('[DEBUG] Writing aggregated results to: %s\n', outputExcelPath);
    writetable(aggregatedResults, outputExcelPath);
    fprintf('[DEBUG] Processing complete for subject %d, region %s.\n', subjectNum, regionName);
end

%% ----- Function: loadSessionBetas_real -----
function localBetasMap = loadSessionBetas_real(sessionFilePath, roiMaskPath, sessionTrialIds)
    allBetas = h5read(sessionFilePath, '/betas');
    roiMask = niftiread(roiMaskPath);
    roiMask = (roiMask > 0);
    sizeBetas = size(allBetas);
    sizeMask  = size(roiMask);
    if length(sizeBetas) < 4
        error('[ERROR] /betas dataset not 4D. Found size: %s', mat2str(sizeBetas));
    end
    if any(sizeBetas(1:3) ~= sizeMask)
        warning('[WARNING] Mismatch in mask vs betas dims: betas(1:3)=%s, mask=%s', mat2str(sizeBetas(1:3)), mat2str(sizeMask));
    end
    localBetasMap = containers.Map('KeyType','double','ValueType','any');
    [~, fname, ~] = fileparts(sessionFilePath);
    tokens = regexp(fname, 'betas_session(\d+)$', 'tokens', 'once');
    sessNum = str2double(tokens{1});
    startID = (sessNum - 1)*750 + 1;
    nTrialsInThisSession = sizeBetas(4);
    for i = 1:numel(sessionTrialIds)
        tID = sessionTrialIds(i);
        trialIndex = tID - startID + 1;
        if trialIndex < 1 || trialIndex > nTrialsInThisSession
            warning('[WARNING] Trial %d out of range for file %s; skipping.', tID, sessionFilePath);
            continue;
        end
        volume3D = allBetas(:, :, :, trialIndex);
        maskedVoxels = volume3D(roiMask);
        localBetasMap(tID) = double(maskedVoxels(:));
    end
end
