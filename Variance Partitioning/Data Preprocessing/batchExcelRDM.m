%% ========= CONFIGURATION SECTION =========
% Set the base directory where subject folders (named like "subj*") are located.
% For example, if your subject folders "subj1", "subj2", etc. are located in 
% D:\ML_project\RDM_results\final then set BASE_SUBJECT_DIR to that folder.
BASE_SUBJECT_DIR = 'D:\ML_project\RDM_results\final';  % <-- Change this if needed.

% Optionally, limit processing to specific subject numbers.
% For example, to run only subject 3, set: subjectIndices = [3];
subjectIndices = [1];  % Leave empty ([]) to process all subject folders.
%% ===========================================

%% ========= BATCH PROCESSING CODE =========
% This script will:
%   1. Recursively search each subject folder for .xlsx files.
%   2. For each file, determine the “region” as the name of the file’s parent folder.
%   3. If the Excel file’s base name is long (>= 10 characters), assume it’s an
%      original aggregated file:
%         a. If a formatted file ([region].xlsx) does not exist, process it via processExcelData.
%         b. Then, if a PNG ([region_RDM.png]) does not exist, run createRDMFromExcel on the formatted file.
%   4. Otherwise (if the file’s base name is short), assume it is already formatted and
%      use it for RDM creation if needed.
%   5. Record summary info for each step.
%
% NOTE: This version “feeds” the RDM stage from the formatting stage so that files
% processed in step 1 automatically are used for step 2.

% Initialize summary storage.
allStats = struct('subjectNum', {}, 'regionName', {}, 'step', {}, 'status', {}, 'inputFile', {}, 'outputFile', {});

% List all subject folders (assumed to be named like "subj*") in BASE_SUBJECT_DIR.
subjectFolders = dir(fullfile(BASE_SUBJECT_DIR, 'subj*'));
fprintf('[DEBUG] Found %d subject folders in %s\n', length(subjectFolders), BASE_SUBJECT_DIR);

for k = 1:length(subjectFolders)
    if subjectFolders(k).isdir
        % Extract subject number from folder name (e.g., "subj3" -> 3)
        subjectNum = sscanf(subjectFolders(k).name, 'subj%d');
        if ~isempty(subjectIndices) && ~ismember(subjectNum, subjectIndices)
            fprintf('[DEBUG] Skipping subject %d as it is not in subjectIndices.\n', subjectNum);
            continue;
        end
        
        % Full path of the subject folder.
        subjectFolderPath = fullfile(subjectFolders(k).folder, subjectFolders(k).name);
        fprintf('[DEBUG] Processing subject %d in folder: %s\n', subjectNum, subjectFolderPath);
        
        % Recursively find all .xlsx files within this subject folder.
        xlsxFiles = dir(fullfile(subjectFolderPath, '**', '*.xlsx'));
        fprintf('[DEBUG] Found %d .xlsx files in subject %d\n', length(xlsxFiles), subjectNum);
        
        for i = 1:length(xlsxFiles)
            inputFile = fullfile(xlsxFiles(i).folder, xlsxFiles(i).name);
            [~, baseName, ~] = fileparts(inputFile);
            % Use the current folder's name as the region name.
            [~, region] = fileparts(xlsxFiles(i).folder);
            fprintf('[DEBUG] Considering file: %s (base name: "%s", length: %d) in folder: %s (region: %s)\n', ...
                inputFile, baseName, length(baseName), xlsxFiles(i).folder, region);
            
            % Define output names:
            % Formatted Excel file will be named [region].xlsx in the same folder.
            formattedFile = fullfile(xlsxFiles(i).folder, [region, '.xlsx']);
            % RDM PNG will be named [region_RDM.png]
            pngFile = fullfile(xlsxFiles(i).folder, [region, '_RDM.png']);
            
            if length(baseName) >= 10
                % -------------------------------
                % Formatting branch (original aggregated file)
                % -------------------------------
                % Check if a "sub10" folder exists in the current folder; if so, skip formatting.
                currentFolder = xlsxFiles(i).folder;
                subFolders = dir(currentFolder);
                isDir = [subFolders.isdir];
                subFolderNames = {subFolders(isDir).name};
                subFolderNames = subFolderNames(~ismember(subFolderNames, {'.', '..'}));
                if any(strcmp(subFolderNames, 'sub10'))
                    fprintf('[DEBUG] Skipping formatting for %s because "sub10" folder exists in %s\n', inputFile, currentFolder);
                    allStats(end+1) = struct('subjectNum', subjectNum, 'regionName', region, 'step', 'formatting', 'status', 'skipped (sub10 exists)', 'inputFile', inputFile, 'outputFile', formattedFile);
                else
                    if ~isfile(formattedFile)
                        fprintf('[DEBUG] Formatting Excel for subject %d, region "%s"\n', subjectNum, region);
                        processExcelData(inputFile, formattedFile);
                        allStats(end+1) = struct('subjectNum', subjectNum, 'regionName', region, 'step', 'formatting', 'status', 'formatted', 'inputFile', inputFile, 'outputFile', formattedFile);
                    else
                        fprintf('[DEBUG] Skipping formatting for %s because processed file already exists: %s\n', inputFile, formattedFile);
                        allStats(end+1) = struct('subjectNum', subjectNum, 'regionName', region, 'step', 'formatting', 'status', 'skipped (output exists)', 'inputFile', inputFile, 'outputFile', formattedFile);
                    end
                end
                % Now, after (or if skipped) formatting, run the RDM stage on the formatted file.
                if ~isfile(pngFile)
                    fprintf('[DEBUG] Creating RDM for subject %d, region "%s" using formatted file: %s\n', subjectNum, region, formattedFile);
                    createRDMFromExcel(formattedFile);
                    allStats(end+1) = struct('subjectNum', subjectNum, 'regionName', region, 'step', 'RDM creation', 'status', 'RDM created', 'inputFile', formattedFile, 'outputFile', pngFile);
                else
                    fprintf('[DEBUG] Skipping RDM creation for %s because PNG already exists: %s\n', formattedFile, pngFile);
                    allStats(end+1) = struct('subjectNum', subjectNum, 'regionName', region, 'step', 'RDM creation', 'status', 'skipped (PNG exists)', 'inputFile', formattedFile, 'outputFile', pngFile);
                end
            else
                % -------------------------------
                % RDM creation branch (file is already short, assumed formatted)
                % -------------------------------
                if ~isfile(pngFile)
                    fprintf('[DEBUG] Creating RDM for subject %d, region "%s" using file: %s\n', subjectNum, region, inputFile);
                    createRDMFromExcel(inputFile);
                    allStats(end+1) = struct('subjectNum', subjectNum, 'regionName', region, 'step', 'RDM creation', 'status', 'RDM created', 'inputFile', inputFile, 'outputFile', pngFile);
                else
                    fprintf('[DEBUG] Skipping RDM creation for %s because PNG already exists: %s\n', inputFile, pngFile);
                    allStats(end+1) = struct('subjectNum', subjectNum, 'regionName', region, 'step', 'RDM creation', 'status', 'skipped (PNG exists)', 'inputFile', inputFile, 'outputFile', pngFile);
                end
            end
        end
    end
end

%% ========= COMBINED SUMMARY REPORT =========
fprintf('\n\n========== COMBINED SUMMARY REPORT ==========\n');
% Group the summary by subject and region.
uniqueKeys = unique(arrayfun(@(s) sprintf('subj%d_%s', s.subjectNum, s.regionName), allStats, 'UniformOutput', false));
for i = 1:length(uniqueKeys)
    key = uniqueKeys{i};
    idx = find(arrayfun(@(s) strcmp(sprintf('subj%d_%s', s.subjectNum, s.regionName), key), allStats));
    % For this subject-region, extract formatting and RDM creation statuses.
    formattingEntries = allStats(idx(strcmp({allStats(idx).step}, 'formatting')));
    rdmEntries = allStats(idx(strcmp({allStats(idx).step}, 'RDM creation')));
    
    if ~isempty(formattingEntries)
        formattingStatus = formattingEntries(1).status;
    else
        formattingStatus = 'N/A';
    end
    if ~isempty(rdmEntries)
        rdmStatus = rdmEntries(1).status;
    else
        rdmStatus = 'N/A';
    end
    
    subj = allStats(idx(1)).subjectNum;
    reg = allStats(idx(1)).regionName;
    fprintf('  subj%d: %s -> Formatting: %s; RDM creation: %s\n', subj, reg, formattingStatus, rdmStatus);
end
fprintf('===============================================\n');

%% ========= FUNCTIONS =========
% ----- Function: processExcelData -----
function processExcelData(inputFile, outputFile)
    tic;  % Start timer.
    fprintf('Reading input Excel file: %s\n', inputFile);
    try
        data = readcell(inputFile);
    catch ME
        error(['Failed to read the Excel file. Error: ', ME.message]);
    end
    if size(data, 1) < 1
        error('The Excel file must have at least one row to remove.');
    end
    data(1, :) = [];  % Remove the first row.
    if isempty(data)
        error('No data left after removing the first row.');
    end
    data = data';  % Transpose.
    for col = 1:size(data, 2)
        entry = data{1, col};
        if ischar(entry) || isstring(entry)
            data{1, col} = ['x' entry];
        elseif isnumeric(entry)
            data{1, col} = ['x' num2str(entry)];
        else
            data{1, col} = ['x' char(string(entry))];
        end
    end
    for row = 2:size(data, 1)
        for col = 1:size(data, 2)
            currentValue = data{row, col};
            if isnumeric(currentValue) && ~isnan(currentValue)
                data{row, col} = currentValue / 300;
            end
        end
    end
    fprintf('Writing processed data to: %s\n', outputFile);
    try
        writecell(data, outputFile);
    catch ME
        error(['Failed to write the Excel file. Error: ', ME.message]);
    end
    elapsedTime = toc;
    fprintf('Processing complete for %s. Elapsed time: %.2f seconds.\n', inputFile, elapsedTime);
end

% ----- Function: createRDMFromExcel -----
function createRDMFromExcel(excelPath)
    fprintf('Reading Excel file for RDM: %s\n', excelPath);
    data = readtable(excelPath);
    vectors = table2array(data);
    numColumns = size(vectors, 2);
    fprintf('Number of vectors: %d\n', numColumns);
    RDM = zeros(numColumns);
    fprintf('Computing RDM...\n');
    for i = 1:numColumns
        for j = i+1:numColumns
            r = corr(vectors(:, i), vectors(:, j));
            RDM(i, j) = 1 - r;
        end
    end
    RDM = RDM + RDM';  % Make RDM symmetric.
    fprintf('Generating dynamic title...\n');
    dynamicTitle = generateTitle(excelPath);
    fprintf('Displaying RDM figure...\n');
    figure;
    imagesc(RDM);
    colorbar;
    title(dynamicTitle, 'Interpreter', 'none');
    set(gca, 'XTick', [], 'YTick', []);
    set(gca, 'Box', 'off');
    set(gca, 'TickLength', [0 0]);
    set(gca, 'XColor', 'none', 'YColor', 'none');
    axis square;
    [filePath, baseName, ~] = fileparts(excelPath);
    pngFilename = fullfile(filePath, [baseName, '_RDM.png']);
    fprintf('Saving RDM PNG as: %s\n', pngFilename);
    saveas(gcf, pngFilename);
    close(gcf);
    fprintf('RDM saved successfully.\n');
end

% ----- Helper Function: generateTitle -----
function titleStr = generateTitle(excelPath)
    subjectNumber = 'Unknown';
    side = '';
    region = 'Unknown';
    subjPattern = 'subj(\d+)';
    subjMatch = regexp(excelPath, subjPattern, 'tokens', 'ignorecase');
    if ~isempty(subjMatch)
        subjectNumber = subjMatch{1}{1};
    end
    [~, fileName, ~] = fileparts(excelPath);
    if startsWith(fileName, {'r', 'R'}, 'IgnoreCase', true)
        side = 'Right';
        region = extractAfter(fileName, 1);
    elseif startsWith(fileName, {'l', 'L'}, 'IgnoreCase', true)
        side = 'Left';
        region = extractAfter(fileName, 1);
    else
        region = fileName;
    end
    titleStr = sprintf('Subject %s %s %s RDM', subjectNumber, side, region);
    titleStr = strtrim(titleStr);
end
