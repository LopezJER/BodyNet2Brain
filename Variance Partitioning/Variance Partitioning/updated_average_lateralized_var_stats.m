%updated_average_lateralized_var_stats

% Aggregated ROI Results - Mean and Variance
% This script loops over a set of Excel files in an input folder,
% groups the data by ROI (from the 'ROI' column),
% computes the average (mean) and variance for each numeric variable (excluding SubjectID),
% and writes the results to an output Excel file.
% The output table will list all average columns first, followed by all variance columns.

tic;  % Start timer
clear; clc;

%% Configuration
inputFolder = 'D:\ML_project\Variance\var_excel\sapiens_allmodels';  % Set your folder path

% Get list of all Excel files in the folder.
files = dir(fullfile(inputFolder, '*.xlsx'));

allTables = {};      % To store each file's table

%% Loop through each Excel file
for k = 1:length(files)
    fileName = files(k).name;
    fullPath = fullfile(inputFolder, fileName);
    
    % Extract subject ID from the filename (assumes the only number in the filename is the subject's ID)
    token = regexp(fileName, '\d+', 'match');
    if isempty(token)
        warning('No subject ID found in file: %s. Skipping.', fileName);
        continue;
    end
    subjID = str2double(token{1});
    
    % Read the Excel file.
    T = readtable(fullPath);
    
    % Ensure the table has an 'ROI' column.
    if ~ismember('ROI', T.Properties.VariableNames)
        warning('File %s does not contain a "ROI" column. Skipping.', fileName);
        continue;
    end
    
    % Add a SubjectID column.
    T.SubjectID = repmat(subjID, height(T), 1);
    
    % Store the table.
    allTables{end+1} = T;
end

if isempty(allTables)
    error('No valid Excel files processed.');
end

% Combine all tables into one.
combinedTable = vertcat(allTables{:});

% Make sure the ROI column is categorical for grouping.
if ~iscategorical(combinedTable.ROI)
    combinedTable.ROI = categorical(combinedTable.ROI);
end

%% Identify Numeric Columns to Process in Original Order
allVarNames = combinedTable.Properties.VariableNames;
numericVars = {};  % This will store numeric variable names in the order they appear.
for idx = 1:length(allVarNames)
    if isnumeric(combinedTable.(allVarNames{idx})) && ~strcmp(allVarNames{idx}, 'SubjectID')
        numericVars{end+1} = allVarNames{idx}; %#ok<AGROW>
    end
end

%% Group the Data by ROI
[grp, roiList] = findgroups(combinedTable.ROI);
numROIs = length(roiList);

% Preallocate cell arrays for means and variances.
meanData = cell(numROIs, length(numericVars));
varData  = cell(numROIs, length(numericVars));

% Loop over each ROI group.
for i = 1:numROIs
    idx = (grp == i);
    Tgroup = combinedTable(idx, :);
    for j = 1:length(numericVars)
        colName = numericVars{j};
        meanData{i,j} = mean(Tgroup.(colName), 'omitnan');
        varData{i,j}  = var(Tgroup.(colName), 'omitnan');
    end
end

%% Create Output Table
resultTable = table;
resultTable.ROI = cellstr(roiList);

% First, add all average columns in the original order.
for j = 1:length(numericVars)
    meanColName = sprintf('Avg_%s', numericVars{j});
    resultTable.(meanColName) = cell2mat(meanData(:, j));
end

% Then, add all variance columns (appended after the avg columns).
for j = 1:length(numericVars)
    varColName = sprintf('Var_%s', numericVars{j});
    resultTable.(varColName) = cell2mat(varData(:, j));
end

%% Write Output Table to Excel
outputFile = fullfile(inputFolder, 'Aggregated_ROI_Results_updated.xlsx');
writetable(resultTable, outputFile);
fprintf('Aggregated results saved to %s\n', outputFile);

toc;  % Display elapsed time
