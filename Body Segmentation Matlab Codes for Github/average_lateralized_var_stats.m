%average_lateralized_var_stats
% Script to aggregate Excel files by ROI.
% For each unique ROI (from the 'ROI' column in each file):
%   - Averages all numeric columns (except SubjectID) across files.
%   - Reports the number of subjects that have that ROI.
%   - Creates one column per subject (e.g. Sub_1, Sub_2, ...) with a 1/0 indicator.

%% Configuration
inputFolder = 'D:\ML_project\Variance\var_excel\updated_sanitized_allmodels\collapsed_randoms';  % <-- Set this to your folder

% Get list of all Excel files in the folder.
files = dir(fullfile(inputFolder, '*.xlsx'));

allTables = {};      % To store each file's table.
allSubjectIDs = [];  % To record subject IDs from filenames.

%% Loop through each Excel file
for k = 1:length(files)
    fileName = files(k).name;
    fullPath = fullfile(inputFolder, fileName);
    
    % Extract subject ID from the filename.
    % This assumes the only number in the filename is the subject's ID.
    token = regexp(fileName, '\d+', 'match');
    if isempty(token)
        warning('No subject ID found in file: %s. Skipping.', fileName);
        continue;
    end
    subjID = str2double(token{1});
    allSubjectIDs(end+1) = subjID; %#ok<AGROW>
    
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

% Get the union of all subject IDs.
uniqueSubjects = unique(combinedTable.SubjectID);
numSubjects = length(uniqueSubjects);

% Make sure the ROI column is categorical for grouping.
if ~iscategorical(combinedTable.ROI)
    combinedTable.ROI = categorical(combinedTable.ROI);
end

%% Identify Numeric Columns to Average
allVarNames = combinedTable.Properties.VariableNames;
% Determine which variables are numeric.
isNumeric = varfun(@isnumeric, combinedTable, 'OutputFormat', 'uniform');
% Exclude 'SubjectID' (and ROI if it were numeric, but it should be categorical).
numericVars = setdiff(allVarNames(isNumeric), {'SubjectID'});

%% Group the Data by ROI
[grp, roiList] = findgroups(combinedTable.ROI);
numROIs = length(roiList);

% Preallocate for averaged numeric data and for counting subjects.
avgData = cell(numROIs, length(numericVars));
nSubjectsCol = zeros(numROIs, 1);
subjectIndicator = zeros(numROIs, numSubjects);  % Each row for an ROI, each column for one subject.

% Loop over each unique ROI group.
for i = 1:numROIs
    idx = (grp == i);
    Tgroup = combinedTable(idx, :);
    
    % Average each numeric variable for this ROI.
    for j = 1:length(numericVars)
        avgData{i,j} = mean(Tgroup.(numericVars{j}));
    end
    
    % Get the unique subject IDs that contributed to this ROI.
    subjIDsForROI = unique(Tgroup.SubjectID);
    nSubjectsCol(i) = length(subjIDsForROI);
    
    % Create a binary indicator for each subject.
    for s = 1:numSubjects
        if ismember(uniqueSubjects(s), subjIDsForROI)
            subjectIndicator(i, s) = 1;
        else
            subjectIndicator(i, s) = 0;
        end
    end
end

%% Create Output Table
resultTable = table;
% ROI column as character array.
resultTable.ROI = cellstr(roiList);

% For each numeric variable, add a column with the averaged value.
for j = 1:length(numericVars)
    colName = sprintf('Avg_%s', numericVars{j});
    resultTable.(colName) = cell2mat(avgData(:, j));
end

% Add a column for the number of subjects that provided data for the ROI.
resultTable.nSubjects = nSubjectsCol;

% Add one column per subject with binary indicators.
for s = 1:numSubjects
    subColName = sprintf('Sub_%d', uniqueSubjects(s));
    resultTable.(subColName) = subjectIndicator(:, s);
end

%% Write Output Table to Excel
outputFile = fullfile(inputFolder, 'Aggregated_ROI_Results.xlsx');
writetable(resultTable, outputFile);
fprintf('Aggregated results saved to %s\n', outputFile);
