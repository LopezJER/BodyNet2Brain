function splitExcelBySubject(inputExcelPath, outputFolderPath)
% splitExcelBySubject: Splits an input Excel file into separate files for each subject.
%   Keeps only relevant columns (subject repetition and cocoId) and removes rows
%   with no data in repetition columns for each subject.
% 
% Inputs:
%   inputExcelPath - Path to the input Excel file.
%   outputFolderPath - Path to the folder where subject-specific Excel files will be saved.

    % Read the input Excel file
    data = readtable(inputExcelPath);

    % Extract subject and repetition columns
    columnNames = data.Properties.VariableNames;
    cocoIDColumn = 'cocoId';

    % Ensure the cocoId column exists
    if ~ismember(cocoIDColumn, columnNames)
        error('The input Excel file must contain a column named "%s".', cocoIDColumn);
    end

    % Determine unique subject numbers from repetition columns
    repColumns = columnNames(contains(columnNames, '_rep'));
    subjectNumbers = unique(cellfun(@(x) sscanf(x, 'subject%d_rep%d', 1), repColumns, 'UniformOutput', true));

    % Create the output folder if it doesn't exist
    if ~exist(outputFolderPath, 'dir')
        mkdir(outputFolderPath);
    end

    % Process each subject
    for subjectID = subjectNumbers
        fprintf('Processing Subject %d...\n', subjectID);

        % Identify columns for the current subject
        subjectRepColumns = repColumns(contains(repColumns, sprintf('subject%d_rep', subjectID)));

        % Select cocoId and subject-specific repetition columns
        relevantColumns = [cocoIDColumn, subjectRepColumns];
        subjectData = data(:, relevantColumns);

        % Remove rows where all repetition columns are NaN
        validRows = any(~isnan(table2array(subjectData(:, subjectRepColumns))), 2);
        subjectData = subjectData(validRows, :);

        % Save the processed data to a new Excel file
        outputFileName = sprintf('subject%d_data.xlsx', subjectID);
        outputFilePath = fullfile(outputFolderPath, outputFileName);
        writetable(subjectData, outputFilePath);

        fprintf('  Saved Subject %d data to %s\n', subjectID, outputFilePath);
    end

    fprintf('Processing complete. Files saved in %s\n', outputFolderPath);
end
