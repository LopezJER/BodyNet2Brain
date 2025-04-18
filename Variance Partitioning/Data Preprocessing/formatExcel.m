% ============================================================
% MATLAB Function: processExcelData.m
% ============================================================
% This function processes an Excel file by performing the following steps:
% 1. Removes the entire first row.
% 2. Transposes the remaining data (rows become columns and vice versa).
% 3. Divides all numeric values (except those in the first row) by 300.
% 4. Adds an 'x' prefix to all entries in the first row.
% 5. Saves the processed data to a new Excel file.
%
% Usage:
%   processExcelData('input.xlsx', 'output.xlsx');
%
% ============================================================

function processExcelData(inputFile, outputFile)
    % ==================================================
    % Step 0: Start Timing (Optional)
    % ==================================================
    tic;  % Start the timer

    % ==================================================
    % Step 1: Validate Input Arguments
    % ==================================================
    if nargin < 2
        error('Usage: processExcelData(inputFile, outputFile)');
    end

    % ==================================================
    % Step 2: Read the Excel File as a Cell Array
    % ==================================================
    fprintf('Reading input Excel file: %s\n', inputFile);
    try
        data = readcell(inputFile);
    catch ME
        error(['Failed to read the Excel file. Ensure the file exists and is not open in another program.\n', ...
               'Error: ', ME.message]);
    end

    % Verify that the Excel file has at least one row
    if size(data, 1) < 1
        error('The Excel file must have at least one row to remove.');
    end

    % ==================================================
    % Step 3: Remove the First Row
    % ==================================================
    fprintf('Removing the first row of the Excel data.\n');
    data(1, :) = [];  % Remove the first row

    % Verify that there is data left after removal
    if isempty(data)
        error('No data left after removing the first row.');
    end

    % ==================================================
    % Step 4: Transpose the Data
    % ==================================================
    fprintf('Transposing the data (rows become columns and vice versa).\n');
    data = data';  % Transpose the data

    % ==================================================
    % Step 5: Add 'x' Prefix to the First Row
    % ==================================================
    fprintf('Adding an ''x'' prefix to all entries in the first row.\n');
    for col = 1:size(data, 2)
        entry = data{1, col};
        if ischar(entry) || isstring(entry)
            data{1, col} = ['x' entry];
        elseif isnumeric(entry)
            data{1, col} = ['x' num2str(entry)];
        else
            % Handle other data types (e.g., logical, datetime)
            data{1, col} = ['x' char(string(entry))];
        end
    end

    % ==================================================
    % Step 6: Divide All Numeric Values (Except First Row) by 300
    % ==================================================
    fprintf('Dividing all numeric values (except those in the first row) by 300.\n');
    for row = 2:size(data, 1)  % Start from the second row
        for col = 1:size(data, 2)
            currentValue = data{row, col};
            % Check if the current cell contains a numeric value
            if isnumeric(currentValue) && ~isnan(currentValue)
                data{row, col} = currentValue / 300;
            end
            % If the cell contains text or other types, leave it unchanged
        end
    end

    % ==================================================
    % Step 7: Save the Processed Data to a New Excel File
    % ==================================================
    fprintf('Writing the processed data to the output Excel file: %s\n', outputFile);
    try
        writecell(data, outputFile);
    catch ME
        error(['Failed to write to the Excel file. Ensure you have write permissions and the file is not open.\n', ...
               'Error: ', ME.message]);
    end

    % ==================================================
    % Step 8: Stop Timing and Display Elapsed Time (Optional)
    % ==================================================
    elapsedTime = toc;  % Stop the timer
    fprintf('Processing complete. The output file has been saved as %s.\n', outputFile);
    fprintf('Elapsed time: %.2f seconds.\n', elapsedTime);
end

% ============================================================
% Example Usage:
% ============================================================
% To use this function, simply call it with your input and output filenames.
%
% Example:
%   processExcelData('input.xlsx', 'output.xlsx');
%
% ============================================================

% Uncomment the following line to run the function with example filenames:
% processExcelData('input.xlsx', 'output.xlsx');
