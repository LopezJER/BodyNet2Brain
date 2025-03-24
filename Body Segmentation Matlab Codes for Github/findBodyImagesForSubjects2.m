function findBodyImagesForSubjects2(excelFilePath, bodyImageFilePath, outputFilePath)
    % Load the data from the Excel file
    data = readtable(excelFilePath);

    % Load the body image IDs from the text file
    bodyImageIDs = readmatrix(bodyImageFilePath);

    % Initialize an empty table for output
    outputData = table;

    % Add the body image IDs as the first column, filtered by the order in cocoID
    cocoIDsInOrder = intersect(data.cocoId, bodyImageIDs, 'stable');
    outputData.cocoId = cocoIDsInOrder;

    % Process each subject (1 to 8)
    for subjectNumber = 1:8
        % Generate the column names for the subject
        subjectColumn = sprintf('subject%d', subjectNumber);
        rep0Column = sprintf('subject%d_rep0', subjectNumber);
        rep1Column = sprintf('subject%d_rep1', subjectNumber);
        rep2Column = sprintf('subject%d_rep2', subjectNumber);

        % Check if the required columns exist in the data
        requiredColumns = {subjectColumn, rep0Column, rep1Column, rep2Column};
        if ~all(ismember(requiredColumns, data.Properties.VariableNames))
            error('The specified columns for subject %d do not exist in the dataset.', subjectNumber);
        end

        % Initialize empty vectors for storing trial IDs shown to the subject
        trial_rep0 = NaN(height(outputData), 1);
        trial_rep1 = NaN(height(outputData), 1);
        trial_rep2 = NaN(height(outputData), 1);

        % Iterate through the filtered cocoIDs and find matches
        for i = 1:height(outputData)
            imageID = outputData.cocoId(i); % Get the cocoId in the desired order

            % Find the row in the original data that matches this imageID
            matchingRow = find(data.cocoId == imageID & data.(subjectColumn) == 1);

            if ~isempty(matchingRow) % Ensure there's a valid match
                % Store the trial IDs from rep0, rep1, and rep2
                trial_rep0(i) = data.(rep0Column)(matchingRow);
                trial_rep1(i) = data.(rep1Column)(matchingRow);
                trial_rep2(i) = data.(rep2Column)(matchingRow);
            end
        end

        % Add columns for the trial IDs to the output table
        outputData.(rep0Column) = trial_rep0;
        outputData.(rep1Column) = trial_rep1;
        outputData.(rep2Column) = trial_rep2;
    end

    % Write the updated table back to the spreadsheet
    writetable(outputData, outputFilePath);

    % Inform the user of the results
    fprintf('Processing complete for Subjects 1-8. Results saved to %s\n', outputFilePath);
end
