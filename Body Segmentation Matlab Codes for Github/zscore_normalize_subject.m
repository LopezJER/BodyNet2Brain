function zscore_normalize_subject()

    % === HARDCODED PATHS ===
    inputFolder = 'D:\ML_project\RDM_results\final\subj8';         % <-- CHANGE THIS
    outputFolder = 'D:\ML_project\RDM_results\zscore\subj8';       % <-- CHANGE THIS

    % Create output folder if it doesn't exist
    if ~exist(outputFolder, 'dir')
        mkdir(outputFolder);
        fprintf('Created output folder: %s\n', outputFolder);
    end

    roiFolders = dir(inputFolder);
    roiFolders = roiFolders([roiFolders.isdir] & ~startsWith({roiFolders.name}, '.'));

    allData = [];

    % === First pass: collect all numeric data across ROIs ===
    for i = 1:length(roiFolders)
        roiName = roiFolders(i).name;
        roiPath = fullfile(inputFolder, roiName);
        filePath = fullfile(roiPath, [roiName, '.xlsx']);

        if exist(filePath, 'file')
            try
                T = readtable(filePath);
                if size(T,1) > 0
                    numericData = T{:, :};
                    allData = [allData; numericData]; % stack rows
                else
                    fprintf('File %s has no data rows.\n', filePath);
                end
            catch
                fprintf('Error reading file: %s (is it really a table of numbers?)\n', filePath);
            end
        else
            fprintf('Skipping ROI "%s": no Excel file named "%s.xlsx" found.\n', roiName, roiName);
        end
    end

    % === Compute global stats ===
    dataMean = mean(allData(:), 'omitnan');
    dataStd  = std(allData(:), 'omitnan');

    fprintf('\nGlobal mean: %.4f\n', dataMean);
    fprintf('Global std : %.4f\n\n', dataStd);

    % === Second pass: normalize and save with headers ===
    for i = 1:length(roiFolders)
        roiName = roiFolders(i).name;
        roiPath = fullfile(inputFolder, roiName);
        filePath = fullfile(roiPath, [roiName, '.xlsx']);

        if exist(filePath, 'file')
            try
                T = readtable(filePath);

                numericData = T{:, :};
                zdata = (numericData - dataMean) / dataStd;

                % Build new table with same headers
                T_z = array2table(zdata, 'VariableNames', T.Properties.VariableNames);

                % Save
                outputName = [roiName, '_zscore.xlsx'];
                outputPath = fullfile(outputFolder, outputName);
                writetable(T_z, outputPath);

                fprintf('Z-scored file saved: %s\n', outputPath);

            catch ME
                fprintf('Error processing %s: %s\n', filePath, ME.message);
            end
        end
    end

    disp('All done! Z-score normalization complete.');
end
