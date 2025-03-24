% Clear workspace and command window
clear;
clc;

% Configuration: Base directory for input files
BASE_INPUT_DIR = 'D:\ML_project\RDM_results\final\subj8';  % Change this to your input folder

fprintf('=== Vectorizing RDMs from Excel Files ===\n');

% Recursively find all .xlsx files in the input directory
xlsxFiles = dir(fullfile(BASE_INPUT_DIR, '**', '*.xlsx'));
numFiles = length(xlsxFiles);

fprintf('Found %d Excel files before filtering.\n', numFiles);

% Process each Excel file that has a base name shorter than 10 characters
for i = 1:numFiles
    inputFile = fullfile(xlsxFiles(i).folder, xlsxFiles(i).name);
    [filePath, baseName, ~] = fileparts(inputFile);
    
    if length(baseName) < 10
        fprintf('\n[INFO] Processing Excel file: %s\n', inputFile);
        
        % Load data from Excel, treating the file as not having variable names.
        % This way, the first row (image names) is read as data.
        data = readtable(inputFile, 'ReadVariableNames', false);
        
        % Remove the first row (which contains image names)
        if size(data,1) < 2
            fprintf('[SKIP] Not enough rows after removing header. Skipping %s\n', inputFile);
            continue;
        end
        vectors = table2array(data(2:end, :));
        
        % Now, columns represent images.
        numImages = size(vectors, 2);
        fprintf('[INFO] Number of images (columns): %d\n', numImages);
        
        % (Optional) You could check if the number of images is what you expect.
        % For example, if you expect 613 images, you might want to error if not.
        
        % Compute RDM using 1 - Pearson's Correlation
        % Initialize RDM matrix (numImages x numImages)
        RDM = zeros(numImages);
        fprintf('[INFO] Computing RDM...\n');
        for m = 1:numImages
            for n = m+1:numImages
                r = corr(vectors(:, m), vectors(:, n));
                if isnan(r)
                    r = 0;  % Replace any NaN correlation with 0
                end
                RDM(m, n) = 1 - r;
            end
        end
        RDM = RDM + RDM';  % Make RDM symmetric
        
        % Vectorize the RDM: extract the upper triangle (excluding the diagonal)
        upper_triangle = triu(true(size(RDM)), 1);
        rdm_vec = RDM(upper_triangle);
        vec_length = length(rdm_vec);
        fprintf('[DEBUG] Vectorized RDM Length: %d\n', vec_length);
        
        % Before saving, delete any existing file with the same name
        outputFile = fullfile(filePath, [baseName, '_rdm_vec.xlsx']);
        if isfile(outputFile)
            fprintf('[INFO] Deleting existing file: %s\n', outputFile);
            delete(outputFile);
        end
        
        % Save the vectorized RDM to the same folder with _rdm_vec appended
        fprintf('[INFO] Saving Vectorized RDM to: %s\n', outputFile);
        writematrix(rdm_vec, outputFile);
        
        fprintf('[SUCCESS] Vectorized RDM saved successfully for: %s\n', inputFile);
    else
        fprintf('[SKIP] Skipping %s (Name too long: %d characters)\n', inputFile, length(baseName));
    end
end

fprintf('\n=== All Excel files processed successfully. ===\n');
