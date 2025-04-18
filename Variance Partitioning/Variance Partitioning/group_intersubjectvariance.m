%% Clear Workspace, Command Window, and Figures
clc; clear; close all;

%% Part 1: Data Reading and Filtering
% Folder path and file selection (only files with "variance_partitioning" in their name)
inputFolder = 'D:\ML_project\Variance\var_excel\updated_sanitized_allmodels';  
files = dir(fullfile(inputFolder, '*.xlsx'));

allTables = {};  % To store each file's table

for k = 1:length(files)
    fileName = files(k).name;
    
    % Process only files that contain 'variance_partitioning'
    if ~contains(lower(fileName), 'variance_partitioning')
        fprintf('Skipping file (does not match pattern): %s\n', fileName);
        continue;
    end
    
    fullPath = fullfile(inputFolder, fileName);
    
    % Extract subject ID (assumes the first number in the filename is the subject ID)
    token = regexp(fileName, '\d+', 'match');
    if isempty(token)
        warning('No subject ID found in file: %s. Skipping.', fileName);
        continue;
    end
    subjID = str2double(token{1});
    
    % Read the Excel file.
    T = readtable(fullPath);
    
    % Check for 'ROI' column.
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
    error('No valid Excel files were processed.');
end

% Combine tables and ensure ROI is categorical.
combinedTable = vertcat(allTables{:});
if ~iscategorical(combinedTable.ROI)
    combinedTable.ROI = categorical(combinedTable.ROI);
end

% Filter for ROIs that start with 'm' (case-insensitive)
roiStrings = cellstr(combinedTable.ROI);
mask = startsWith(roiStrings, 'm', 'IgnoreCase', true);
filteredTable = combinedTable(mask, :);

% Restrict analysis to only the three predictors of interest.
predictorNames = {'Full_R2', 'unique_seg', 'unique_pose'};
% Define corresponding display labels.
predictorDisplay = {'Full Model R^2', 'Unique Segmentation Variance', 'Unique Pose Estimation Variance'};

%% Part 2: Variability Metrics Calculation (Justifying Averaging)
% Since there is only one measurement per subject for each ROI, a formal ANOVA of 
% inter-subject variability isn’t possible. Instead, we calculate descriptive 
% metrics such as standard deviation (SD), the coefficient of variation (CV), and
% the intraclass correlation coefficient (ICC) computed across ROIs (for each predictor).
%
% First, compute summary statistics per ROI:

[grp, roiList] = findgroups(filteredTable.ROI);
numROIs = length(roiList);
summaryTable = table;
summaryTable.ROI = roiList;
for p = 1:length(predictorNames)
    predName = predictorNames{p};
    % Mean and standard deviation across subjects per ROI
    meanVals = splitapply(@(x) mean(x, 'omitnan'), filteredTable.(predName), grp);
    stdVals  = splitapply(@(x) std(x, 'omitnan'), filteredTable.(predName), grp);
    cvVals   = stdVals ./ meanVals;  % may be negative if mean is negative
    
    summaryTable.(sprintf('%s_mean', predName)) = meanVals;
    summaryTable.(sprintf('%s_std', predName))  = stdVals;
    summaryTable.(sprintf('%s_CV', predName))   = cvVals;
end
disp('Summary Table (per ROI, averaged across subjects):');
disp(summaryTable);

% Next, compute the ICC for each predictor across ROIs.
uniqueSubjects = unique(filteredTable.SubjectID);
uniqueROIs = unique(filteredTable.ROI);
numSubj = length(uniqueSubjects);
numROI = length(uniqueROIs);

% Initialize structure to hold ICC values for each predictor.
iccValues = struct();
for p = 1:length(predictorNames)
    predictor = predictorNames{p};
    % Preallocate a matrix for the predictor data
    dataMatrix = NaN(numSubj, numROI);
    for i = 1:numSubj
        for j = 1:numROI
            subj = uniqueSubjects(i);
            roi  = uniqueROIs(j);
            idx = (filteredTable.SubjectID == subj) & (filteredTable.ROI == roi);
            if any(idx)
                dataMatrix(i,j) = filteredTable.(predictor)(idx);
            end
        end
    end
    % Compute ICC for the dataMatrix using the custom function below.
    icc_val = compute_icc(dataMatrix);
    iccValues.(predictor) = icc_val;
    fprintf('ICC for %s: %.2f\n', predictorDisplay{p}, icc_val);
end

%% Verbal Conclusions per ROI and Predictor
fprintf('\nVerbal Conclusions (per ROI and Predictor):\n');
cvThreshold = 0.20;
sdThresholdFraction = 0.10;  % SD must be less than 10% of absolute mean to be considered low variability

for i = 1:height(summaryTable)
    roiName = string(summaryTable.ROI(i));
    fprintf('ROI %s:\n', roiName);
    for p = 1:length(predictorNames)
        predName = predictorNames{p};
        dispLabel = predictorDisplay{p};
        mean_val = summaryTable.(sprintf('%s_mean', predName))(i);
        std_val = summaryTable.(sprintf('%s_std', predName))(i);
        cv_val = summaryTable.(sprintf('%s_CV', predName))(i);
        
        abs_cv = abs(cv_val);
        if abs(mean_val) > 0
            sd_fraction = std_val / abs(mean_val);
        else
            sd_fraction = NaN;
        end
        
        conclusion = sprintf('  For %s: Mean = %.2f, SD = %.2f, CV = %.2f (|CV| = %.2f, SD/|Mean| = %.2f).', ...
            dispLabel, mean_val, std_val, cv_val, abs_cv, sd_fraction);
        
        if abs_cv < cvThreshold && sd_fraction < sdThresholdFraction
            conclusion = [conclusion, ' Variability is low.'];
        else
            conclusion = [conclusion, ' Variability is high.'];
        end
        
        overallICC = iccValues.(predName);
        if overallICC > 0.70
            conclusion = [conclusion, sprintf(' Overall ICC = %.2f indicates high reliability across ROIs.', overallICC)];
        else
            conclusion = [conclusion, sprintf(' Overall ICC = %.2f indicates moderate to low reliability across ROIs.', overallICC)];
        end
        
        if (abs_cv < cvThreshold && sd_fraction < sdThresholdFraction && overallICC > 0.70)
            conclusion = [conclusion, ' Averaging across subjects is justified.'];
        else
            conclusion = [conclusion, ' Averaging across subjects may obscure individual differences.'];
        end
        
        fprintf('%s\n', conclusion);
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Part 3: Scatter Plots for Each Predictor (Variance Plots)
% This section produces scatter plots displaying individual subject data points per ROI
% along with the average and error bars (±1 SD) for each predictor.

roiCategories = unique(filteredTable.ROI);
nRois = length(roiCategories);
subjIDs = unique(filteredTable.SubjectID);
nSubj = length(subjIDs);

% Define a custom set of 8 distinct colors.
customColorsPlot = [1 0 0;       % red
                    0 0 1;       % blue
                    0 1 0;       % green
                    1 0.5 0;     % orange
                    0.5 0 0.5;   % purple
                    0 1 1;       % cyan
                    1 0 1;       % magenta
                    0.6 0.4 0.2];% brown

if nSubj <= 8
    colors = customColorsPlot(1:nSubj, :);
else
    colors = lines(nSubj); % fallback if more than 8 subjects
end

for p = 1:length(predictorNames)
    predictor = predictorNames{p};
    dispLabel = predictorDisplay{p};
    figure('Units','normalized','OuterPosition',[0 0 1 1]);
    hold on;
    
    for i = 1:nRois
        currentROI = roiCategories(i);
        roiIdx = filteredTable.ROI == currentROI;
        dataRoi = filteredTable.(predictor)(roiIdx);
        subjRoi = filteredTable.SubjectID(roiIdx);
        
        for s = 1:nSubj
            subjIdx = subjRoi == subjIDs(s);
            if any(subjIdx)
                values = dataRoi(subjIdx);
                scatter(i * ones(size(values)), values, 50, colors(s,:), 'filled');
            end
        end
        
        if ~isempty(dataRoi)
            mVal = mean(dataRoi, 'omitnan');
            sVal = std(dataRoi, 'omitnan');
            errorbar(i, mVal, sVal, 'k', 'LineWidth', 2, 'CapSize', 10);
            scatter(i, mVal, 120, 'k', 'filled');  % Bold average point.
        end
    end
    
    hLegend = gobjects(nSubj,1);
    for s = 1:nSubj
        hLegend(s) = scatter(NaN, NaN, 50, colors(s,:), 'filled', 'DisplayName', sprintf('Subject %d', subjIDs(s)));
    end
    legend(hLegend, 'Location', 'bestoutside');
    
    xlim([0 nRois+1]);
    set(gca, 'XTick', 1:nRois, 'XTickLabel', cellstr(roiCategories), 'XTickLabelRotation', 45, 'TickLabelInterpreter', 'none');
    
    % Set constant y-axis limits and ticks.
    ylim([-0.1, 0.12]);
    set(gca, 'YTick', -0.1:0.02:0.12);
    
    title(sprintf('Distribution of %s Across ROIs (Averaged Across Subjects)', dispLabel), 'FontSize', 14, 'FontWeight', 'bold', 'Interpreter', 'none');
    xlabel('ROI', 'Interpreter', 'none');
    ylabel(dispLabel, 'Interpreter', 'none');
    grid on;
    hold off;
    
    % Save the current figure as a PNG with a timestamp.
    timestamp = datestr(now, 'yyyymmdd_HHMMSS');
    safeLabel = regexprep(dispLabel, '\s+', '_');  % Replace spaces with underscores
    filename = fullfile(inputFolder, sprintf('Scatter_%s_%s.png', safeLabel, timestamp));
    saveas(gcf, filename);
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Supporting Function: Compute ICC (ICC(1,1) one-way random effects model)
function icc = compute_icc(X)
    % X is a matrix with rows = subjects and columns = ROIs (or repeated measurements).
    % Remove any rows with NaNs (requires complete data for ICC computation)
    X = X(~any(isnan(X),2),:);
    [n, k] = size(X);
    if n < 2 || k < 2
        warning('Not enough data to compute ICC.');
        icc = NaN;
        return;
    end

    % Compute the subject means and the grand mean.
    subjMeans = mean(X,2);
    grandMean = mean(X(:));
    
    % Between-subjects sum of squares.
    SS_between = k * sum((subjMeans - grandMean).^2);
    df_between = n - 1;
    BMS = SS_between / df_between;
    
    % Within-subjects sum of squares.
    SS_within = sum(sum((X - subjMeans).^2));
    df_within = n*(k-1);
    WMS = SS_within / df_within;
    
    % ICC(1,1)
    icc = (BMS - WMS) / (BMS + (k-1)*WMS);
end
