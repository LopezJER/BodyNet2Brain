% Before and after taking out low level features - comparison of averaged
% excels

%% Configuration
% Hardcoded input file paths:
beforeFile = 'D:\ML_project\Variance\var_excel\unsanitized and unmerged\Aggregated_ROI_Results_uu_only_always.xlsx';  % "Before" Excel file
afterFile  = 'D:\ML_project\Variance\var_excel\updated_sanitized_allmodels\Aggregated_ROI_Results.xlsx';              % "After" Excel file

% Hardcoded output folder (will be created if it doesn't exist)
outputFolder = 'D:\ML_project\Variance\var_excel';
if ~exist(outputFolder, 'dir')
    mkdir(outputFolder);
end

% Columns to compare (exact names in both Excel files)
columnsToPlot = {'Avg_Full_R2', 'Avg_unique_pose', 'Avg_unique_seg'};

% Create a minimal timestamp for output filenames (format: yyyymmdd_HHMMSS)
timestamp = datestr(now, 'yyyymmdd_HHMMSS');

%% Read the Data
% We assume both files share the same structure and row order.
beforeTable = readtable(beforeFile);
afterTable  = readtable(afterFile);

% Assume there is a column named ROI in both.
roiList = beforeTable.ROI;
roiListStr = cellstr(roiList);
nROIs = length(roiListStr);

%% Create Color and Marker Mapping for ROIs
% Each ROI name is assumed to start with 'l' or 'r'. ROIs sharing the same base name
% (i.e. the characters after the first letter) will share the same color.
baseNames = cellfun(@(s) s(2:end), roiListStr, 'UniformOutput', false);
uniqueBaseNames = unique(baseNames);
% Use a colormap with at least 10 distinct colors.
colorMat = lines(max(10, length(uniqueBaseNames)));
colorMap = containers.Map;
for i = 1:length(uniqueBaseNames)
    colorMap(uniqueBaseNames{i}) = colorMat(i,:);
end

%% Reorder Data by Base Name (so similar ROIs appear together)
% Sort the ROI names by their base names (ignoring case)
[~, sortIdx] = sort(lower(cellfun(@(s) s(2:end), roiListStr, 'UniformOutput', false)));
roiListStrSorted = roiListStr(sortIdx);

%% Loop over Each Column to Compare and Plot Graphs
for k = 1:length(columnsToPlot)
    colName = columnsToPlot{k};  % e.g., 'Avg_Full_R2'
    
    % Extract the before and after values from the respective tables.
    rawBefore = beforeTable.(colName);
    rawAfter  = afterTable.(colName);
    
    % Reorder the values based on the sorted order.
    beforeValuesSorted = rawBefore(sortIdx);
    afterValuesSorted  = rawAfter(sortIdx);
    roiListSorted = roiListStrSorted;  % Sorted ROI names
    
    % Define x positions: before at 1, after at 1.2.
    xBefore = 1;
    xAfter  = 1.2;
    
    % Compute overall means and standard deviations (ignoring NaNs).
    meanBefore = mean(beforeValuesSorted, 'omitnan');
    stdBefore  = std(beforeValuesSorted, 'omitnan');
    meanAfter  = mean(afterValuesSorted, 'omitnan');
    stdAfter   = std(afterValuesSorted, 'omitnan');
    
    % Compute differences and Cohen's d (absolute value for two-tailed test)
    diffValues = beforeValuesSorted - afterValuesSorted;
    d_cohen = abs(mean(diffValues, 'omitnan') / std(diffValues, 'omitnan'));
    
    % Interpret Cohen's d using common benchmarks.
    if d_cohen < 0.2
        effectInterp = 'Negligible';
    elseif d_cohen < 0.5
        effectInterp = 'Small';
    elseif d_cohen < 0.8
        effectInterp = 'Medium';
    else
        effectInterp = 'Large';
    end
    
    % Perform a two-tailed paired t-test.
    [~, p, ~, stats] = ttest(beforeValuesSorted, afterValuesSorted);
    % Bonferroni correction for 3 comparisons: multiply p by 3 (capped at 1).
    p_adj = min(p * 3, 1);
    if p_adj < 0.05
        sigStr = 'Significant';
    else
        sigStr = 'Not Significant';
    end
    
    % Create a large figure.
    fig = figure('Position', [100, 100, 1500, 1000]);
    hold on;
    
    % Set x-axis limits to leave extra margin on left and right.
    xlim([0.8, 1.4]);
    
    % Preallocate arrays for legend handles and entries.
    handles = gobjects(nROIs, 1);
    legendEntries = cell(nROIs, 1);
    
    % Loop over each ROI (in sorted order).
    for i = 1:nROIs
        roiStr = roiListSorted{i};  % e.g., 'rEBA' or 'lEBA'
        % Determine marker style based on the first letter.
        if lower(roiStr(1)) == 'l'
            marker = 's';      % filled square for left
            markerFace = 'auto';
        else
            marker = 'o';      % open circle for right
            markerFace = 'none';
        end
        % Use the base name (everything after the first character) for color.
        baseName = roiStr(2:end);
        col = colorMap(baseName);
        
        % Plot the ROI's paired data.
        handles(i) = plot([xBefore, xAfter], [beforeValuesSorted(i), afterValuesSorted(i)], '-', ...
            'Color', col, 'LineWidth', 1.5, 'Marker', marker, ...
            'MarkerFaceColor', markerFace, 'MarkerSize', 8);
        
        % For the legend, keep the full ROI name (including prefix).
        legendEntries{i} = roiStr;
    end
    
    % Plot overall average error bars and capture their handles.
    hErrorBefore = errorbar(xBefore, meanBefore, stdBefore, 'ok', 'LineWidth', 3, 'MarkerSize', 10, 'CapSize', 10);
    hErrorAfter  = errorbar(xAfter, meanAfter, stdAfter, 'ok', 'LineWidth', 3, 'MarkerSize', 10, 'CapSize', 10);
    
    % Plot the overall mean line as a thick, dashed black line.
    hMeanLine = plot([xBefore, xAfter], [meanBefore, meanAfter], '--k', 'LineWidth', 3);
    
    % Append overall mean info to the legend.
    %legendEntries{end+1} = sprintf('Before Mean and Std: %.4f, %.4f', meanBefore, stdBefore);
    %legendEntries{end+1} = sprintf('After Mean and Std: %.4f, %.4f', meanAfter, stdAfter);
    
    % Set x-axis tick labels.
    set(gca, 'XTick', [xBefore, xAfter], 'XTickLabel', {'Before','After'});
    
    % Set the y-axis limits and constant tick marks across graphs.
    ylim([-0.02, 0.06]);
    set(gca, 'YTick', -0.02:0.01:0.06);
    
    % Create the legend with all handles and corresponding entries.
    lgd = legend([handles; hErrorBefore; hErrorAfter; hMeanLine], legendEntries, 'Location', 'northeast', 'FontSize', 8);
    lgd.Interpreter = 'none';
    drawnow;  % Update legend position
    
    % Get the legend's position in normalized units.
    lgdPos = lgd.Position;  % [x, y, width, height]
    % Place the annotation box immediately below the legend.
    annX = lgdPos(1);
    annW = lgdPos(3);
    annH = lgdPos(4);  % same width as legend; height equal to legend's height
    annY = lgdPos(2) - annH - 0.01;  % just below the legend
    
    % Prepare an annotation string including the t-test results and effect size.
    annotationStr = sprintf(['Two-Tailed Paired t-Test (p x 3 for Bonferroni correction):\n' ...
                             'Uncorrected p = %.4f, Corrected p = %.4f\nResult: %s\n' ...
                             't(%d) = %.2f\nEffect Size (Cohen''s d) = %.2f (%s)'], ...
                             p, p_adj, sigStr, stats.df, stats.tstat, d_cohen, effectInterp);
    annotation('textbox', [annX, annY, annW, annH], 'String', annotationStr, ...
               'EdgeColor', 'black', 'FontSize', 10, 'BackgroundColor', 'w', 'Interpreter', 'none');
    
    % Change y-axis label based on the column name.
    switch colName
        case 'Avg_Full_R2'
            ylabel('Average Full model R^2 values', 'FontSize', 14, 'Interpreter', 'none');
        case 'Avg_unique_pose'
            ylabel('Average Unique Pose Estimation Variance', 'FontSize', 14, 'Interpreter', 'none');
        case 'Avg_unique_seg'
            ylabel('Average Unique Body Segmentation Variance', 'FontSize', 14, 'Interpreter', 'none');
        otherwise
            ylabel([colName ' Values'], 'FontSize', 14, 'Interpreter', 'none');
    end
    
    % Set main title based on the column name.
    switch colName
        case 'Avg_Full_R2'
            mainTitle = 'Average Full model R^2 values Comparison: Before and After Controlling For Low Level Features';
        case 'Avg_unique_pose'
            mainTitle = 'Average Unique Pose Estimation Variance Comparison: Before and After Controlling For Low Level Features';
        case 'Avg_unique_seg'
            mainTitle = 'Average Unique Body Segmentation Variance Comparison: Before and After Controlling For Low Level Features';
        otherwise
            mainTitle = sprintf('%s Comparison: Before and After Controlling For Low Level Features', colName);
    end
    title(mainTitle, 'FontSize', 14, 'Interpreter', 'none');
    
    xlabel('Condition', 'FontSize', 14, 'Interpreter', 'none');
    
    hold off;
    
    % Append the timestamp to the output filename.
    outFile = fullfile(outputFolder, sprintf('%s_Comparison_%s.png', colName, timestamp));
    saveas(fig, outFile);
    close(fig);
end
