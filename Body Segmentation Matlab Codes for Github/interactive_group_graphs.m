%% Full MATLAB Script for Group-Level Graphs with Menu
clear; clc;

%% --- MENU SELECTION ---
% Menu for graph category (hardcoded binary choices)
choice = menu('Select Graph Category',...
    'Merged ROI Graphs (Statistical Version)',...
    'All ROI Graphs (New Version: Unique Laterality Metrics)',...
    'Group Laterality Graphs',...
    'Per-ROI Laterality Graphs',...
    'All Categories');

% Menu for output mode
outputChoice = menu('Select Output Mode',...
    'Display Only (Do Not Save)',...
    'Save Figures to Folder');

if outputChoice == 2
    saveFigures = true;
else
    saveFigures = false;
end

%% --- COMMON CONFIGURATION ---
subjects = 1:8;
nSubjects = numel(subjects);
dataPathFormat = 'D:\\ML_project\\Variance\\var_excel\\updated_sanitized_allmodels\\subject_%d_variance_partitioning.xlsx';

% Fixed colors
brightGreen  = [0, 0.8, 0];
strongOrange = [1, 0.6, 0];
blackColor   = [0, 0, 0];
redColor     = [1, 0, 0];
greyColor    = [0.5, 0.5, 0.5];

% Custom subject colors (8 distinct colors)
customColors = [1 0 0;
                0 0 1;
                0 1 0;
                1 0.5 0;
                0.5 0 0.5;
                0 1 1;
                1 0 1;
                0.6 0.4 0.2];

% Fixed legend entries (always the same 8 subjects)
subjectLegendEntries = cell(nSubjects,1);
for s = 1:nSubjects
    subjectLegendEntries{s} = sprintf('Subject %d', s);
end

% Uniform y-axis parameters
uniformYTicks_G1 = -0.1:0.02:0.2;
uniformYLim_G1 = [min(uniformYTicks_G1), max(uniformYTicks_G1)];
uniformYTicks_G2 = -0.1:0.02:0.2;
uniformYLim_G2 = [min(uniformYTicks_G2), max(uniformYTicks_G2)];
uniformYTicks_G3 = -0.1:0.02:0.2;
uniformYLim_G3 = [min(uniformYTicks_G3), max(uniformYTicks_G3)];

%% --- OUTPUT FOLDER SETUP (if saving) ---
if saveFigures
    baseOutputFolder = 'GraphOutputs';
    if ~exist(baseOutputFolder, 'dir'), mkdir(baseOutputFolder); end
    mergedStatFolder = fullfile(baseOutputFolder, 'Merged_ROIs_Statistical');
    if ~exist(mergedStatFolder, 'dir'), mkdir(mergedStatFolder); end
    newROI_Folder = fullfile(baseOutputFolder, 'All_ROIs_New_UniqueLaterality');
    if ~exist(newROI_Folder, 'dir'), mkdir(newROI_Folder); end
    groupFolder = fullfile(baseOutputFolder, 'Group_Laterality');
    if ~exist(groupFolder, 'dir'), mkdir(groupFolder); end
    perROIFolder = fullfile(baseOutputFolder, 'PerROI_Laterality');
    if ~exist(perROIFolder, 'dir'), mkdir(perROIFolder); end
end

%% --- Section 1: Merged ROI Graphs (Statistical Version) ---
if choice == 1 || choice == 5
    % Use merged ROIs from first subject (assumes merged ROIs start with a specific letter, e.g., 'm')
    tmpTable = readtable(sprintf(dataPathFormat, subjects(1)));
    mergedIdx = startsWith(tmpTable.ROI, 'm');
    mergedROIs = unique(cellfun(@(x) x(2:end), tmpTable.ROI(mergedIdx), 'UniformOutput', false));
    
    % Preallocate groupData for merged ROIs
    groupData = struct();
    for i = 1:numel(mergedROIs)
        groupData(i).name = mergedROIs{i};
        groupData(i).unique_pose = nan(nSubjects,1);
        groupData(i).unique_seg  = nan(nSubjects,1);
        groupData(i).unique_rpose = nan(nSubjects,1);
        groupData(i).unique_rseg  = nan(nSubjects,1);
        groupData(i).full_R2     = nan(nSubjects,1);
    end
    
    for s = 1:nSubjects
        tbl = readtable(sprintf(dataPathFormat, subjects(s)));
        for i = 1:numel(mergedROIs)
            roiLabel = ['m' mergedROIs{i}];
            idx = strcmp(tbl.ROI, roiLabel);
            if any(idx)
                groupData(i).unique_pose(s) = tbl.unique_pose(find(idx,1));
                groupData(i).unique_seg(s)  = tbl.unique_seg(find(idx,1));
                groupData(i).unique_rpose(s)= tbl.unique_rpose(find(idx,1));
                groupData(i).unique_rseg(s) = tbl.unique_rseg(find(idx,1));
                groupData(i).full_R2(s)     = tbl.Full_R2(find(idx,1));
            end
        end
    end
    
    % Compute paired t-test (Unique Pose vs. Unique Seg) for each merged ROI
    rawP_merged = nan(numel(mergedROIs),1);
    for i = 1:numel(mergedROIs)
        [~, p] = ttest(groupData(i).unique_pose, groupData(i).unique_seg);
        rawP_merged(i) = p;
    end
    % FDR correction across merged ROIs
    N_merged = numel(rawP_merged);
    [sortedP, sortIdx] = sort(rawP_merged);
    q = sortedP .* N_merged ./ (1:N_merged)';
    q = cummin(q(end:-1:1));
    q = q(end:-1:1);
    p_fdr_merged = nan(size(rawP_merged));
    p_fdr_merged(sortIdx) = q;
    clear sortedP sortIdx q;
    
    % Common x-axis for Graph 1 (5 columns)
    xPositions = [1, 1.8, 2.6, 3.4, 4.2];
    x_labels = {'Pose Estimation','Body Segmentation','Random Pose Estimation','Random Body Segmentation','Full Model'};
    xlim_G1 = [min(xPositions)-0.5, max(xPositions)+0.5];
    
    for i = 1:numel(mergedROIs)
        fig = figure('units','normalized','outerposition',[0 0 1 1]);
        hold on;
        for s = 1:nSubjects
            scatter(xPositions(1), groupData(i).unique_pose(s), 50, customColors(s,:), 'filled');
            scatter(xPositions(2), groupData(i).unique_seg(s), 50, customColors(s,:), 'filled');
            scatter(xPositions(3), groupData(i).unique_rpose(s), 50, customColors(s,:), 'filled');
            scatter(xPositions(4), groupData(i).unique_rseg(s), 50, customColors(s,:), 'filled');
            scatter(xPositions(5), groupData(i).full_R2(s),     50, customColors(s,:), 'filled');
        end
        m_vals = [nanmean(groupData(i).unique_pose), nanmean(groupData(i).unique_seg), ...
                  nanmean(groupData(i).unique_rpose), nanmean(groupData(i).unique_rseg), ...
                  nanmean(groupData(i).full_R2)];
        se_vals = [nanstd(groupData(i).unique_pose), nanstd(groupData(i).unique_seg), ...
                   nanstd(groupData(i).unique_rpose), nanstd(groupData(i).unique_rseg), ...
                   nanstd(groupData(i).full_R2)] / sqrt(nSubjects);
        for k = 1:length(xPositions)
            errorbar(xPositions(k), m_vals(k), se_vals(k), 'o', 'Color', blackColor, ...
                'MarkerFaceColor', blackColor, 'MarkerSize',10, 'LineWidth',2, 'CapSize',10);
        end
        set(gca, 'XTick', xPositions, 'XTickLabel', x_labels, 'XLim', xlim_G1, ...
            'YLim', uniformYLim_G1, 'YTick', uniformYTicks_G1);
        ylabel('Variance Explained (R^2)');
        title(sprintf('Merged ROI: %s (Paired t-test: Unique Pose vs. Unique Seg)', mergedROIs{i}));
        
        hSubjects = gobjects(nSubjects,1);
        for s = 1:nSubjects
            hSubjects(s) = scatter(nan, nan, 50, customColors(s,:), 'filled');
        end
        legH = legend(hSubjects, subjectLegendEntries, 'Location','northeast');
        
        if p_fdr_merged(i) < 0.05
            y_sig = max([max(groupData(i).unique_pose), max(groupData(i).unique_seg)]) + 0.01;
            hLine = line([xPositions(1) xPositions(2)], [y_sig y_sig], 'Color', 'k', 'LineWidth', 2);
            set(hLine, 'HandleVisibility','off');
            text(mean([xPositions(1) xPositions(2)]), y_sig+0.005, '*', 'HorizontalAlignment','center', 'FontSize',14, 'Color','k');
        end
        
        annotation('textbox', [0.05, 0.85, 0.4, 0.1], 'String', ...
            sprintf('Paired t-test: uncorrected p = %.4f, FDR p = %.4f (N = %d comparisons)',...
            rawP_merged(i), p_fdr_merged(i), N_merged), ...
            'FitBoxToText', 'on', 'BackgroundColor', 'white', ...
            'EdgeColor', 'black', 'FontSize',10, 'LineStyle', '-', 'Interpreter','none');
        
        legPos = get(legH, 'Position');
        annotation('textbox', [legPos(1), legPos(2)-legPos(4)-0.01, legPos(3), legPos(4)], ...
            'String', 'Significance Legend: *: p < 0.05', ...
            'FitBoxToText', 'on', 'BackgroundColor', 'white', ...
            'EdgeColor', 'black', 'FontSize',10, 'LineStyle', '-', 'Interpreter','none');
        
        hold off;
        if saveFigures
            filename = sprintf('MergedROI_Stat_%s.png', mergedROIs{i});
            saveas(fig, fullfile(mergedStatFolder, filename));
            close(fig);
        else
            pause(1);
        end
    end
end

%% --- Section 2: All ROI Graphs (New Version: Unique Laterality Metrics) ---
if choice == 2 || choice == 5
    % Get all ROI names from the first subject (use ROI names as-is)
    tmpTable = readtable(sprintf(dataPathFormat, subjects(1)));
    allROIs = unique(tmpTable.ROI);
    
    % Define new x-axis for 4 measures in the required order:
    % 1. Unique Pose (Pose Estimation)
    % 2. Unique Segmentation (Body Segmentation)
    % 3. Random Pose Estimation
    % 4. Random Body Segmentation
    xPositions_new = [1, 2, 3, 4];
    x_labels_new = {'Pose Estimation','Body Segmentation','Random Pose Estimation','Random Body Segmentation'};
    xlim_new = [min(xPositions_new)-0.5, max(xPositions_new)+0.5];
    
    for i = 1:numel(allROIs)
        % Preallocate arrays for each measure for this ROI
        poseVals = nan(nSubjects,1);
        segVals  = nan(nSubjects,1);
        rposeVals = nan(nSubjects,1);
        rsegVals  = nan(nSubjects,1);
        
        for s = 1:nSubjects
            tbl = readtable(sprintf(dataPathFormat, subjects(s)));
            idx = strcmp(tbl.ROI, allROIs{i});
            if any(idx)
                poseVals(s)  = tbl.unique_pose(find(idx,1));
                segVals(s)   = tbl.unique_seg(find(idx,1));
                rposeVals(s) = tbl.unique_rpose(find(idx,1));
                rsegVals(s)  = tbl.unique_rseg(find(idx,1));
            end
        end
        
        % Compute group means and standard errors for each measure
        m_vals_new = [nanmean(poseVals), nanmean(segVals), nanmean(rposeVals), nanmean(rsegVals)];
        se_vals_new = [nanstd(poseVals)/sqrt(sum(~isnan(poseVals))), ...
                       nanstd(segVals)/sqrt(sum(~isnan(segVals))), ...
                       nanstd(rposeVals)/sqrt(sum(~isnan(rposeVals))), ...
                       nanstd(rsegVals)/sqrt(sum(~isnan(rsegVals)))];
        
        fig = figure('units','normalized','outerposition',[0 0 1 1]);
        hold on;
        for s = 1:nSubjects
            scatter(xPositions_new(1), poseVals(s), 50, customColors(s,:), 'filled');
            scatter(xPositions_new(2), segVals(s), 50, customColors(s,:), 'filled');
            scatter(xPositions_new(3), rposeVals(s), 50, customColors(s,:), 'filled');
            scatter(xPositions_new(4), rsegVals(s), 50, customColors(s,:), 'filled');
        end
        errorbar(xPositions_new(1), m_vals_new(1), se_vals_new(1), 'o', 'Color', blackColor, 'MarkerFaceColor', blackColor, 'MarkerSize',10, 'LineWidth',2, 'CapSize',10);
        errorbar(xPositions_new(2), m_vals_new(2), se_vals_new(2), 'o', 'Color', blackColor, 'MarkerFaceColor', blackColor, 'MarkerSize',10, 'LineWidth',2, 'CapSize',10);
        errorbar(xPositions_new(3), m_vals_new(3), se_vals_new(3), 'o', 'Color', blackColor, 'MarkerFaceColor', blackColor, 'MarkerSize',10, 'LineWidth',2, 'CapSize',10);
        errorbar(xPositions_new(4), m_vals_new(4), se_vals_new(4), 'o', 'Color', blackColor, 'MarkerFaceColor', blackColor, 'MarkerSize',10, 'LineWidth',2, 'CapSize',10);
        
        yline(0, 'k', 'LineWidth', 1);
        
        set(gca, 'XTick', xPositions_new, 'XTickLabel', x_labels_new, 'XLim', xlim_new, ...
            'YLim', uniformYLim_G1, 'YTick', uniformYTicks_G1);
        ylabel('Unique Variance Explained (R^2)');
        title(sprintf('Unique Variance Explained Across Predictors For Averaged ROI: %s', allROIs{i}), 'Interpreter','none');
        
        hSubjects = gobjects(nSubjects,1);
        for s = 1:nSubjects
            hSubjects(s) = scatter(nan, nan, 50, customColors(s,:), 'filled');
        end
        legend(hSubjects, subjectLegendEntries, 'Location','northeast');
        hold off;
        
        if saveFigures
            filename = sprintf('AllROIs_New_%s.png', allROIs{i});
            saveas(fig, fullfile(newROI_Folder, filename));
            close(fig);
        else
            pause(1);
        end
    end
end

%% --- Section 3: Group Laterality Graphs ---
if choice == 3 || choice == 5
    % (A) Group Laterality: Full Model Only
    leftVals_full = nan(nSubjects,1);
    rightVals_full = nan(nSubjects,1);
    for s = 1:nSubjects
        tbl = readtable(sprintf(dataPathFormat, subjects(s)));
        lIdx = startsWith(tbl.ROI, 'l');
        rIdx = startsWith(tbl.ROI, 'r');
        if any(lIdx)
            leftVals_full(s) = nanmean(tbl.Full_R2(lIdx));
        end
        if any(rIdx)
            rightVals_full(s) = nanmean(tbl.Full_R2(rIdx));
        end
    end
    fig = figure('units','normalized','outerposition',[0 0 1 1]);
    hold on;
    for s = 1:nSubjects
        scatter(1, leftVals_full(s), 60, customColors(s,:), 'filled');
        scatter(2, rightVals_full(s), 60, customColors(s,:), 'filled');
        plot([1 2], [leftVals_full(s) rightVals_full(s)], 'Color', [0.7 0.7 0.7], 'HandleVisibility','off');
    end
    m_left_full = nanmean(leftVals_full);
    m_right_full = nanmean(rightVals_full);
    se_left_full = nanstd(leftVals_full) / sqrt(nSubjects);
    se_right_full = nanstd(rightVals_full) / sqrt(nSubjects);
    errorbar(1, m_left_full, se_left_full, 'o', 'Color', blackColor, 'MarkerFaceColor', blackColor, 'MarkerSize',10, 'LineWidth',2, 'CapSize',10);
    errorbar(2, m_right_full, se_right_full, 'o', 'Color', blackColor, 'MarkerFaceColor', blackColor, 'MarkerSize',10, 'LineWidth',2, 'CapSize',10);
    line([1,2],[m_left_full, m_right_full], 'Color', blackColor, 'LineWidth', 3, 'HandleVisibility','off');
    set(gca, 'XTick', [1 2], 'XTickLabel', {'Left ROIs','Right ROIs'}, 'XLim', [0.8, 2.2], ...
        'YLim', uniformYLim_G2, 'YTick', uniformYTicks_G2);
    ylabel('Mean Full Model R^2');
    title('Group Laterality Test (Full Model Only)', 'Interpreter','none');
    hSubjects = gobjects(nSubjects,1);
    for s = 1:nSubjects
        hSubjects(s) = scatter(nan, nan, 60, customColors(s,:), 'filled');
    end
    legH = legend(hSubjects, subjectLegendEntries, 'Location','northeast');
    [~, p_val_full, ~, stats_full] = ttest(leftVals_full, rightVals_full);
    if p_val_full < 0.05
        significanceText_full = 'statistically significant';
    else
        significanceText_full = 'not statistically significant';
    end
    statReport_full = sprintf('Full Model Laterality Test:\nPaired t-test: t(%d)=%.2f, uncorrected p = %.4f.\nThe difference is %s.', ...
        stats_full.df, stats_full.tstat, p_val_full, significanceText_full);
    annotation('textbox', [0.05, 0.85, 0.4, 0.1], 'String', statReport_full, ...
        'FitBoxToText', 'on', 'BackgroundColor', 'white', 'EdgeColor', 'black', 'FontSize',10, 'LineStyle', '-', 'Interpreter','none');
    legPos = get(legH, 'Position');
    annotation('textbox', [legPos(1), legPos(2)-legPos(4)-0.01, legPos(3), legPos(4)], ...
        'String', 'Significance Legend: *: p < 0.05', ...
        'FitBoxToText', 'on', 'BackgroundColor', 'white', 'EdgeColor', 'black', 'FontSize',10, 'LineStyle', '-', 'Interpreter','none');
    hold off;
    if saveFigures
        saveas(fig, fullfile(groupFolder, 'GroupLaterality_FullModel.png'));
        close(fig);
    else
        pause(1);
    end
    
    % (B) Group Laterality: Unique Pose Only
    leftVals_pose = nan(nSubjects,1);
    rightVals_pose = nan(nSubjects,1);
    for s = 1:nSubjects
        tbl = readtable(sprintf(dataPathFormat, subjects(s)));
        lIdx = startsWith(tbl.ROI, 'l');
        rIdx = startsWith(tbl.ROI, 'r');
        if any(lIdx)
            leftVals_pose(s) = nanmean(tbl.unique_pose(lIdx));
        end
        if any(rIdx)
            rightVals_pose(s) = nanmean(tbl.unique_pose(rIdx));
        end
    end
    fig = figure('units','normalized','outerposition',[0 0 1 1]);
    hold on;
    for s = 1:nSubjects
        scatter(1, leftVals_pose(s), 60, customColors(s,:), 'filled');
        scatter(2, rightVals_pose(s), 60, customColors(s,:), 'filled');
        plot([1 2], [leftVals_pose(s) rightVals_pose(s)], 'Color', [0.7 0.7 0.7], 'HandleVisibility','off');
    end
    m_left_pose = nanmean(leftVals_pose);
    m_right_pose = nanmean(rightVals_pose);
    se_left_pose = nanstd(leftVals_pose) / sqrt(nSubjects);
    se_right_pose = nanstd(rightVals_pose) / sqrt(nSubjects);
    errorbar(1, m_left_pose, se_left_pose, 'o', 'Color', blackColor, 'MarkerFaceColor', blackColor, 'MarkerSize',10, 'LineWidth',2, 'CapSize',10);
    errorbar(2, m_right_pose, se_right_pose, 'o', 'Color', blackColor, 'MarkerFaceColor', blackColor, 'MarkerSize',10, 'LineWidth',2, 'CapSize',10);
    line([1,2],[m_left_pose, m_right_pose], 'Color', blackColor, 'LineWidth', 3, 'HandleVisibility','off');
    set(gca, 'XTick', [1 2], 'XTickLabel', {'Left ROIs','Right ROIs'}, 'XLim', [0.8, 2.2], ...
        'YLim', uniformYLim_G2, 'YTick', uniformYTicks_G2);
    ylabel('Mean Unique Pose Estimation Variance');
    title('Group Laterality Test (Unique Pose Only)', 'Interpreter','none');
    hSubjects = gobjects(nSubjects,1);
    for s = 1:nSubjects
        hSubjects(s) = scatter(nan, nan, 60, customColors(s,:), 'filled');
    end
    legH = legend(hSubjects, subjectLegendEntries, 'Location','northeast');
    [~, p_val_pose, ~, stats_pose] = ttest(leftVals_pose, rightVals_pose);
    if p_val_pose < 0.05
        significanceText_pose = 'statistically significant';
    else
        significanceText_pose = 'not statistically significant';
    end
    statReport_pose = sprintf('Unique Pose Laterality Test:\nPaired t-test: t(%d)=%.2f, uncorrected p = %.4f.\nThe difference is %s.', ...
        stats_pose.df, stats_pose.tstat, p_val_pose, significanceText_pose);
    annotation('textbox', [0.05, 0.85, 0.4, 0.1], 'String', statReport_pose, ...
        'FitBoxToText', 'on', 'BackgroundColor', 'white', 'EdgeColor', 'black', 'FontSize',10, 'LineStyle', '-', 'Interpreter','none');
    legPos = get(legH, 'Position');
    annotation('textbox', [legPos(1), legPos(2)-legPos(4)-0.01, legPos(3), legPos(4)], ...
        'String', 'Significance Legend: *: p < 0.05', ...
        'FitBoxToText', 'on', 'BackgroundColor', 'white', 'EdgeColor', 'black', 'FontSize',10, 'LineStyle', '-', 'Interpreter','none');
    hold off;
    if saveFigures
        saveas(fig, fullfile(groupFolder, 'GroupLaterality_UniquePose.png'));
        close(fig);
    else
        pause(1);
    end
    
    % (C) Group Laterality: Unique Segmentation Only
    leftVals_seg = nan(nSubjects,1);
    rightVals_seg = nan(nSubjects,1);
    for s = 1:nSubjects
        tbl = readtable(sprintf(dataPathFormat, subjects(s)));
        lIdx = startsWith(tbl.ROI, 'l');
        rIdx = startsWith(tbl.ROI, 'r');
        if any(lIdx)
            leftVals_seg(s) = nanmean(tbl.unique_seg(lIdx));
        end
        if any(rIdx)
            rightVals_seg(s) = nanmean(tbl.unique_seg(rIdx));
        end
    end
    fig = figure('units','normalized','outerposition',[0 0 1 1]);
    hold on;
    for s = 1:nSubjects
        scatter(1, leftVals_seg(s), 60, customColors(s,:), 'filled');
        scatter(2, rightVals_seg(s), 60, customColors(s,:), 'filled');
        plot([1 2], [leftVals_seg(s) rightVals_seg(s)], 'Color', [0.7 0.7 0.7], 'HandleVisibility','off');
    end
    m_left_seg = nanmean(leftVals_seg);
    m_right_seg = nanmean(rightVals_seg);
    se_left_seg = nanstd(leftVals_seg) / sqrt(nSubjects);
    se_right_seg = nanstd(rightVals_seg) / sqrt(nSubjects);
    errorbar(1, m_left_seg, se_left_seg, 'o', 'Color', blackColor, 'MarkerFaceColor', blackColor, 'MarkerSize',10, 'LineWidth',2, 'CapSize',10);
    errorbar(2, m_right_seg, se_right_seg, 'o', 'Color', blackColor, 'MarkerFaceColor', blackColor, 'MarkerSize',10, 'LineWidth',2, 'CapSize',10);
    line([1,2],[m_left_seg, m_right_seg], 'Color', blackColor, 'LineWidth', 3, 'HandleVisibility','off');
    set(gca, 'XTick', [1 2], 'XTickLabel', {'Left ROIs','Right ROIs'}, 'XLim', [0.8, 2.2], ...
        'YLim', uniformYLim_G2, 'YTick', uniformYTicks_G2);
    ylabel('Mean Unique Body Segmentation Variance');
    title('Group Laterality Test (Unique Segmentation Only)', 'Interpreter','none');
    hSubjects = gobjects(nSubjects,1);
    for s = 1:nSubjects
        hSubjects(s) = scatter(nan, nan, 60, customColors(s,:), 'filled');
    end
    legH = legend(hSubjects, subjectLegendEntries, 'Location','northeast');
    [~, p_val_seg, ~, stats_seg] = ttest(leftVals_seg, rightVals_seg);
    if p_val_seg < 0.05
        significanceText_seg = 'statistically significant';
    else
        significanceText_seg = 'not statistically significant';
    end
    statReport_seg = sprintf('Unique Segmentation Laterality Test:\nPaired t-test: t(%d)=%.2f, uncorrected p = %.4f.\nThe difference is %s.', ...
        stats_seg.df, stats_seg.tstat, p_val_seg, significanceText_seg);
    annotation('textbox', [0.05, 0.85, 0.4, 0.1], 'String', statReport_seg, ...
        'FitBoxToText', 'on', 'BackgroundColor', 'white', 'EdgeColor', 'black', 'FontSize',10, 'LineStyle', '-', 'Interpreter','none');
    legPos = get(legH, 'Position');
    annotation('textbox', [legPos(1), legPos(2)-legPos(4)-0.01, legPos(3), legPos(4)], ...
        'String', 'Significance Legend: *: p < 0.05', ...
        'FitBoxToText', 'on', 'BackgroundColor', 'white', 'EdgeColor', 'black', 'FontSize',10, 'LineStyle', '-', 'Interpreter','none');
    hold off;
    if saveFigures
        saveas(fig, fullfile(groupFolder, 'GroupLaterality_UniqueSeg.png'));
        close(fig);
    else
        pause(1);
    end
end

%% --- Section 4: Per-ROI Laterality Graphs ---
if choice == 4 || choice == 5
    tbl0 = readtable(sprintf(dataPathFormat, subjects(1)));
    lrIdx = startsWith(tbl0.ROI, 'l') | startsWith(tbl0.ROI, 'r');
    roiNames = cellfun(@(x) x(2:end), tbl0.ROI(lrIdx), 'UniformOutput', false);
    uniqueROIs = unique(roiNames);
    nComparisons = numel(uniqueROIs);

    rawP_full = nan(nComparisons,1);
    roiData_full = struct();
    for i = 1:nComparisons
        roiData_full(i).name = uniqueROIs{i};
        roiData_full(i).left = nan(nSubjects,1);
        roiData_full(i).right = nan(nSubjects,1);
        for s = 1:nSubjects
            tbl = readtable(sprintf(dataPathFormat, subjects(s)));
            idx_left = strcmp(tbl.ROI, ['l' uniqueROIs{i}]);
            idx_right = strcmp(tbl.ROI, ['r' uniqueROIs{i}]);
            if any(idx_left)
                roiData_full(i).left(s) = tbl.Full_R2(find(idx_left,1));
            end
            if any(idx_right)
                roiData_full(i).right(s) = tbl.Full_R2(find(idx_right,1));
            end
        end
        valid = ~isnan(roiData_full(i).left) & ~isnan(roiData_full(i).right);
        if sum(valid) > 1
            [~, p, ~, stats_i] = ttest(roiData_full(i).left(valid), roiData_full(i).right(valid));
            rawP_full(i) = p;
            testStats_full{i} = sprintf('t(%d)=%.2f', stats_i.df, stats_i.tstat);
        else
            rawP_full(i) = NaN;
            testStats_full{i} = 'n/a';
        end
    end
    [sortedP_full, sortIdx] = sort(rawP_full);
    q_full = sortedP_full .* nComparisons ./ (1:nComparisons)';
    q_full = cummin(q_full(end:-1:1));
    q_full = q_full(end:-1:1);
    p_fdr_full = nan(size(rawP_full));
    p_fdr_full(sortIdx) = q_full;
    clear sortedP_full sortIdx q_full

    for i = 1:nComparisons
        fig = figure('units','normalized','outerposition',[0 0 1 1]);
        hold on;
        for s = 1:nSubjects
            if ~isnan(roiData_full(i).left(s)) && ~isnan(roiData_full(i).right(s))
                scatter(1, roiData_full(i).left(s), 60, customColors(s,:), 'filled');
                scatter(2, roiData_full(i).right(s), 60, customColors(s,:), 'filled');
                plot([1 2], [roiData_full(i).left(s) roiData_full(i).right(s)], 'Color', [0.7 0.7 0.7], 'HandleVisibility','off');
            end
        end
        m_left_full_roi = nanmean(roiData_full(i).left);
        m_right_full_roi = nanmean(roiData_full(i).right);
        se_left_full_roi = nanstd(roiData_full(i).left) / sqrt(sum(~isnan(roiData_full(i).left)));
        se_right_full_roi = nanstd(roiData_full(i).right) / sqrt(sum(~isnan(roiData_full(i).right)));
        errorbar(1, m_left_full_roi, se_left_full_roi, 'o', 'Color', blackColor, 'MarkerFaceColor', blackColor, 'MarkerSize',10, 'LineWidth',2, 'CapSize',10);
        errorbar(2, m_right_full_roi, se_right_full_roi, 'o', 'Color', blackColor, 'MarkerFaceColor', blackColor, 'MarkerSize',10, 'LineWidth',2, 'CapSize',10);
        line([1 2], [m_left_full_roi, m_right_full_roi], 'Color', blackColor, 'LineWidth', 3, 'HandleVisibility','off');
        set(gca, 'XTick', [1 2], 'XTickLabel', {'Left','Right'}, 'XLim', [0.8, 2.2], ...
            'YLim', uniformYLim_G3, 'YTick', uniformYTicks_G3);
        ylabel('Full Model R^2');
        title(sprintf('Per-ROI Laterality (Full Model) for ROI: %s', uniqueROIs{i}), 'Interpreter','none');
        hSubjects = gobjects(nSubjects,1);
        for s = 1:nSubjects
            hSubjects(s) = scatter(nan, nan, 60, customColors(s,:), 'filled');
        end
        legH = legend(hSubjects, subjectLegendEntries, 'Location','northeast');
        statReport_full_roi = sprintf('For ROI %s, paired t-test (Full Model): %s, uncorrected p = %.4f, FDR p = %.4f (N = %d)', ...
            uniqueROIs{i}, testStats_full{i}, rawP_full(i), p_fdr_full(i), nComparisons);
        annotation('textbox', [0.05, 0.85, 0.4, 0.1], 'String', statReport_full_roi, ...
            'FitBoxToText', 'on', 'BackgroundColor', 'white', 'EdgeColor', 'black', 'FontSize',10, 'LineStyle', '-', 'Interpreter','none');
        legPos = get(legH, 'Position');
        annotation('textbox', [legPos(1), legPos(2)-legPos(4)-0.01, legPos(3), legPos(4)], ...
            'String', 'Significance Legend: *: p < 0.05', ...
            'FitBoxToText', 'on', 'BackgroundColor', 'white', 'EdgeColor', 'black', 'FontSize',10, 'LineStyle', '-', 'Interpreter','none');
        hold off;
        filename = sprintf('PerROI_FullModel_%s.png', uniqueROIs{i});
        saveas(fig, fullfile(perROIFolder, filename));
        close(fig);
    end
end

%% --- End of Script ---
