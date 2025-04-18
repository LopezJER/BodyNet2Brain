%% Full MATLAB Script for Group-Level Graphs (Full Screen Figures)
clear; clc;

%% Configuration
subjects = 1:8;
nSubjects = numel(subjects);
dataPathFormat = 'D:\\ML_project\\Variance\\var_excel\\allmodels_sanitized\\subject_%d_variance_partitioning.xlsx';

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

% Fixed legend entries for subjects (always the same)
subjectLegendEntries = cell(nSubjects,1);
for s = 1:nSubjects
    subjectLegendEntries{s} = sprintf('Subject %d', s);
end

% Define uniform y-ticks and limits for graphs
uniformYTicks_G1 = -0.1:0.02:0.2;
uniformYLim_G1 = [min(uniformYTicks_G1), max(uniformYTicks_G1)];
uniformYTicks_G2 = -0.1:0.02:0.2;
uniformYLim_G2 = [min(uniformYTicks_G2), max(uniformYTicks_G2)];
uniformYTicks_G3 = -0.1:0.02:0.2;
uniformYLim_G3 = [min(uniformYTicks_G3), max(uniformYTicks_G3)];

%% Create Base Output Folder and Subfolders
baseOutputFolder = 'GraphOutputs';
if ~exist(baseOutputFolder, 'dir')
    mkdir(baseOutputFolder);
end
% For merged ROI comparisons (Graph 1)
mergedFolder = fullfile(baseOutputFolder, 'Merged_ROIs');
if ~exist(mergedFolder, 'dir')
    mkdir(mergedFolder);
end
% For group laterality tests:
groupFullFolder = fullfile(baseOutputFolder, 'Group_Laterality_FullModel');
if ~exist(groupFullFolder, 'dir')
    mkdir(groupFullFolder);
end
groupSegFolder = fullfile(baseOutputFolder, 'Group_Laterality_UniqueSeg');
if ~exist(groupSegFolder, 'dir')
    mkdir(groupSegFolder);
end
groupPoseFolder = fullfile(baseOutputFolder, 'Group_Laterality_UniquePose');
if ~exist(groupPoseFolder, 'dir')
    mkdir(groupPoseFolder);
end
% For per-ROI laterality tests:
perROIFullFolder = fullfile(baseOutputFolder, 'PerROI_Laterality_FullModel');
if ~exist(perROIFullFolder, 'dir')
    mkdir(perROIFullFolder);
end
perROISegFolder = fullfile(baseOutputFolder, 'PerROI_Laterality_UniqueSeg');
if ~exist(perROISegFolder, 'dir')
    mkdir(perROISegFolder);
end
perROIPoseFolder = fullfile(baseOutputFolder, 'PerROI_Laterality_UniquePose');
if ~exist(perROIPoseFolder, 'dir')
    mkdir(perROIPoseFolder);
end

%% --- PRE-COMPUTE MERGED ROI DATA (for Graph 1) ---
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

% For Graph 1: perform a paired t-test between Unique Pose and Unique Seg for each merged ROI.
rawP_merged = nan(numel(mergedROIs),1);
for i = 1:numel(mergedROIs)
    [~, p] = ttest(groupData(i).unique_pose, groupData(i).unique_seg);
    rawP_merged(i) = p;
end
% FDR correction across merged ROIs (Benjamini-Hochberg)
N_merged = numel(rawP_merged);
[sortedP, sortIdx] = sort(rawP_merged);
q = sortedP .* N_merged ./ (1:N_merged)';
q = cummin(q(end:-1:1));
q = q(end:-1:1);
p_fdr_merged = nan(size(rawP_merged));
p_fdr_merged(sortIdx) = q;
clear sortedP sortIdx q;

% Common x-axis parameters for Graph 1
xPositions = [1, 1.8, 2.6, 3.4, 4.2];
x_labels = {'Pose Estimation','Body Segmentation','Random Pose Estimation','Random Body Segmentation','Full Model'};
xlim_G1 = [min(xPositions)-0.5, max(xPositions)+0.5];

%% Graph 1: Merged ROI Scatter Plot (Unique Pose vs. Unique Seg)
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
    title(sprintf('Merged ROI: %s (Unique Pose vs. Unique Seg, t-test)', mergedROIs{i}));
    
    % Fixed legend (only subjects)
    hSubjects = gobjects(nSubjects,1);
    for s = 1:nSubjects
        hSubjects(s) = scatter(nan, nan, 50, customColors(s,:), 'filled');
    end
    legH = legend(hSubjects, subjectLegendEntries, 'Location','northeast');
    
    % If FDR-corrected p < 0.05, draw significance line between xPositions(1) and (2)
    if p_fdr_merged(i) < 0.05
        y_sig = max([max(groupData(i).unique_pose), max(groupData(i).unique_seg)]) + 0.01;
        hLine = line([xPositions(1) xPositions(2)], [y_sig y_sig], 'Color', 'k', 'LineWidth', 2);
        set(hLine, 'HandleVisibility','off');
        text(mean([xPositions(1) xPositions(2)]), y_sig+0.005, '*', 'HorizontalAlignment','center', 'FontSize',14, 'Color','k');
    end
    
    % Annotation: report that a paired t-test was done, showing uncorrected and FDR p-values.
    % For the p_fdr part, use 'Interpreter','none' so that the text is written literally.
    annotation('textbox', [0.05, 0.85, 0.4, 0.1], 'String', ...
        sprintf('Paired t-test: uncorrected p = %.4f, FDR p = %.4f (N = %d comparisons)', rawP_merged(i), p_fdr_merged(i), N_merged), ...
        'FitBoxToText', 'on', 'BackgroundColor', 'white', 'EdgeColor', 'black', 'FontSize',10, 'LineStyle', '-', 'Interpreter','none');
    
    % Fixed significance legend (if any ROI is significant)
    if any(p_fdr_merged < 0.05)
        legPos = get(legH, 'Position');
        annotation('textbox', [legPos(1), legPos(2)-legPos(4)-0.01, legPos(3), legPos(4)], ...
            'String', 'Significance Legend: *: p < 0.05', ...
            'FitBoxToText', 'on', 'BackgroundColor', 'white', 'EdgeColor', 'black', 'FontSize',10, 'LineStyle', '-', 'Interpreter','none');
    end
    
    hold off;
    filename = sprintf('MergedROI_%s.png', mergedROIs{i});
    saveas(fig, fullfile(mergedFolder, filename));
    close(fig);
end

%% Graph 2a: Group Laterality Test (Full Model Only)
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
title('Group Laterality Test (Full Model Only)');

hSubjects = gobjects(nSubjects,1);
for s = 1:nSubjects
    hSubjects(s) = scatter(nan, nan, 60, customColors(s,:), 'filled');
end
legH = legend(hSubjects, subjectLegendEntries, 'Location','northeast');

[~, p_val_full, ~, stats_full] = ttest(leftVals_full, rightVals_full);
significanceText_full = ternary(p_val_full < 0.05, 'statistically significant', 'not statistically significant');
statReport_full = sprintf('Full Model Laterality Test:\nPaired t-test: t(%d)=%.2f, uncorrected p = %.4f.\nThe difference is %s.', ...
    stats_full.df, stats_full.tstat, p_val_full, significanceText_full);
annotation('textbox', [0.05, 0.85, 0.4, 0.1], 'String', statReport_full, ...
    'FitBoxToText', 'on', 'BackgroundColor', 'white', 'EdgeColor', 'black', 'FontSize',10, 'LineStyle', '-');
legPos = get(legH, 'Position');
annotation('textbox', [legPos(1), legPos(2)-legPos(4)-0.01, legPos(3), legPos(4)], ...
    'String', 'Significance Legend: *: p < 0.05', ...
    'FitBoxToText', 'on', 'BackgroundColor', 'white', 'EdgeColor', 'black', 'FontSize',10, 'LineStyle', '-', 'Interpreter','none');
hold off;
filename = 'GroupLaterality_FullModel.png';
saveas(fig, fullfile(groupFullFolder, filename));
close(fig);

%% Graph 2b: Group Laterality Test (Unique Segmentation Only)
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
ylabel('Mean Unique Segmentation R^2');
title('Group Laterality Test (Unique Segmentation Only)');

hSubjects = gobjects(nSubjects,1);
for s = 1:nSubjects
    hSubjects(s) = scatter(nan, nan, 60, customColors(s,:), 'filled');
end
legH = legend(hSubjects, subjectLegendEntries, 'Location','northeast');

[~, p_val_seg, ~, stats_seg] = ttest(leftVals_seg, rightVals_seg);
significanceText_seg = ternary(p_val_seg < 0.05, 'statistically significant', 'not statistically significant');
statReport_seg = sprintf('Unique Segmentation Laterality Test:\nPaired t-test: t(%d)=%.2f, uncorrected p = %.4f.\nThe difference is %s.', ...
    stats_seg.df, stats_seg.tstat, p_val_seg, significanceText_seg);
annotation('textbox', [0.05, 0.85, 0.4, 0.1], 'String', statReport_seg, ...
    'FitBoxToText', 'on', 'BackgroundColor', 'white', 'EdgeColor', 'black', 'FontSize',10, 'LineStyle', '-');
legPos = get(legH, 'Position');
annotation('textbox', [legPos(1), legPos(2)-legPos(4)-0.01, legPos(3), legPos(4)], ...
    'String', 'Significance Legend: *: p < 0.05', ...
    'FitBoxToText', 'on', 'BackgroundColor', 'white', 'EdgeColor', 'black', 'FontSize',10, 'LineStyle', '-', 'Interpreter','none');
hold off;
filename = 'GroupLaterality_UniqueSeg.png';
saveas(fig, fullfile(groupSegFolder, filename));
close(fig);

%% Graph 2c: Group Laterality Test (Unique Pose Only)
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
ylabel('Mean Unique Pose R^2');
title('Group Laterality Test (Unique Pose Only)');

hSubjects = gobjects(nSubjects,1);
for s = 1:nSubjects
    hSubjects(s) = scatter(nan, nan, 60, customColors(s,:), 'filled');
end
legH = legend(hSubjects, subjectLegendEntries, 'Location','northeast');

[~, p_val_pose, ~, stats_pose] = ttest(leftVals_pose, rightVals_pose);
significanceText_pose = ternary(p_val_pose < 0.05, 'statistically significant', 'not statistically significant');
statReport_pose = sprintf('Unique Pose Laterality Test:\nPaired t-test: t(%d)=%.2f, uncorrected p = %.4f.\nThe difference is %s.', ...
    stats_pose.df, stats_pose.tstat, p_val_pose, significanceText_pose);
annotation('textbox', [0.05, 0.85, 0.4, 0.1], 'String', statReport_pose, ...
    'FitBoxToText', 'on', 'BackgroundColor', 'white', 'EdgeColor', 'black', 'FontSize',10, 'LineStyle', '-');
legPos = get(legH, 'Position');
annotation('textbox', [legPos(1), legPos(2)-legPos(4)-0.01, legPos(3), legPos(4)], ...
    'String', 'Significance Legend: *: p < 0.05', ...
    'FitBoxToText', 'on', 'BackgroundColor', 'white', 'EdgeColor', 'black', 'FontSize',10, 'LineStyle', '-', 'Interpreter','none');
hold off;
filename = 'GroupLaterality_UniquePose.png';
saveas(fig, fullfile(groupPoseFolder, filename));
close(fig);

%% Graph 3: Per-ROI Laterality Comparisons
tbl0 = readtable(sprintf(dataPathFormat, subjects(1)));
lrIdx = startsWith(tbl0.ROI, 'l') | startsWith(tbl0.ROI, 'r');
roiNames = cellfun(@(x) x(2:end), tbl0.ROI(lrIdx), 'UniformOutput', false);
uniqueROIs = unique(roiNames);
nComparisons = numel(uniqueROIs);

%% (A) Per-ROI Laterality for Full Model
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
    if p_fdr_full(i) < 0.05
        y_sig = max([max(roiData_full(i).left), max(roiData_full(i).right)]) + 0.01;
        hLine = line([1 2], [y_sig y_sig], 'Color', 'k', 'LineWidth', 2);
        set(hLine, 'HandleVisibility','off');
        text(1.5, y_sig+0.005, '*', 'HorizontalAlignment','center', 'FontSize',14, 'Color','k');
    end
    set(gca, 'XTick', [1 2], 'XTickLabel', {'Left','Right'}, 'XLim', [0.8, 2.2], ...
        'YLim', uniformYLim_G3, 'YTick', uniformYTicks_G3);
    ylabel('Full Model R^2');
    title(sprintf('Per-ROI Laterality (Full Model) for ROI: %s', uniqueROIs{i}));
    hSubjects = gobjects(nSubjects,1);
    for s = 1:nSubjects
        hSubjects(s) = scatter(nan, nan, 60, customColors(s,:), 'filled');
    end
    legH = legend(hSubjects, subjectLegendEntries, 'Location','northeast');
    statReport_full_roi = sprintf('For ROI %s, paired t-test (Full Model): %s, uncorrected p = %.4f, FDR p = %.4f (N = %d)', ...
        uniqueROIs{i}, testStats_full{i}, rawP_full(i), p_fdr_full(i), nComparisons);
    annotation('textbox', [0.05, 0.85, 0.4, 0.1], 'String', statReport_full_roi, ...
        'FitBoxToText', 'on', 'BackgroundColor', 'white', 'EdgeColor', 'black', 'FontSize',10, 'LineStyle', '-');
    legPos = get(legH, 'Position');
    annotation('textbox', [legPos(1), legPos(2)-legPos(4)-0.01, legPos(3), legPos(4)], ...
        'String', 'Significance Legend: *: p < 0.05', ...
        'FitBoxToText', 'on', 'BackgroundColor', 'white', 'EdgeColor', 'black', 'FontSize',10, 'LineStyle', '-', 'Interpreter','none');
    hold off;
    filename = sprintf('PerROI_FullModel_%s.png', uniqueROIs{i});
    saveas(fig, fullfile(perROIFullFolder, filename));
    close(fig);
end

%% (B) Per-ROI Laterality for Unique Segmentation
rawP_seg = nan(nComparisons,1);
roiData_seg = struct();
for i = 1:nComparisons
    roiData_seg(i).name = uniqueROIs{i};
    roiData_seg(i).left = nan(nSubjects,1);
    roiData_seg(i).right = nan(nSubjects,1);
    for s = 1:nSubjects
        tbl = readtable(sprintf(dataPathFormat, subjects(s)));
        idx_left = strcmp(tbl.ROI, ['l' uniqueROIs{i}]);
        idx_right = strcmp(tbl.ROI, ['r' uniqueROIs{i}]);
        if any(idx_left)
            roiData_seg(i).left(s) = tbl.unique_seg(find(idx_left,1));
        end
        if any(idx_right)
            roiData_seg(i).right(s) = tbl.unique_seg(find(idx_right,1));
        end
    end
    valid = ~isnan(roiData_seg(i).left) & ~isnan(roiData_seg(i).right);
    if sum(valid) > 1
        [~, p, ~, stats_i] = ttest(roiData_seg(i).left(valid), roiData_seg(i).right(valid));
        rawP_seg(i) = p;
        testStats_seg{i} = sprintf('t(%d)=%.2f', stats_i.df, stats_i.tstat);
    else
        rawP_seg(i) = NaN;
        testStats_seg{i} = 'n/a';
    end
end
[sortedP_seg, sortIdx] = sort(rawP_seg);
q_seg = sortedP_seg .* nComparisons ./ (1:nComparisons)';
q_seg = cummin(q_seg(end:-1:1));
q_seg = q_seg(end:-1:1);
p_fdr_seg = nan(size(rawP_seg));
p_fdr_seg(sortIdx) = q_seg;
clear sortedP_seg sortIdx q_seg

for i = 1:nComparisons
    fig = figure('units','normalized','outerposition',[0 0 1 1]);
    hold on;
    for s = 1:nSubjects
        if ~isnan(roiData_seg(i).left(s)) && ~isnan(roiData_seg(i).right(s))
            scatter(1, roiData_seg(i).left(s), 60, customColors(s,:), 'filled');
            scatter(2, roiData_seg(i).right(s), 60, customColors(s,:), 'filled');
            plot([1 2], [roiData_seg(i).left(s) roiData_seg(i).right(s)], 'Color', [0.7 0.7 0.7], 'HandleVisibility','off');
        end
    end
    m_left_seg_roi = nanmean(roiData_seg(i).left);
    m_right_seg_roi = nanmean(roiData_seg(i).right);
    se_left_seg_roi = nanstd(roiData_seg(i).left) / sqrt(sum(~isnan(roiData_seg(i).left)));
    se_right_seg_roi = nanstd(roiData_seg(i).right) / sqrt(sum(~isnan(roiData_seg(i).right)));
    errorbar(1, m_left_seg_roi, se_left_seg_roi, 'o', 'Color', blackColor, 'MarkerFaceColor', blackColor, 'MarkerSize',10, 'LineWidth',2, 'CapSize',10);
    errorbar(2, m_right_seg_roi, se_right_seg_roi, 'o', 'Color', blackColor, 'MarkerFaceColor', blackColor, 'MarkerSize',10, 'LineWidth',2, 'CapSize',10);
    line([1 2], [m_left_seg_roi, m_right_seg_roi], 'Color', blackColor, 'LineWidth', 3, 'HandleVisibility','off');
    if p_fdr_seg(i) < 0.05
        y_sig = max([max(roiData_seg(i).left), max(roiData_seg(i).right)]) + 0.01;
        hLine = line([1 2], [y_sig y_sig], 'Color', 'k', 'LineWidth', 2);
        set(hLine, 'HandleVisibility','off');
        text(1.5, y_sig+0.005, '*', 'HorizontalAlignment','center', 'FontSize',14, 'Color','k');
    end
    set(gca, 'XTick', [1 2], 'XTickLabel', {'Left','Right'}, 'XLim', [0.8, 2.2], ...
        'YLim', uniformYLim_G3, 'YTick', uniformYTicks_G3);
    ylabel('Unique Segmentation R^2');
    title(sprintf('Per-ROI Laterality (Unique Segmentation) for ROI: %s', uniqueROIs{i}));
    hSubjects = gobjects(nSubjects,1);
    for s = 1:nSubjects
        hSubjects(s) = scatter(nan, nan, 60, customColors(s,:), 'filled');
    end
    legH = legend(hSubjects, subjectLegendEntries, 'Location','northeast');
    statReport_seg_roi = sprintf('For ROI %s, paired t-test (Unique Segmentation): %s, uncorrected p = %.4f, FDR p = %.4f (N = %d)', ...
        uniqueROIs{i}, testStats_seg{i}, rawP_seg(i), p_fdr_seg(i), nComparisons);
    annotation('textbox', [0.05, 0.85, 0.4, 0.1], 'String', statReport_seg_roi, ...
        'FitBoxToText', 'on', 'BackgroundColor', 'white', 'EdgeColor', 'black', 'FontSize',10, 'LineStyle', '-');
    legPos = get(legH, 'Position');
    annotation('textbox', [legPos(1), legPos(2)-legPos(4)-0.01, legPos(3), legPos(4)], ...
        'String', 'Significance Legend: *: p < 0.05', ...
        'FitBoxToText', 'on', 'BackgroundColor', 'white', 'EdgeColor', 'black', 'FontSize',10, 'LineStyle', '-', 'Interpreter','none');
    hold off;
    filename = sprintf('PerROI_UniqueSeg_%s.png', uniqueROIs{i});
    saveas(fig, fullfile(perROISegFolder, filename));
    close(fig);
end

%% (C) Per-ROI Laterality for Unique Pose
rawP_pose = nan(nComparisons,1);
roiData_pose = struct();
for i = 1:nComparisons
    roiData_pose(i).name = uniqueROIs{i};
    roiData_pose(i).left = nan(nSubjects,1);
    roiData_pose(i).right = nan(nSubjects,1);
    for s = 1:nSubjects
        tbl = readtable(sprintf(dataPathFormat, subjects(s)));
        idx_left = strcmp(tbl.ROI, ['l' uniqueROIs{i}]);
        idx_right = strcmp(tbl.ROI, ['r' uniqueROIs{i}]);
        if any(idx_left)
            roiData_pose(i).left(s) = tbl.unique_pose(find(idx_left,1));
        end
        if any(idx_right)
            roiData_pose(i).right(s) = tbl.unique_pose(find(idx_right,1));
        end
    end
    valid = ~isnan(roiData_pose(i).left) & ~isnan(roiData_pose(i).right);
    if sum(valid) > 1
        [~, p, ~, stats_i] = ttest(roiData_pose(i).left(valid), roiData_pose(i).right(valid));
        rawP_pose(i) = p;
        testStats_pose{i} = sprintf('t(%d)=%.2f', stats_i.df, stats_i.tstat);
    else
        rawP_pose(i) = NaN;
        testStats_pose{i} = 'n/a';
    end
end
[sortedP_pose, sortIdx] = sort(rawP_pose);
q_pose = sortedP_pose .* nComparisons ./ (1:nComparisons)';
q_pose = cummin(q_pose(end:-1:1));
q_pose = q_pose(end:-1:1);
p_fdr_pose = nan(size(rawP_pose));
p_fdr_pose(sortIdx) = q_pose;
clear sortedP_pose sortIdx q_pose

for i = 1:nComparisons
    fig = figure('units','normalized','outerposition',[0 0 1 1]);
    hold on;
    for s = 1:nSubjects
        if ~isnan(roiData_pose(i).left(s)) && ~isnan(roiData_pose(i).right(s))
            scatter(1, roiData_pose(i).left(s), 60, customColors(s,:), 'filled');
            scatter(2, roiData_pose(i).right(s), 60, customColors(s,:), 'filled');
            plot([1 2], [roiData_pose(i).left(s) roiData_pose(i).right(s)], 'Color', [0.7 0.7 0.7], 'HandleVisibility','off');
        end
    end
    m_left_pose_roi = nanmean(roiData_pose(i).left);
    m_right_pose_roi = nanmean(roiData_pose(i).right);
    se_left_pose_roi = nanstd(roiData_pose(i).left) / sqrt(sum(~isnan(roiData_pose(i).left)));
    se_right_pose_roi = nanstd(roiData_pose(i).right) / sqrt(sum(~isnan(roiData_pose(i).right)));
    errorbar(1, m_left_pose_roi, se_left_pose_roi, 'o', 'Color', blackColor, 'MarkerFaceColor', blackColor, 'MarkerSize',10, 'LineWidth',2, 'CapSize',10);
    errorbar(2, m_right_pose_roi, se_right_pose_roi, 'o', 'Color', blackColor, 'MarkerFaceColor', blackColor, 'MarkerSize',10, 'LineWidth',2, 'CapSize',10);
    line([1 2], [m_left_pose_roi, m_right_pose_roi], 'Color', blackColor, 'LineWidth', 3, 'HandleVisibility','off');
    if p_fdr_pose(i) < 0.05
        y_sig = max([max(roiData_pose(i).left), max(roiData_pose(i).right)]) + 0.01;
        hLine = line([1 2], [y_sig y_sig], 'Color', 'k', 'LineWidth', 2);
        set(hLine, 'HandleVisibility','off');
        text(1.5, y_sig+0.005, '*', 'HorizontalAlignment','center', 'FontSize',14, 'Color','k');
    end
    set(gca, 'XTick', [1 2], 'XTickLabel', {'Left','Right'}, 'XLim', [0.8, 2.2], ...
        'YLim', uniformYLim_G3, 'YTick', uniformYTicks_G3);
    ylabel('Unique Pose R^2');
    title(sprintf('Per-ROI Laterality (Unique Pose) for ROI: %s', uniqueROIs{i}));
    hSubjects = gobjects(nSubjects,1);
    for s = 1:nSubjects
        hSubjects(s) = scatter(nan, nan, 60, customColors(s,:), 'filled');
    end
    legH = legend(hSubjects, subjectLegendEntries, 'Location','northeast');
    statReport_pose_roi = sprintf('For ROI %s, paired t-test (Unique Pose): %s, uncorrected p = %.4f, FDR p = %.4f (N = %d)', ...
        uniqueROIs{i}, testStats_pose{i}, rawP_pose(i), p_fdr_pose(i), nComparisons);
    annotation('textbox', [0.05, 0.85, 0.4, 0.1], 'String', statReport_pose_roi, ...
        'FitBoxToText', 'on', 'BackgroundColor', 'white', 'EdgeColor', 'black', 'FontSize',10, 'LineStyle', '-');
    legPos = get(legH, 'Position');
    annotation('textbox', [legPos(1), legPos(2)-legPos(4)-0.01, legPos(3), legPos(4)], ...
        'String', 'Significance Legend: *: p < 0.05', ...
        'FitBoxToText', 'on', 'BackgroundColor', 'white', 'EdgeColor', 'black', 'FontSize',10, 'LineStyle', '-', 'Interpreter','none');
    hold off;
    filename = sprintf('PerROI_UniquePose_%s.png', uniqueROIs{i});
    saveas(fig, fullfile(perROIPoseFolder, filename));
    close(fig);
end

disp('All graphs generated.');

% --- Helper function for inline conditional text ---
function out = ternary(cond, trueText, falseText)
    if cond
        out = trueText;
    else
        out = falseText;
    end
end
