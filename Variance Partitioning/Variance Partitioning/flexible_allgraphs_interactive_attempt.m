%% Full MATLAB Script for Group-Level Graphs with Menu (Fully Expanded, Flexible)
clear; clc;

%% --- MENU SELECTION ---
choice = menu('Select Graph Category',...
    'Merged ROI Graphs (Statistical Version)',...
    'All ROI Graphs (New Version: Unique Laterality Metrics)',...
    'Group Laterality Graphs',...
    'Per-ROI Laterality Graphs',...
    'All Categories');

outputChoice = menu('Select Output Mode',...
    'Display Only (Do Not Save)',...
    'Save Figures to Folder');

saveFigures = (outputChoice == 2);

%% --- COMMON CONFIGURATION ---
subjects = 1:8;
nSubjects = numel(subjects);
dataPathFormat = 'D:\\ML_project\\Variance\\var_excel\\sapiens\\lowercase_excels\\subject_%d_variance_partitioning.xlsx';

brightGreen  = [0, 0.8, 0];
strongOrange = [1, 0.6, 0];
blackColor   = [0, 0, 0];
redColor     = [1, 0, 0];
greyColor    = [0.5, 0.5, 0.5];

customColors = [1 0 0;
                0 0 1;
                0 1 0;
                1 0.5 0;
                0.5 0 0.5;
                0 1 1;
                1 0 1;
                0.6 0.4 0.2];

subjectLegendEntries = arrayfun(@(s) sprintf('Subject %d', s), subjects, 'UniformOutput', false);

uniformYTicks_G1 = -0.1:0.02:0.2;
uniformYLim_G1 = [min(uniformYTicks_G1), max(uniformYTicks_G1)];
uniformYTicks_G2 = uniformYTicks_G1;
uniformYLim_G2 = uniformYLim_G1;
uniformYTicks_G3 = uniformYTicks_G1;
uniformYLim_G3 = uniformYLim_G1;

%% --- OUTPUT FOLDER SETUP ---
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

%% --- SECTION 1: Merged ROI Graphs (Statistical Version with Flexibility) ---
if choice == 1 || choice == 5
    tmpTable = readtable(sprintf(dataPathFormat, subjects(1)));
    mergedIdx = startsWith(tmpTable.ROI, 'merged');
    mergedROIs = unique(cellfun(@(x) erase(x, 'merged'), tmpTable.ROI(mergedIdx), 'UniformOutput', false));
    
    has_rpose = any(strcmpi(tmpTable.Properties.VariableNames, 'unique_rpose'));
    has_rseg  = any(strcmpi(tmpTable.Properties.VariableNames, 'unique_rseg'));

    groupData = struct();
    for i = 1:numel(mergedROIs)
        groupData(i).name = mergedROIs{i};
        groupData(i).unique_pose = nan(nSubjects,1);
        groupData(i).unique_seg  = nan(nSubjects,1);
        if has_rpose, groupData(i).unique_rpose = nan(nSubjects,1); end
        if has_rseg,  groupData(i).unique_rseg  = nan(nSubjects,1); end
        groupData(i).full_R2     = nan(nSubjects,1);
    end

    for s = 1:nSubjects
        tbl = readtable(sprintf(dataPathFormat, subjects(s)));
        for i = 1:numel(mergedROIs)
            roiLabel = ['merged' mergedROIs{i}];
            idx = strcmp(tbl.ROI, roiLabel);
            if any(idx)
                groupData(i).unique_pose(s) = tbl.unique_pose(find(idx,1));
                groupData(i).unique_seg(s)  = tbl.unique_seg(find(idx,1));
                if has_rpose, groupData(i).unique_rpose(s) = tbl.unique_rpose(find(idx,1)); end
                if has_rseg,  groupData(i).unique_rseg(s)  = tbl.unique_rseg(find(idx,1));  end
                groupData(i).full_R2(s)     = tbl.Full_R2(find(idx,1));
            end
        end
    end

    rawP_merged = nan(numel(mergedROIs),1);
    for i = 1:numel(mergedROIs)
        [~, p] = ttest(groupData(i).unique_pose, groupData(i).unique_seg);
        rawP_merged(i) = p;
    end

    N_merged = numel(rawP_merged);
    [sortedP, sortIdx] = sort(rawP_merged);
    q = sortedP .* N_merged ./ (1:N_merged)';
    q = cummin(q(end:-1:1));
    q = q(end:-1:1);
    p_fdr_merged = nan(size(rawP_merged));
    p_fdr_merged(sortIdx) = q;

    for i = 1:numel(mergedROIs)
        xPositions = [1, 1.8]; xLabels = {'Pose Estimation','Body Segmentation'};
        if has_rpose, xPositions(end+1) = xPositions(end)+0.8; xLabels{end+1} = 'Random Pose'; end
        if has_rseg,  xPositions(end+1) = xPositions(end)+0.8; xLabels{end+1} = 'Random Seg'; end
        xPositions(end+1) = xPositions(end)+0.8; xLabels{end+1} = 'Full Model';

        fig = figure('units','normalized','outerposition',[0 0 1 1]); hold on;
  
        for s = 1:nSubjects
            scatter(xPositions(1), groupData(i).unique_pose(s), 50, customColors(s,:), 'filled');
            scatter(xPositions(2), groupData(i).unique_seg(s), 50, customColors(s,:), 'filled');
            pIdx = 3;
            if has_rpose, scatter(xPositions(pIdx), groupData(i).unique_rpose(s), 50, customColors(s,:), 'filled'); pIdx = pIdx+1; end
            if has_rseg,  scatter(xPositions(pIdx), groupData(i).unique_rseg(s), 50, customColors(s,:), 'filled');  pIdx = pIdx+1; end
            scatter(xPositions(pIdx), groupData(i).full_R2(s), 50, customColors(s,:), 'filled');
        end

        m_vals = [nanmean(groupData(i).unique_pose), nanmean(groupData(i).unique_seg)];
        se_vals = [nanstd(groupData(i).unique_pose), nanstd(groupData(i).unique_seg)] / sqrt(nSubjects);
        if has_rpose, m_vals(end+1) = nanmean(groupData(i).unique_rpose); se_vals(end+1) = nanstd(groupData(i).unique_rpose)/sqrt(nSubjects); end
        if has_rseg,  m_vals(end+1) = nanmean(groupData(i).unique_rseg);  se_vals(end+1) = nanstd(groupData(i).unique_rseg)/sqrt(nSubjects); end
        m_vals(end+1) = nanmean(groupData(i).full_R2);
        se_vals(end+1) = nanstd(groupData(i).full_R2)/sqrt(nSubjects);

        for k = 1:length(xPositions)
            errorbar(xPositions(k), m_vals(k), se_vals(k), 'o', 'Color', blackColor, ...
                'MarkerFaceColor', 'none', 'MarkerSize',10, 'LineWidth',2, 'CapSize',10);
        end

        set(gca, 'XTick', xPositions, 'XTickLabel', xLabels, ...
            'XLim', [xPositions(1)-0.5, xPositions(end)+0.5], ...
            'YLim', uniformYLim_G1, 'YTick', uniformYTicks_G1);
        ylabel('Variance Explained (R^2)');
        title(sprintf('Merged ROI: %s', mergedROIs{i}));

        hSubjects = gobjects(nSubjects,1);
        for s = 1:nSubjects
            hSubjects(s) = scatter(nan, nan, 50, customColors(s,:), 'filled');
        end
        legend(hSubjects, subjectLegendEntries, 'Location','northeast');

        if p_fdr_merged(i) < 0.05
            y_sig = max([groupData(i).unique_pose; groupData(i).unique_seg]) + 0.01;
            line([xPositions(1) xPositions(2)], [y_sig y_sig], 'Color', 'k', 'LineWidth', 2);
            text(mean([xPositions(1), xPositions(2)]), y_sig+0.005, '*', ...
                'HorizontalAlignment','center', 'FontSize',14, 'Color','k');
        end

        annotation('textbox', [0.05, 0.85, 0.4, 0.1], 'String', ...
            sprintf('Paired t-test: p = %.4f, FDR p = %.4f', rawP_merged(i), p_fdr_merged(i)), ...
            'FitBoxToText', 'on', 'BackgroundColor', 'white', 'EdgeColor', 'black');

        
        hline = line(xlim(), [0 0], 'LineStyle', '--', 'Color', [0 0 0], 'LineWidth', 1);
        set(get(get(hline,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
        hold off;

        if saveFigures
            saveas(fig, fullfile(mergedStatFolder, sprintf('MergedROI_Stat_%s.png', mergedROIs{i})));
            close(fig);
        else
            pause(1);
        end
    end
end

%% --- SECTION 2: All ROI Graphs (Flexible) ---
if choice == 2 || choice == 5
    tmpTable = readtable(sprintf(dataPathFormat, subjects(1)));
    allROIs = unique(tmpTable.ROI);

    has_rpose = any(strcmpi(tmpTable.Properties.VariableNames, 'unique_rpose'));
    has_rseg  = any(strcmpi(tmpTable.Properties.VariableNames, 'unique_rseg'));

    xPositions = [1, 2];
    xLabels = {'Pose Estimation', 'Body Segmentation'};
    if has_rpose, xPositions(end+1) = xPositions(end)+1; xLabels{end+1} = 'Random Pose'; end
    if has_rseg,  xPositions(end+1) = xPositions(end)+1; xLabels{end+1} = 'Random Seg'; end
    xlim_new = [min(xPositions)-0.5, max(xPositions)+0.5];

    for i = 1:numel(allROIs)
        poseVals = nan(nSubjects,1);
        segVals  = nan(nSubjects,1);
        if has_rpose, rposeVals = nan(nSubjects,1); end
        if has_rseg,  rsegVals  = nan(nSubjects,1); end

        for s = 1:nSubjects
            tbl = readtable(sprintf(dataPathFormat, subjects(s)));
            idx = strcmp(tbl.ROI, allROIs{i});
            if any(idx)
                poseVals(s)  = tbl.unique_pose(find(idx,1));
                segVals(s)   = tbl.unique_seg(find(idx,1));
                if has_rpose, rposeVals(s) = tbl.unique_rpose(find(idx,1)); end
                if has_rseg,  rsegVals(s)  = tbl.unique_rseg(find(idx,1));  end
            end
        end

        fig = figure('units','normalized','outerposition',[0 0 1 1]); hold on;
        for s = 1:nSubjects
            scatter(xPositions(1), poseVals(s), 50, customColors(s,:), 'filled');
            scatter(xPositions(2), segVals(s), 50, customColors(s,:), 'filled');
            pIdx = 3;
            if has_rpose
                scatter(xPositions(pIdx), rposeVals(s), 50, customColors(s,:), 'filled'); 
                pIdx = pIdx+1; 
            end
            if has_rseg
                scatter(xPositions(pIdx), rsegVals(s), 50, customColors(s,:), 'filled'); 
            end
        end

        m_vals = [nanmean(poseVals), nanmean(segVals)];
        se_vals = [nanstd(poseVals)/sqrt(sum(~isnan(poseVals))), nanstd(segVals)/sqrt(sum(~isnan(segVals)))];
        if has_rpose
            m_vals(end+1) = nanmean(rposeVals);
            se_vals(end+1) = nanstd(rposeVals)/sqrt(sum(~isnan(rposeVals)));
        end
        if has_rseg
            m_vals(end+1) = nanmean(rsegVals);
            se_vals(end+1) = nanstd(rsegVals)/sqrt(sum(~isnan(rsegVals)));
        end

        for k = 1:length(xPositions)
            errorbar(xPositions(k), m_vals(k), se_vals(k), 'o', 'Color', blackColor, ...
                'MarkerFaceColor', 'none', 'MarkerSize',10, 'LineWidth',2, 'CapSize',10);
        end

        yline(0, 'k', 'LineWidth', 1);
        set(gca, 'XTick', xPositions, 'XTickLabel', xLabels, ...
            'XLim', xlim_new, 'YLim', uniformYLim_G1, 'YTick', uniformYTicks_G1);
        ylabel('Unique Variance Explained (R^2)');
        title(sprintf('Unique Variance Explained for ROI: %s', allROIs{i}), 'Interpreter','none');

        hSubjects = gobjects(nSubjects,1);
        for s = 1:nSubjects
            hSubjects(s) = scatter(nan, nan, 50, customColors(s,:), 'filled');
        end
        legend(hSubjects, subjectLegendEntries, 'Location','northeast');


        hline = line(xlim(), [0 0], 'LineStyle', '--', 'Color', [0 0 0], 'LineWidth', 1);
        set(get(get(hline,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
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

%% --- SECTION 3: Group Laterality Graphs (Flexible) ---
if choice == 3 || choice == 5
    measures = {'Full_R2', 'unique_pose', 'unique_seg'};
    measureLabels = {'Full Model', 'Unique Pose', 'Unique Segmentation'};

    for m = 1:length(measures)
        leftVals = nan(nSubjects,1);
        rightVals = nan(nSubjects,1);
        for s = 1:nSubjects
            tbl = readtable(sprintf(dataPathFormat, subjects(s)));
            lIdx = startsWith(tbl.ROI, 'l');
            rIdx = startsWith(tbl.ROI, 'r');
            if any(lIdx)
                leftVals(s) = nanmean(tbl.(measures{m})(lIdx));
            end
            if any(rIdx)
                rightVals(s) = nanmean(tbl.(measures{m})(rIdx));
            end
        end

        fig = figure('units','normalized','outerposition',[0 0 1 1]); hold on;
        for s = 1:nSubjects
            scatter(1, leftVals(s), 60, customColors(s,:), 'filled');
            scatter(2, rightVals(s), 60, customColors(s,:), 'filled');
            plot([1 2], [leftVals(s) rightVals(s)], 'Color', [0.7 0.7 0.7]);
        end

        m_left = nanmean(leftVals);
        m_right = nanmean(rightVals);
        se_left = nanstd(leftVals)/sqrt(sum(~isnan(leftVals)));
        se_right = nanstd(rightVals)/sqrt(sum(~isnan(rightVals)));

        errorbar(1, m_left, se_left, 'o', 'Color', blackColor, ...
            'MarkerFaceColor', 'none', 'MarkerSize',10, 'LineWidth',2, 'CapSize',10);
        errorbar(2, m_right, se_right, 'o', 'Color', blackColor, ...
            'MarkerFaceColor', 'none', 'MarkerSize',10, 'LineWidth',2, 'CapSize',10);
        line([1,2],[m_left, m_right], 'Color', blackColor, 'LineWidth', 3);

        set(gca, 'XTick', [1 2], 'XTickLabel', {'Left ROIs','Right ROIs'}, ...
            'XLim', [0.8, 2.2], 'YLim', uniformYLim_G2, 'YTick', uniformYTicks_G2);
        ylabel(sprintf('Mean %s (R^2)', measureLabels{m}));
        title(sprintf('Group Laterality Test (%s)', measureLabels{m}), 'Interpreter','none');

        hSubjects = gobjects(nSubjects,1);
        for s = 1:nSubjects
            hSubjects(s) = scatter(nan, nan, 60, customColors(s,:), 'filled');
        end
        legend(hSubjects, subjectLegendEntries, 'Location','northeast');

        [~, p_val, ~, stats] = ttest(leftVals, rightVals);
        sigText = 'n.s.';
        if p_val < 0.05, sigText = '*'; end

        annotation('textbox', [0.05, 0.85, 0.4, 0.1], 'String', ...
            sprintf('t(%d)=%.2f, p=%.4f %s', stats.df, stats.tstat, p_val, sigText),...
            'FitBoxToText', 'on', 'BackgroundColor', 'white', 'EdgeColor', 'black');

        hline = line(xlim(), [0 0], 'LineStyle', '--', 'Color', [0 0 0], 'LineWidth', 1);
        set(get(get(hline,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
        hold off;

        if saveFigures
            filename = sprintf('GroupLaterality_%s.png', measureLabels{m});
            saveas(fig, fullfile(groupFolder, filename));
            close(fig);
        else
            pause(1);
        end
    end
end

%% --- SECTION 4: Per-ROI Laterality Graphs (Flexible) ---
if choice == 4 || choice == 5
    tmpTable = readtable(sprintf(dataPathFormat, subjects(1)));
    lrIdx = startsWith(tmpTable.ROI, 'l') | startsWith(tmpTable.ROI, 'r');
    roiNames = cellfun(@(x) x(2:end), tmpTable.ROI(lrIdx), 'UniformOutput', false);
    uniqueROIs = unique(roiNames);
    nComparisons = numel(uniqueROIs);
    
    measures = {'Full_R2', 'unique_pose', 'unique_seg'};
    measureLabels = {'Full Model', 'Unique Pose', 'Unique Segmentation'};
    
    for m = 1:length(measures)
        rawP = nan(nComparisons,1);
        roiData = struct();
        for i = 1:nComparisons
            roiData(i).name = uniqueROIs{i};
            roiData(i).left = nan(nSubjects,1);
            roiData(i).right = nan(nSubjects,1);
            for s = 1:nSubjects
                tbl = readtable(sprintf(dataPathFormat, subjects(s)));
                idx_left = strcmp(tbl.ROI, ['l' uniqueROIs{i}]);
                idx_right = strcmp(tbl.ROI, ['r' uniqueROIs{i}]);
                if any(idx_left)
                    roiData(i).left(s) = tbl.(measures{m})(find(idx_left,1));
                end
                if any(idx_right)
                    roiData(i).right(s) = tbl.(measures{m})(find(idx_right,1));
                end
            end
            
            valid = ~isnan(roiData(i).left) & ~isnan(roiData(i).right);
            if sum(valid) > 1
                [~, p, ~, stats_i] = ttest(roiData(i).left(valid), roiData(i).right(valid));
                rawP(i) = p;
                statsText{i} = sprintf('t(%d)=%.2f', stats_i.df, stats_i.tstat);
            else
                rawP(i) = NaN;
                statsText{i} = 'n/a';
            end
        end
        
        [sortedP, sortIdx] = sort(rawP);
        q = sortedP .* nComparisons ./ (1:nComparisons)';
        q = cummin(q(end:-1:1));
        q = q(end:-1:1);
        p_fdr = nan(size(rawP));
        p_fdr(sortIdx) = q;

        for i = 1:nComparisons
            fig = figure('units','normalized','outerposition',[0 0 1 1]); hold on;
            for s = 1:nSubjects
                if ~isnan(roiData(i).left(s)) && ~isnan(roiData(i).right(s))
                    scatter(1, roiData(i).left(s), 60, customColors(s,:), 'filled');
                    scatter(2, roiData(i).right(s), 60, customColors(s,:), 'filled');
                    plot([1 2], [roiData(i).left(s) roiData(i).right(s)], 'Color', [0.7 0.7 0.7]);
                end
            end
            
            m_left = nanmean(roiData(i).left);
            m_right = nanmean(roiData(i).right);
            se_left = nanstd(roiData(i).left)/sqrt(sum(~isnan(roiData(i).left)));
            se_right = nanstd(roiData(i).right)/sqrt(sum(~isnan(roiData(i).right)));
            
            errorbar(1, m_left, se_left, 'o', 'Color', blackColor, ...
                'MarkerFaceColor', 'none', 'MarkerSize',10, 'LineWidth',2, 'CapSize',10);
            errorbar(2, m_right, se_right, 'o', 'Color', blackColor, ...
                'MarkerFaceColor', 'none', 'MarkerSize',10, 'LineWidth',2, 'CapSize',10);
            line([1,2],[m_left, m_right], 'Color', blackColor, 'LineWidth', 3);
            
            set(gca, 'XTick', [1 2], 'XTickLabel', {'Left','Right'}, ...
                'XLim', [0.8, 2.2], 'YLim', uniformYLim_G3, 'YTick', uniformYTicks_G3);
            ylabel(sprintf('%s R^2', measureLabels{m}));
            title(sprintf('Per-ROI Laterality (%s) for ROI: %s', measureLabels{m}, uniqueROIs{i}), 'Interpreter','none');
            
            hSubjects = gobjects(nSubjects,1);
            for s = 1:nSubjects
                hSubjects(s) = scatter(nan, nan, 60, customColors(s,:), 'filled');
            end
            legend(hSubjects, subjectLegendEntries, 'Location','northeast');
            
            annotation('textbox', [0.05, 0.85, 0.4, 0.1], 'String', ...
                sprintf('%s: %s, p=%.4f, FDR p=%.4f', measureLabels{m}, statsText{i}, rawP(i), p_fdr(i)),...
                'FitBoxToText', 'on', 'BackgroundColor', 'white', 'EdgeColor', 'black');

            hline = line(xlim(), [0 0], 'LineStyle', '--', 'Color', [0 0 0], 'LineWidth', 1);
            set(get(get(hline,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
            hold off;
            
            if saveFigures
                filename = sprintf('PerROI_%s_%s.png', measureLabels{m}, uniqueROIs{i});
                saveas(fig, fullfile(perROIFolder, filename));
                close(fig);
            else
                pause(1);
            end
        end
    end
end
