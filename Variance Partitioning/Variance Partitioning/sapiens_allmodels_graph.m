% sapiens_allmodels_graph_robust_final.m
clear; clc;

%% --- Output Mode ---
choice = menu('Output Mode', 'Display Only (No Save)', 'Save Figures to Folder');
saveFigures = (choice == 2);

%% --- Config ---
subjects = 1:8;
nSubjects = numel(subjects);
dataPathFormat = 'D:\\ML_project\\Variance\\var_excel\\sapiens_allmodels\\subject_%d_variance_partitioning_updated.xlsx';
outputFolder = 'GraphOutputs\\Merged_ROIs_Statistical_4Models';
if saveFigures && ~exist(outputFolder, 'dir'), mkdir(outputFolder); end

modelLabels = {'Pose Estimation', 'Body Segmentation', 'Depth Estimation', 'Surface Normal Estimation'};
modelFields = {'unique_pose', 'unique_seg', 'unique_depth', 'unique_normals'};
nModels = numel(modelLabels);

customColors = [1 0 0;
                0 0 1;
                0 1 0;
                1 0.5 0;
                0.5 0 0.5;
                0 1 1;
                1 0 1;
                0.6 0.4 0.2];
blackColor = [0 0 0];

uniformYTicks = -0.05:0.05:0.25;
uniformYLim = [min(uniformYTicks), max(uniformYTicks)];

%% --- Get ROIs
tmp = readtable(sprintf(dataPathFormat, subjects(1)));
mergedROIs = unique(cellstr(tmp.ROI));

%% --- Loop over ROIs
for i = 1:numel(mergedROIs)
    data = nan(nSubjects, nModels);
    fullModelData = nan(nSubjects,1);

    for s = 1:nSubjects
        tbl = readtable(sprintf(dataPathFormat, subjects(s)));
        idx = strcmp(tbl.ROI, mergedROIs{i});
        if any(idx)
            for m = 1:nModels
                data(s, m) = tbl.(modelFields{m})(idx);
            end
            fullModelData(s) = tbl.Full_R2(idx);
        end
    end

    if all(isnan(data), 'all')
        continue;
    end

    %% --- ANOVA (only 4 models)
    T = array2table(data, 'VariableNames', {'Pose','Seg','Depth','Normals'});
    T.Subject = categorical((1:nSubjects)');
    rm = fitrm(T, 'Pose-Normals ~ 1', 'WithinDesign', ...
        table(modelLabels', 'VariableNames', {'Model'}));
    ranova_tbl = ranova(rm);
    p_anova = ranova_tbl.pValue(1);
    df1 = ranova_tbl.DF(1); df2 = ranova_tbl.DF(2);

    %% --- Post Hoc with FDR
    pairwiseStr = "";
    nComparisons = 0;
    sigPosthoc = [];

    if p_anova < 0.05
        posthoc = multcompare(rm, 'Model', 'ComparisonType', 'tukey-kramer');
        raw_p = posthoc.pValue;
        nComparisons = numel(raw_p);

        % FDR correction
        [sortedP, sortIdx] = sort(raw_p);
        q = sortedP .* nComparisons ./ (1:nComparisons)';
        q = cummin(q(end:-1:1));
        q = q(end:-1:1);
        p_fdr = nan(size(raw_p)); p_fdr(sortIdx) = q;

        % Filter significant comparisons safely
        for c = 1:length(p_fdr)
            if ~isnan(p_fdr(c)) && p_fdr(c) < 0.05
                row = posthoc.Row{c};
                col = posthoc.Column{c};
                diff = posthoc.Difference(c);
                sigPosthoc(end+1,:) = [find(strcmp(modelLabels,row)), find(strcmp(modelLabels,col))]; %#ok<AGROW>
                pairwiseStr = pairwiseStr + sprintf('%s > %s: q = %.4f\n', row, col, p_fdr(c));
            end
        end
    end

    %% --- Plot
    fig = figure('units','normalized','outerposition',[0.1 0.1 0.75 0.85]); hold on;
    xPos = 1:(nModels + 1);
    xPadding = 0.3;

    % Plot subject values
    for s = 1:nSubjects
        for m = 1:nModels
            scatter(xPos(m), data(s,m), 50, customColors(s,:), 'filled');
        end
        scatter(xPos(end), fullModelData(s), 50, customColors(s,:), 'filled');
    end

    % Mean ± SEM
    allMeans = [mean(data, 'omitnan'), mean(fullModelData, 'omitnan')];
    allSEMs = [std(data, 'omitnan') ./ sqrt(nSubjects), std(fullModelData, 'omitnan')/sqrt(nSubjects)];
    for m = 1:length(xPos)
        errorbar(xPos(m), allMeans(m), allSEMs(m), 'o', ...
            'Color', blackColor, 'MarkerFaceColor', 'none', ...
            'MarkerSize', 10, 'LineWidth', 2, 'CapSize', 10);
    end

    % Significance brackets for post hoc
    yMax = max([data(:); fullModelData(:)]) + 0.02;
    yOffset = 0.015;
    for c = 1:size(sigPosthoc,1)
        x1 = xPos(sigPosthoc(c,1));
        x2 = xPos(sigPosthoc(c,2));
        y = yMax + (c-1)*yOffset;
        line([x1 x2], [y y], 'Color','k', 'LineWidth',1.5);
        text(mean([x1,x2]), y + 0.003, '*', 'FontSize', 12, ...
            'HorizontalAlignment','center', 'Color','k');
    end

    % Format
    xLabels = [modelLabels, {'Full Model'}];
    xlim([xPos(1)-xPadding, xPos(end)+xPadding]);
    ylabel('Unique Variance Explained (R^2)');
    title(sprintf('Variance Partitioning - %s', mergedROIs{i}), 'Interpreter','none');
    set(gca, 'XTick', xPos, 'XTickLabel', xLabels, ...
        'YTick', uniformYTicks, 'YLim', uniformYLim, 'FontSize', 10);

    % Legend
    hSubjects = gobjects(nSubjects,1);
    for s = 1:nSubjects
        hSubjects(s) = scatter(nan, nan, 50, customColors(s,:), 'filled');
    end
    legend(hSubjects, arrayfun(@(x) sprintf('Subject %d', x), subjects, 'UniformOutput', false), ...
        'Location', 'northeastoutside');

    % Stat box outside axes (right side)
    statStr = sprintf(['Repeated-measures ANOVA:\n', ...
        'df = %d, %d | p = %.4f (uncorrected)\n', ...
        'Subjects = %d | Models = %d\n', ...
        'Post-hoc: Tukey-Kramer + FDR (%d comparisons)\n%s'], ...
        df1, df2, p_anova, nSubjects, nModels, nComparisons, pairwiseStr);
    
    annotation('textbox', [0.76, 0.1, 0.22, 0.3], 'String', statStr, ...
        'FitBoxToText', 'on', 'BackgroundColor', 'white', ...
        'EdgeColor', 'black', 'FontSize', 9);

    % Horizontal line (not in legend)
    hline = yline(0, '--', 'Color', [0.5 0.5 0.5], 'LineWidth', 1);
    set(get(get(hline,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');

    hold off;

    %% Save or display
    if saveFigures
        saveas(fig, fullfile(outputFolder, sprintf('ROI_%s_GroupBarplot.png', mergedROIs{i})));
        close(fig);
    else
        pause(1);
    end
end

fprintf('✅ All graphs finished without crashing.\n');
