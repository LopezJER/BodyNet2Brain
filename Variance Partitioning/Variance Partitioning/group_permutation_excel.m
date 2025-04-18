%% -------------------- GROUP LEVEL SUMMARY BUILDER --------------------
% Aggregates all subject summary Excel files into one master table
% Computes group mean & std for each stat across ROIs

clear; clc;

%% --- CONFIGURATION ---
nSubjects = 8;
input_base = 'D:/ML_project/Variance/var_excel/sapiens/subject_excels/';
summary_files = dir(fullfile(input_base, '**', 'subject_*_summary.xlsx'));

% Initialize storage
allData = table();

for i = 1:length(summary_files)
    fname = fullfile(summary_files(i).folder, summary_files(i).name);
    subj_data = readtable(fname);

    % Extract subject number from filename
    s = regexp(fname, 'subject_(\d+)_summary.xlsx', 'tokens');
    subject_id = str2double(s{1}{1});

    % Add subject column
    subj_col = repmat(subject_id, height(subj_data), 1);
    subj_data.Subject = subj_col;

    allData = [allData; subj_data];
end

% Reorder columns: Subject, ROI, then stats
allData = movevars(allData, 'Subject', 'Before', 1);

%% --- SAVE COMBINED DATA ---
group_dir = fullfile(input_base, 'group_level_2');
if ~exist(group_dir, 'dir'), mkdir(group_dir); end

combined_file = fullfile(group_dir, 'all_subjects_summary.xlsx');
writetable(allData, combined_file);

%% --- STATS PER ROI ---
rois = unique(allData.ROI);
stats = {'Unique_Pose','Unique_Seg','Full_R2','Shared'};
pvals = {'p_Unique_Pose','p_Unique_Seg','p_Full_R2','p_Shared'};

summary_stats = cell(length(rois), 1 + 2*length(stats) + 2*length(pvals));
header = {'ROI'};
for s = 1:length(stats)
    header{end+1} = [stats{s} '_Mean'];
    header{end+1} = [stats{s} '_Std'];
end
for s = 1:length(pvals)
    header{end+1} = [pvals{s} '_Mean'];
    header{end+1} = [pvals{s} '_Std'];
end

for r = 1:length(rois)
    roi = rois{r};
    summary_stats{r,1} = roi;
    roi_data = allData(strcmp(allData.ROI, roi), :);
    
    for s = 1:length(stats)
        vals = roi_data.(stats{s});
        summary_stats{r, 2*s}   = mean(vals, 'omitnan');
        summary_stats{r, 2*s+1} = std(vals, 'omitnan');
    end

    offset = 2*length(stats);
    for s = 1:length(pvals)
        vals = roi_data.(pvals{s});
        summary_stats{r, offset + 2*s}   = mean(vals, 'omitnan');
        summary_stats{r, offset + 2*s+1} = std(vals, 'omitnan');
    end
end

% Save stats table
summaryT = cell2table(summary_stats, 'VariableNames', header);
summary_file = fullfile(group_dir, 'group_stats_mean_std.xlsx');
writetable(summaryT, summary_file);

fprintf('\n--- Group-Level Summary Complete ---\nSaved: %s\nSaved: %s\n', combined_file, summary_file);
