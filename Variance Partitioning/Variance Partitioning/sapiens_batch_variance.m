% 12.4 version of variance partitioning, for sapiens data

tic;  % Start timer
clear; clc;

%% Configuration
subject = 8;

% List of ROIs (no laterality)
rois = {'EBA', 'FBA_1', 'FBA_2', 'FFA_1', 'FFA_2', 'hV4', 'MTL', 'OFA', ...
        'V1d', 'V1v', 'V2d', 'V2v', 'V3d', 'V3v'};

% Base paths
brain_base = sprintf('D:\\ML_project\\Variance\\new_brain_results\\subject_%d', subject);
dnn_base   = sprintf('D:\\ML_project\\Variance\\sapiens_results\\subject_%d', subject);

% Prepare results table
headers = {'ROI', 'Full_R2', 'Full_RMSE', 'Pose_R2', 'Seg_R2', 'Unique_Pose', 'Unique_Seg', 'Shared_Variance'};
results = [headers];

%% Loop over ROIs
for i = 1:length(rois)
    roi = rois{i};
    roi_label = sprintf('merged%s', roi);  % e.g., mergedEBA
    
    % File paths
    brain_train_path = fullfile(brain_base, 'train', roi, sprintf('%s_merged_data.xlsx', roi));
    brain_val_path   = fullfile(brain_base, 'val', roi, sprintf('%s_merged_data.xlsx', roi));
    
    pose_train_path = fullfile(dnn_base, 'pose', 'train', sprintf('best_rdm_for_%s_concatenated_data.xlsx', roi));
    pose_val_path   = fullfile(dnn_base, 'pose', 'val', sprintf('best_rdm_for_%s_concatenated_data.xlsx', roi));
    
    seg_train_path  = fullfile(dnn_base, 'seg', 'train', sprintf('best_rdm_for_%s_concatenated_data.xlsx', roi));
    seg_val_path    = fullfile(dnn_base, 'seg', 'val', sprintf('best_rdm_for_%s_concatenated_data.xlsx', roi));
    
    % Skip if any files are missing
    if exist(brain_train_path, 'file') ~= 2 || exist(brain_val_path, 'file') ~= 2 || ...
       exist(pose_train_path, 'file') ~= 2 || exist(pose_val_path, 'file') ~= 2 || ...
       exist(seg_train_path, 'file') ~= 2 || exist(seg_val_path, 'file') ~= 2
        fprintf('Skipping %s due to missing file(s).\n', roi_label);
        continue;
    end
    
    fprintf('Processing %s...\n', roi_label);
    
    %% Load Data
    brain_train = readmatrix(brain_train_path);
    pose_train  = readmatrix(pose_train_path);
    seg_train   = readmatrix(seg_train_path);
    
    brain_val = readmatrix(brain_val_path);
    pose_val  = readmatrix(pose_val_path);
    seg_val   = readmatrix(seg_val_path);
    
    % Check for size mismatches
    if length(pose_train) ~= length(seg_train) || length(pose_train) ~= length(brain_train)
        fprintf('Skipping %s due to mismatched training data dimensions.\n', roi_label);
        continue;
    end
    if length(pose_val) ~= length(seg_val) || length(pose_val) ~= length(brain_val)
        fprintf('Skipping %s due to mismatched validation data dimensions.\n', roi_label);
        continue;
    end
    
    %% Fit Models
    X_train = [pose_train, seg_train];
    Y_train = brain_train;
    X_val   = [pose_val, seg_val];
    Y_val   = brain_val;
    
    % Full model
    full_model = fitlm(X_train, Y_train);
    Y_pred_full = predict(full_model, X_val);
    r2_full = 1 - sum((Y_val - Y_pred_full).^2) / sum((Y_val - mean(Y_val)).^2);
    rmse_full = sqrt(mean((Y_val - Y_pred_full).^2));
    
    % Pose-only
    model_pose = fitlm(pose_train, Y_train);
    Y_pred_pose = predict(model_pose, pose_val);
    r2_pose = 1 - sum((Y_val - Y_pred_pose).^2) / sum((Y_val - mean(Y_val)).^2);
    
    % Seg-only
    model_seg = fitlm(seg_train, Y_train);
    Y_pred_seg = predict(model_seg, seg_val);
    r2_seg = 1 - sum((Y_val - Y_pred_seg).^2) / sum((Y_val - mean(Y_val)).^2);
    
    % Variance components
    unique_pose = r2_full - r2_seg;
    unique_seg  = r2_full - r2_pose;
    shared_variance = r2_full - unique_pose - unique_seg;
    
    % Output (optional)
    fprintf('%s: Full R2 = %.4f, RMSE = %.4f, Unique Pose = %.4f, Unique Seg = %.4f, Shared = %.4f\n', ...
            roi_label, r2_full, rmse_full, unique_pose, unique_seg, shared_variance);
    
    %% Save Row
    newRow = {roi_label, r2_full, rmse_full, r2_pose, r2_seg, unique_pose, unique_seg, shared_variance};
    results = [results; newRow];
end

%% Save to Excel
resultsTable = cell2table(results(2:end,:), 'VariableNames', results(1,:));

output_folder = 'D:\ML_project\Variance\var_excel';
if ~exist(output_folder, 'dir')
    mkdir(output_folder);
end

output_file = fullfile(output_folder, sprintf('subject_%d_variance_partitioning_merged.xlsx', subject));
writetable(resultsTable, output_file);
fprintf('Results saved to %s\n', output_file);

toc;
