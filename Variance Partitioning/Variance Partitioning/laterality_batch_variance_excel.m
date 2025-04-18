%separated lateralities - variance partitioning results saved to excel

%% Extract Variance Partitioning Values and Save to Excel
% This script loops over a list of ROIs and both laterality options for a given
% subject, computes variance partitioning values from regression models, and saves
% the results to an Excel file. Each row is labeled by ROI with laterality (e.g., rEBA).
% If the necessary files for a brain region are not found, that combination is skipped.
% Execution time is displayed via tic/toc.

tic;  % Start timer
clear; clc;

%% Configuration
% Set the subject number (change this to process a different subject)
subject = 1;

% Define list of ROIs and laterality options
rois = {'EBA', 'FBA_1', 'FBA_2', 'FFA_1', 'FFA_2', 'hV4', 'MTL', 'OFA', ...
        'V1d', 'V1v', 'V2d', 'V2v', 'V3d', 'V3v'};
laterality_options = {'r', 'l'};  % 'r' for right, 'l' for left

% Define base directories for brain data and DNN data (pose and seg)
brain_base = sprintf('D:\\ML_project\\Variance\\brain_results\\subject_%d', subject);
dnn_base   = sprintf('D:\\ML_project\\Variance\\sapiens_results\\subject_%d', subject);

% Prepare output results cell array with header row.
results = {};
headers = {'ROI', 'Full_R2', 'Full_RMSE', 'Pose_R2', 'Seg_R2', 'Unique_Pose', 'Unique_Seg', 'Shared_Variance'};
results = [headers];

%% Loop over each ROI and laterality
for i = 1:length(rois)
    roi = rois{i};
    for j = 1:length(laterality_options)
        laterality_input = laterality_options{j};
        
        % Convert laterality letter to full string ('right' or 'left')
        if lower(laterality_input) == 'r'
            laterality_str = 'right';
        elseif lower(laterality_input) == 'l'
            laterality_str = 'left';
        else
            error('Laterality must be ''r'' or ''l''.');
        end
        
        % Create ROI label (e.g., rEBA or lEBA)
        roi_label = sprintf('%s%s', lower(laterality_input), roi);
        
        % Construct file paths:
        % Brain data (training and validation)
        brain_train_path = fullfile(brain_base, 'train', roi, sprintf('%s_%s_data.xlsx', roi, laterality_str));
        brain_val_path   = fullfile(brain_base, 'val', roi, sprintf('%s_%s_data.xlsx', roi, laterality_str));
        
        % Pose (dnn1) data (training and validation)
        pose_train_path = fullfile(dnn_base, 'pose', 'train', sprintf('best_rdm_for_%s_%s_data.xlsx', roi, laterality_str));
        pose_val_path   = fullfile(dnn_base, 'pose', 'val', sprintf('best_rdm_for_%s_%s_data.xlsx', roi, laterality_str));
        
        % Seg (dnn2) data (training and validation)
        seg_train_path  = fullfile(dnn_base, 'seg', 'train', sprintf('best_rdm_for_%s_%s_data.xlsx', roi, laterality_str));
        seg_val_path    = fullfile(dnn_base, 'seg', 'val', sprintf('best_rdm_for_%s_%s_data.xlsx', roi, laterality_str));
        
        % Check if the essential brain files exist; if not, skip this combination.
        if exist(brain_train_path, 'file') ~= 2 || exist(brain_val_path, 'file') ~= 2
            fprintf('Skipping %s because brain data file not found.\n', roi_label);
            continue;
        end
        
        % Optionally check for pose and seg data; skip if any are missing.
        if exist(pose_train_path, 'file') ~= 2 || exist(pose_val_path, 'file') ~= 2 || ...
           exist(seg_train_path, 'file') ~= 2 || exist(seg_val_path, 'file') ~= 2
            fprintf('Skipping %s because pose or seg data file not found.\n', roi_label);
            continue;
        end
        
        fprintf('Processing %s...\n', roi_label);
        
        %% Load Data
        % Training data
        brain_train = readmatrix(brain_train_path);
        pose_train  = readmatrix(pose_train_path);
        seg_train   = readmatrix(seg_train_path);
        
        % Validation data
        brain_val = readmatrix(brain_val_path);
        pose_val  = readmatrix(pose_val_path);
        seg_val   = readmatrix(seg_val_path);
        
        % Check that dimensions match (assumes data stored as vectors)
        if length(pose_train) ~= length(seg_train) || length(pose_train) ~= length(brain_train)
            fprintf('Skipping %s due to mismatched training data dimensions.\n', roi_label);
            continue;
        end
        if length(pose_val) ~= length(seg_val) || length(pose_val) ~= length(brain_val)
            fprintf('Skipping %s due to mismatched validation data dimensions.\n', roi_label);
            continue;
        end
        
        %% Organize Data for Regression
        X_train = [pose_train, seg_train];
        Y_train = brain_train;
        X_val   = [pose_val, seg_val];
        Y_val   = brain_val;
        
        %% Fit Models and Compute Statistics
        % Full model (both predictors)
        full_model = fitlm(X_train, Y_train);
        Y_pred_full = predict(full_model, X_val);
        r2_full = 1 - sum((Y_val - Y_pred_full).^2) / sum((Y_val - mean(Y_val)).^2);
        rmse_full = sqrt(mean((Y_val - Y_pred_full).^2));
        
        % Pose-only model
        model_pose = fitlm(pose_train, Y_train);
        Y_pred_pose = predict(model_pose, pose_val);
        r2_pose = 1 - sum((Y_val - Y_pred_pose).^2) / sum((Y_val - mean(Y_val)).^2);
        
        % Seg-only model
        model_seg = fitlm(seg_train, Y_train);
        Y_pred_seg = predict(model_seg, seg_val);
        r2_seg = 1 - sum((Y_val - Y_pred_seg).^2) / sum((Y_val - mean(Y_val)).^2);
        
        % Calculate unique contributions and shared variance
        unique_pose = r2_full - r2_seg;   % Unique variance for Pose
        unique_seg  = r2_full - r2_pose;   % Unique variance for Seg
        shared_variance = r2_full - unique_pose - unique_seg;
        
        % Debug output (optional)
        fprintf('%s: Full R2 = %.4f, RMSE = %.4f, Unique Pose = %.4f, Unique Seg = %.4f, Shared = %.4f\n', ...
                roi_label, r2_full, rmse_full, unique_pose, unique_seg, shared_variance);
        
        %% Append Results
        newRow = {roi_label, r2_full, rmse_full, r2_pose, r2_seg, unique_pose, unique_seg, shared_variance};
        results = [results; newRow];
    end
end

%% Save Results to Excel
% Convert results (excluding header row) into a table
resultsTable = cell2table(results(2:end,:), 'VariableNames', results(1,:));

% Define output folder and file name (creates folder if needed)
output_folder = 'D:\ML_project\Variance\var_excel';
if ~exist(output_folder, 'dir')
    mkdir(output_folder);
end
output_file = fullfile(output_folder, sprintf('subject_%d_variance_partitioning.xlsx', subject));

% Write the table to an Excel file
writetable(resultsTable, output_file);
fprintf('Results saved to %s\n', output_file);

toc;  % Display elapsed time
