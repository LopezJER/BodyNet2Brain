% Updated Variance Partitioning Code - 4 DNN models (Pose, Seg, Depth, Normals), merged brain ROIs only

tic;  % Start timer
clear; clc;

%% Configuration
subject = 8;

rois = {'EBA', 'FBA_1', 'FBA_2', 'FFA_1', 'FFA_2', 'hV4', 'MTL', 'OFA', ...
        'V1d', 'V1v', 'V2d', 'V2v', 'V3d', 'V3v'};

% Updated base directories
brain_base = sprintf('D:\\ML_project\\Variance\\new_brain_results\\subject_%d', subject);
dnn_base   = sprintf('D:\\ML_project\\Variance\\sapiens_results2\\subject_%d', subject);

% Headers
headers = {'ROI', 'Full_R2', 'Full_RMSE', ...
    'r2_pose', 'r2_seg', 'r2_depth', 'r2_normals', ...
    'unique_pose', 'unique_seg', 'unique_depth', 'unique_normals', ...
    'shared_pose_seg', 'shared_pose_depth', 'shared_pose_normals', 'shared_seg_depth', 'shared_seg_normals', 'shared_depth_normals', ...
    'shared_pose_seg_depth', 'shared_pose_seg_normals', 'shared_pose_depth_normals', 'shared_seg_depth_normals', ...
    'shared_all'};
results = [headers];

%% Loop over each ROI
for i = 1:length(rois)
    roi = rois{i};
    roi_label = roi;

    % File paths
    brain_train_path = fullfile(brain_base, 'train', roi, sprintf('%s_merged_data.xlsx', roi));
    brain_val_path   = fullfile(brain_base, 'val', roi, sprintf('%s_merged_data.xlsx', roi));

    dnn_types = {'pose', 'seg', 'depth', 'normal'};
    dnn_train_paths = cell(1,4);
    dnn_val_paths = cell(1,4);
    for d = 1:4
        dnn_train_paths{d} = fullfile(dnn_base, dnn_types{d}, 'train', sprintf('best_rdm_for_%s_concatenated_data.xlsx', roi));
        dnn_val_paths{d}   = fullfile(dnn_base, dnn_types{d}, 'val', sprintf('best_rdm_for_%s_concatenated_data.xlsx', roi));
    end

    % Skip if any file is missing
    if exist(brain_train_path, 'file') ~= 2 || exist(brain_val_path, 'file') ~= 2 || ...
       any(cellfun(@(x) exist(x, 'file') ~= 2, dnn_train_paths)) || ...
       any(cellfun(@(x) exist(x, 'file') ~= 2, dnn_val_paths))
        fprintf('Skipping %s due to missing files.\n', roi_label);
        continue;
    end

    fprintf('Processing %s...\n', roi_label);

    %% Load data
    brain_train = readmatrix(brain_train_path);
    brain_val   = readmatrix(brain_val_path);

    dnn_train = cellfun(@readmatrix, dnn_train_paths, 'UniformOutput', false);
    dnn_val   = cellfun(@readmatrix, dnn_val_paths, 'UniformOutput', false);

    % Sanity check
    if any(cellfun(@(x) size(x,1), dnn_train) ~= length(brain_train)) || ...
       any(cellfun(@(x) size(x,1), dnn_val) ~= length(brain_val))
        fprintf('Skipping %s due to size mismatch.\n', roi_label);
        continue;
    end

    %% Full model
    X_train_all = horzcat(dnn_train{:});
    X_val_all   = horzcat(dnn_val{:});
    Y_train = brain_train;
    Y_val   = brain_val;

    full_model = fitlm(X_train_all, Y_train);
    Y_pred_full = predict(full_model, X_val_all);
    r2_full = 1 - sum((Y_val - Y_pred_full).^2) / sum((Y_val - mean(Y_val)).^2);
    rmse_full = sqrt(mean((Y_val - Y_pred_full).^2));

    %% R² for all 15 predictor subsets
    numPredictors = 4;
    nSubsets = 2^numPredictors - 1;
    R2_sub = zeros(1, nSubsets);
    for m = 1:nSubsets
        idx = find(bitget(m, 1:numPredictors));
        X_train_sub = X_train_all(:, idx);
        X_val_sub = X_val_all(:, idx);
        model = fitlm(X_train_sub, Y_train);
        Y_pred = predict(model, X_val_sub);
        R2_sub(m) = 1 - sum((Y_val - Y_pred).^2) / sum((Y_val - mean(Y_val)).^2);
    end

    %% Moebius inversion (commonality coefficients)
    c = zeros(1, nSubsets);
    for m = 1:nSubsets
        A_count = sum(bitget(m, 1:numPredictors));
        ntemp = m;
        while true
            B_count = sum(bitget(ntemp, 1:numPredictors));
            c(m) = c(m) + (-1)^(A_count - B_count) * R2_sub(ntemp);
            if ntemp == 0, break; end
            ntemp = bitand(ntemp-1, m);
            if ntemp == 0, break; end
        end
    end

    %% Rearrange order
    orderIdx = [1, 2, 4, 8, 3, 5, 9, 6, 10, 12, 7, 11, 13, 14, 15];
    commonality = c(orderIdx);

    %% Full R² per predictor
    full_r2 = zeros(1, 4);
    full_r2(1) = sum(commonality([1,5,6,11,7,12,13,15]));  % pose
    full_r2(2) = sum(commonality([2,5,8,11,9,12,14,15]));  % seg
    full_r2(3) = sum(commonality([3,6,8,11,10,13,14,15])); % depth
    full_r2(4) = sum(commonality([4,7,9,10,12,13,14,15])); % normals

    %% Create result row
    newRow = {roi_label, r2_full, rmse_full, ...
        full_r2(1), full_r2(2), full_r2(3), full_r2(4), ...
        commonality(1), commonality(2), commonality(3), commonality(4), ...
        commonality(5), commonality(6), commonality(7), commonality(8), commonality(9), commonality(10), ...
        commonality(11), commonality(12), commonality(13), commonality(14), ...
        commonality(15)};
    results = [results; newRow];
end

%% Save to Excel
resultsTable = cell2table(results(2:end,:), 'VariableNames', results(1,:));
output_folder = 'D:\ML_project\Variance\var_excel';
if ~exist(output_folder, 'dir')
    mkdir(output_folder);
end
output_file = fullfile(output_folder, sprintf('subject_%d_variance_partitioning_updated.xlsx', subject));
writetable(resultsTable, output_file);
fprintf('Results saved to %s\n', output_file);

toc;  % End timer
