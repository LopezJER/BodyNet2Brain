% Separated lateralities - variance partitioning results saved to excel
% Updates: handles merged brain data, concatenated DNN data, computes full commonality analysis,
% and now includes individual (full) R² for each model (i.e. unique + all shared portions).

tic;  % Start timer
clear; clc;

%% Configuration
subject = 8;

rois = {'EBA', 'FBA_1', 'FBA_2', 'FFA_1', 'FFA_2', 'hV4', 'MTL', 'OFA', ...
        'V1d', 'V1v', 'V2d', 'V2v', 'V3d', 'V3v'};
laterality_options = {'r', 'l', 'm'};  % 'r', 'l', 'm'

% Base directories:
brain_base = sprintf('D:\\ML_project\\Variance\\brain_results\\subject_%d', subject);
dnn_base   = sprintf('D:\\ML_project\\Variance\\results_final\\subject_%d', subject);

% Prepare output results cell array with header row.
% Headers (22 columns):
% 1: ROI
% 2: Full_R2 (model with all predictors)
% 3: Full_RMSE
% 4-7: r2_pose, r2_seg, r2_rpose, r2_rseg (full R² for each predictor, i.e. unique + shared portions)
% 8-11: unique contributions for pose, seg, rpose, rseg (from commonality coefficients positions 1,2,? see below)
% 12-17: pairwise shared contributions (6 values)
% 18-21: triplet shared contributions (4 values)
% 22: shared_all (quadruple shared)
headers = {'ROI', 'Full_R2', 'Full_RMSE', ...
    'r2_pose', 'r2_seg', 'r2_rpose', 'r2_rseg', ...
    'unique_pose', 'unique_seg', 'unique_rpose', 'unique_rseg', ...
    'shared_pose_seg', 'shared_pose_rpose', 'shared_pose_rseg', 'shared_seg_rpose', 'shared_seg_rseg', 'shared_rpose_rseg', ...
    'shared_pose_seg_rpose', 'shared_pose_seg_rseg', 'shared_pose_rpose_rseg', 'shared_seg_rpose_rseg', ...
    'shared_all'};
results = [headers];

% (Debug lines below are commented out)
% disp('DEBUG: Header size:'); disp(size(results));

%% Loop over each ROI and laterality
for i = 1:length(rois)
    roi = rois{i};
    for j = 1:length(laterality_options)
        laterality_input = laterality_options{j};
        
        % Set strings for DNN and brain data.
        if lower(laterality_input) == 'r'
            dnn_laterality_str = 'right';
            brain_laterality_str = 'right';
        elseif lower(laterality_input) == 'l'
            dnn_laterality_str = 'left';
            brain_laterality_str = 'left';
        elseif lower(laterality_input) == 'm'
            dnn_laterality_str = 'concatenated';
            brain_laterality_str = 'merged';
        else
            error('Laterality must be ''r'', ''l'', or ''m''.');
        end
        
        roi_label = sprintf('%s%s', lower(laterality_input), roi);
        
        % Construct file paths:
        brain_train_path = fullfile(brain_base, 'train', roi, sprintf('%s_%s_data.xlsx', roi, brain_laterality_str));
        brain_val_path   = fullfile(brain_base, 'val', roi, sprintf('%s_%s_data.xlsx', roi, brain_laterality_str));
        
        dnn_types = {'pose', 'seg', 'random_pose', 'random_seg'};
        dnn_train_paths = cell(1,4);
        dnn_val_paths = cell(1,4);
        for d = 1:4
            dnn_train_paths{d} = fullfile(dnn_base, dnn_types{d}, 'train', sprintf('best_rdm_for_%s_%s_data.xlsx', roi, dnn_laterality_str));
            dnn_val_paths{d}   = fullfile(dnn_base, dnn_types{d}, 'val', sprintf('best_rdm_for_%s_%s_data.xlsx', roi, dnn_laterality_str));
        end
        
        if exist(brain_train_path, 'file') ~= 2 || exist(brain_val_path, 'file') ~= 2
            fprintf('Skipping %s because brain data file not found.\n', roi_label);
            continue;
        end
        if any(cellfun(@(x) exist(x, 'file') ~= 2, dnn_train_paths)) || any(cellfun(@(x) exist(x, 'file') ~= 2, dnn_val_paths))
            fprintf('Skipping %s because one or more DNN data files not found.\n', roi_label);
            continue;
        end
        
        fprintf('Processing %s...\n', roi_label);
        
        %% Load Data
        brain_train = readmatrix(brain_train_path);
        brain_val   = readmatrix(brain_val_path);
        
        pose_train  = readmatrix(dnn_train_paths{1});
        pose_val    = readmatrix(dnn_val_paths{1});
        seg_train   = readmatrix(dnn_train_paths{2});
        seg_val     = readmatrix(dnn_val_paths{2});
        rpose_train = readmatrix(dnn_train_paths{3});
        rpose_val   = readmatrix(dnn_val_paths{3});
        rseg_train  = readmatrix(dnn_train_paths{4});
        rseg_val    = readmatrix(dnn_val_paths{4});
        
        if size(pose_train,1) ~= length(brain_train) || size(pose_val,1) ~= length(brain_val)
            fprintf('Skipping %s due to mismatched dimensions.\n', roi_label);
            continue;
        end
        
        Y_train = brain_train;
        Y_val   = brain_val;
        
        % Build full predictor matrices (columns: pose, seg, random_pose, random_seg)
        X_train_all = [pose_train, seg_train, rpose_train, rseg_train];
        X_val_all   = [pose_val, seg_val, rpose_val, rseg_val];
        
        %% Full Model Statistics
        full_model = fitlm(X_train_all, Y_train);
        Y_pred_full = predict(full_model, X_val_all);
        r2_full = 1 - sum((Y_val - Y_pred_full).^2) / sum((Y_val - mean(Y_val)).^2);
        rmse_full = sqrt(mean((Y_val - Y_pred_full).^2));
        
        %% Individual Models for Full R² (each predictor separately)
        % (We will recompute full predictor R² from commonality coefficients later.)
        % For now, we compute them separately to compare (but these might not be used in output).
        % r2_pose_ind = 1 - sum((Y_val - predict(fitlm(pose_train, Y_train), pose_val)).^2) / sum((Y_val - mean(Y_val)).^2);
        % ... (but we'll calculate full predictor R² from commonality)
        
        %% Compute R² for Every Nonempty Subset of Predictors (15 subsets)
        numPredictors = 4;
        nSubsets = 2^numPredictors - 1;  
        R2_sub = zeros(1, nSubsets);
        for m = 1:nSubsets
            indices = find(bitget(m, 1:numPredictors));
            X_train_sub = X_train_all(:, indices);
            X_val_sub = X_val_all(:, indices);
            model_sub = fitlm(X_train_sub, Y_train);
            Y_pred_sub = predict(model_sub, X_val_sub);
            R2_sub(m) = 1 - sum((Y_val - Y_pred_sub).^2) / sum((Y_val - mean(Y_val)).^2);
        end
        
        %% Compute Commonality Coefficients using Moebius Inversion
        c = zeros(1, nSubsets);
        for m = 1:nSubsets
            A_count = sum(bitget(m, 1:numPredictors));
            ntemp = m;
            c(m) = 0;
            while true
                B_count = sum(bitget(ntemp, 1:numPredictors));
                c(m) = c(m) + (-1)^(A_count - B_count) * R2_sub(ntemp);
                if ntemp == 0, break; end
                ntemp = bitand(ntemp-1, m);
                if ntemp == 0, break; end
            end
        end
        
        %% Reorder Commonality Coefficients to Conventional Order
        % Our chosen conventional order (15 coefficients) is:
        % Unique: {1} (index 1), {2} (2), {3} (4), {4} (8)
        % Pairwise: {1,2} (3), {1,3} (5), {1,4} (9), {2,3} (6), {2,4} (10), {3,4} (12)
        % Triplets: {1,2,3} (7), {1,2,4} (11), {1,3,4} (13), {2,3,4} (14)
        % Quadruple: {1,2,3,4} (15)  -> renamed "shared_all"
        orderIdx = [1, 2, 4, 8, 3, 5, 9, 6, 10, 12, 7, 11, 13, 14, 15];
        commonality = c(orderIdx);
        
        %% Compute Full Predictor R² from Commonality Coefficients
        % For predictor 1 (pose): subsets including predictor 1 are those with original indices: 1,3,5,7,9,11,13,15.
        full_r2_pose = commonality(1) + commonality(5) + commonality(6) + commonality(11) + commonality(7) + commonality(12) + commonality(13) + commonality(15);
        % For predictor 2 (seg): subsets including predictor 2: indices 2,3,6,7,10,11,14,15.
        full_r2_seg = commonality(2) + commonality(5) + commonality(8) + commonality(11) + commonality(9) + commonality(12) + commonality(14) + commonality(15);
        % For predictor 3 (rpose): subsets including predictor 3: indices 4,5,6,7,12,13,14,15.
        full_r2_rpose = commonality(3) + commonality(6) + commonality(8) + commonality(11) + commonality(10) + commonality(13) + commonality(14) + commonality(15);
        % For predictor 4 (rseg): subsets including predictor 4: indices 8,9,10,11,12,13,14,15.
        full_r2_rseg = commonality(4) + commonality(7) + commonality(9) + commonality(10) + commonality(12) + commonality(13) + commonality(14) + commonality(15);
        
        %% Prepare Result Row
        newRow = {roi_label, r2_full, rmse_full, ...
            full_r2_pose, full_r2_seg, full_r2_rpose, full_r2_rseg, ... % Full R² for each predictor
            commonality(1), commonality(2), commonality(3), commonality(4), ...  % unique: pose, seg, rpose, rseg
            commonality(5), commonality(6), commonality(7), commonality(8), commonality(9), commonality(10), ... % pairwise
            commonality(11), commonality(12), commonality(13), commonality(14), ...
            commonality(15)};  % shared_all
        
        % (Optional debug prints are commented out.)
        %{
        disp('DEBUG: Size of results array:'); disp(size(results));
        disp('DEBUG: New row size:'); disp(size(newRow));
        fprintf('DEBUG: Unique random_pose = %.4f, Unique random_seg = %.4f\n', commonality(3), commonality(4));
        %}
        
        if length(newRow) ~= length(headers)
            error('Dimension mismatch: newRow has %d columns but header has %d.', length(newRow), length(headers));
        end
        
        results = [results; newRow];
    end
end

%% Save Results to Excel
resultsTable = cell2table(results(2:end,:), 'VariableNames', results(1,:));

output_folder = 'D:\ML_project\Variance\var_excel';
if ~exist(output_folder, 'dir')
    mkdir(output_folder);
end
output_file = fullfile(output_folder, sprintf('subject_%d_variance_partitioning.xlsx', subject));

writetable(resultsTable, output_file);
fprintf('Results saved to %s\n', output_file);

toc;  % Display elapsed time