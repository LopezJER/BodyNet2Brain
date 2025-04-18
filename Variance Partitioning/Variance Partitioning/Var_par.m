%% Variance Partitioning with Train/Validation Splits
% This script fits regression models on training data and evaluates R² on a separate validation set.
% It decomposes the total explained variance into unique contributions of each DNN model and their shared contribution.
% It then produces a Venn diagram and a bar graph.
%
% Training & Validation:
% - Models are fit on the training set.
% - R² is computed on the validation set.
% Note: Negative unique contributions may occur when predictors are collinear. For plotting,
%       such values are forced to a small positive number so that both circles appear.

clear; clc;
tic;  % Start timer

%% Configuration: Names, File Paths, and Colors
dnn_model_names = {'Pose Estimation', 'Body Segmentation'};  
% User inputs:
laterality_input = 'r';  % Enter 'r' for right or 'l' for left
brain_region = 'V2v';    % Name of the brain region

% Convert laterality letter to full string:
if lower(laterality_input) == 'r'
    laterality_str = 'right';
elseif lower(laterality_input) == 'l'
    laterality_str = 'left';
else
    error('Laterality must be ''r'' or ''l''.');
end

% Create a descriptive brain region name (e.g., 'Right EBA' or 'Left EBA'):
laterality_cap = [upper(laterality_str(1)) laterality_str(2:end)];
brain_region_name = sprintf('%s %s', laterality_cap, brain_region);

% Base directories:
base_dir = 'D:\ML_project\Variance';
dnn_base = fullfile(base_dir, 'sapiens_results', 'subject_1');
brain_base = fullfile(base_dir, 'brain_results', 'subject_1');

% Dynamically construct file paths:
% Training file paths:
dnn1_train_path = fullfile(dnn_base, 'pose', 'train', sprintf('best_rdm_for_%s_%s_data.xlsx', brain_region, laterality_str));
dnn2_train_path = fullfile(dnn_base, 'seg', 'train', sprintf('best_rdm_for_%s_%s_data.xlsx', brain_region, laterality_str));
brain_train_path = fullfile(brain_base, 'train', brain_region, sprintf('%s_%s_data.xlsx', brain_region, laterality_str));

% Validation file paths:
dnn1_val_path = fullfile(dnn_base, 'pose', 'val', sprintf('best_rdm_for_%s_%s_data.xlsx', brain_region, laterality_str));
dnn2_val_path = fullfile(dnn_base, 'seg', 'val', sprintf('best_rdm_for_%s_%s_data.xlsx', brain_region, laterality_str));
brain_val_path = fullfile(brain_base, 'val', brain_region, sprintf('%s_%s_data.xlsx', brain_region, laterality_str));


% Color definitions:
% For Venn diagrams and bar graph.
% For Pose Estimation (right circle): brighter green
brightGreen = [0, 0.8, 0];
% For Body Segmentation (left circle): strong orange
strongOrange = [1, 0.6, 0];
% Neon yellow for shared region
neonYellow = [1, 1, 0];

%% Step 1: Load Training Data
fprintf('Loading training data...\n');
dnn1_train = readmatrix(dnn1_train_path);
dnn2_train = readmatrix(dnn2_train_path);
brain_train = readmatrix(brain_train_path);

%% Step 2: Load Validation Data
fprintf('Loading validation data...\n');
dnn1_val = readmatrix(dnn1_val_path);
dnn2_val = readmatrix(dnn2_val_path);
brain_val = readmatrix(brain_val_path);

% Check dimensions (assuming data is stored as vectors):
if length(dnn1_train) ~= length(dnn2_train) || length(dnn1_train) ~= length(brain_train)
    error('Training data dimensions do not match.');
end
if length(dnn1_val) ~= length(dnn2_val) || length(dnn1_val) ~= length(brain_val)
    error('Validation data dimensions do not match.');
end

%% Step 3: Organize Data for Regression
X_train = [dnn1_train, dnn2_train];
Y_train = brain_train;
X_val = [dnn1_val, dnn2_val];
Y_val = brain_val;

%% Step 4: Train Models and Evaluate on Validation Set
fprintf('Fitting full model on training data...\n');
full_model = fitlm(X_train, Y_train);
Y_pred_full = predict(full_model, X_val);
r2_full = 1 - sum((Y_val - Y_pred_full).^2) / sum((Y_val - mean(Y_val)).^2);
rmse_full = sqrt(mean((Y_val - Y_pred_full).^2));

fprintf('Fitting model with %s only on training data...\n', dnn_model_names{1});
model_dnn1 = fitlm(dnn1_train, Y_train);
Y_pred_dnn1 = predict(model_dnn1, dnn1_val);
r2_dnn1 = 1 - sum((Y_val - Y_pred_dnn1).^2) / sum((Y_val - mean(Y_val)).^2);

fprintf('Fitting model with %s only on training data...\n', dnn_model_names{2});
model_dnn2 = fitlm(dnn2_train, Y_train);
Y_pred_dnn2 = predict(model_dnn2, dnn2_val);
r2_dnn2 = 1 - sum((Y_val - Y_pred_dnn2).^2) / sum((Y_val - mean(Y_val)).^2);

%% Step 5: Calculate Variance Contributions
unique_dnn1 = r2_full - r2_dnn2;  % Pose Estimation (right circle)
unique_dnn2 = r2_full - r2_dnn1;  % Body Segmentation (left circle)
shared_variance = r2_full - unique_dnn1 - unique_dnn2;

% Debug prints:
fprintf('Debug: unique_dnn1 = %.4f, unique_dnn2 = %.4f, shared_variance = %.4f\n', unique_dnn1, unique_dnn2, shared_variance);

%% Step 6: Display Results
fprintf('\n=== Variance Partitioning Results (Validation) for %s ===\n', brain_region_name);
fprintf('Full Model R²: %.4f\n', r2_full);
fprintf('Validation RMSE of Full Model: %.4f\n', rmse_full);
fprintf('Unique Contribution of %s: %.4f\n', dnn_model_names{1}, unique_dnn1);
fprintf('Unique Contribution of %s: %.4f\n', dnn_model_names{2}, unique_dnn2);
fprintf('Shared Contribution: %.4f\n', shared_variance);

%% Step 7: Plot Venn Diagram
% We want circle areas to be proportional to the variance they represent.
% For each circle, the total area should equal (unique + shared).
% Define:
%   Area_right = unique_dnn1 + shared_variance  (Pose Estimation)
%   Area_left  = unique_dnn2 + shared_variance  (Body Segmentation)
%
% We'll use a scaling factor so the circles are large enough for display.
scaling = 100;  
area_right = unique_dnn1 + shared_variance;
area_left  = unique_dnn2 + shared_variance;
% Scaled areas:
A_right = area_right * scaling^2;
A_left  = area_left  * scaling^2;
A_shared = shared_variance * scaling^2;  % desired intersection area

% Compute radii (r = sqrt(area/pi))
r_right = sqrt(A_right/pi);
r_left  = sqrt(A_left/pi);

% Define an inline function for the intersection area of two circles:
clamp = @(x) min(max(x, -1), 1);  % ensures values passed to acos are in [-1,1]
intersection_area = @(d, r1, r2) ...
    r1^2 * acos(clamp((d^2 + r1^2 - r2^2)/(2*d*r1))) + ...
    r2^2 * acos(clamp((d^2 + r2^2 - r1^2)/(2*d*r2))) - ...
    0.5 * sqrt( max(0,(-d+r1+r2)*(d+r1-r2)*(d-r1+r2)*(d+r1+r2) ) );

% Set bounds for the center distance:
lb = max(abs(r_left - r_right), eps);
ub = r_left + r_right - eps;

% Compute maximum possible intersection area at d = lb:
max_int = intersection_area(lb, r_left, r_right);
if A_shared > max_int
    % Cap A_shared if it exceeds the maximum possible intersection area.
    A_shared = max_int;
end

% Define the function for which we solve f(d) = 0.
f = @(d) intersection_area(d, r_left, r_right) - A_shared;

% Check the signs at the bounds:
f_lb = f(lb);
f_ub = f(ub);
if abs(f_lb) < 1e-8
    d = lb;
elseif abs(f_ub) < 1e-8
    d = ub;
elseif f_lb * f_ub > 0
    % If no sign change is found, default to the lower bound (or handle as needed)
    warning('No sign change found in fzero interval; defaulting d to lb.');
    d = lb;
else
    % Solve for d using fzero.
    d = fzero(f, [lb, ub]);
end

% Set centers symmetrically:
center_left = [-d/2, 0];
center_right = [d/2, 0];

% Generate circle coordinates using theta; remove duplicate endpoint.
theta = linspace(0, 2*pi, 200);
theta(end) = [];

% Compute coordinates for each circle:
x_left = center_left(1) + r_left * cos(theta);
y_left = center_left(2) + r_left * sin(theta);
x_right = center_right(1) + r_right * cos(theta);
y_right = center_right(2) + r_right * sin(theta);

% Create polyshape objects.
poly_left = polyshape(x_left, y_left);
poly_right = polyshape(x_right, y_right);
poly_int = intersect(poly_left, poly_right);

% Now plot the circles.
figure;
hold on;
% Plot left circle (Body Segmentation) with full opacity.
plot(poly_left, 'FaceColor', strongOrange, 'FaceAlpha', 1, 'EdgeColor', 'none');
% Plot right circle (Pose Estimation) with full opacity.
plot(poly_right, 'FaceColor', brightGreen, 'FaceAlpha', 1, 'EdgeColor', 'none');
% Plot the intersection (shared variance) on top.
if poly_int.NumRegions > 0
    plot(poly_int, 'FaceColor', neonYellow, 'FaceAlpha', 1, 'EdgeColor', 'none');
end

% Remove any text labels that might have been added automatically.
delete(findobj(gca, 'Type', 'text'));

% Add a title and legend.
title(sprintf('Variance Partitioning Venn Diagram for %s', brain_region_name));
legend({sprintf('Unique %s (%.4f)', dnn_model_names{1}, unique_dnn1), ...
        sprintf('Unique %s (%.4f)', dnn_model_names{2}, unique_dnn2), ...
        sprintf('Shared (%.4f)', shared_variance)}, 'Location', 'BestOutside');
axis equal;
hold off;

%% Step 8: Plot Bar Graph of Unique Variances
figure;
% Use a categorical variable for the brain region.
c = categorical({brain_region_name});
data = [unique_dnn1, unique_dnn2];  % Order: first is Pose Estimation, second is Body Segmentation
% Use grouped bar to get separate objects:
b = bar(c, data, 'grouped');
if numel(b) > 1
    for i = 1:length(b)
        b(i).FaceColor = 'flat';
    end
    % Assign brightGreen to Pose Estimation and strongOrange to Body Segmentation.
    b(1).CData = brightGreen;
    b(2).CData = strongOrange;
else
    b.FaceColor = 'flat';
    b.CData = [brightGreen; strongOrange];
end
ylabel('Unique Variance (R²)');
title(sprintf('Unique Variance Explained by %s and %s in %s', dnn_model_names{1}, dnn_model_names{2}, brain_region_name));
legend({dnn_model_names{1}, dnn_model_names{2}}, 'Location', 'Best');
ax = gca;
ax.YAxis.Exponent = 0;
ytickformat('%.4f');

%% Step 9: Display Detailed Model Information
fprintf('\n=== Model Details ===\n');
disp(full_model);
disp(model_dnn1);
disp(model_dnn2);

%% Step 10: End Timer and Display Execution Time
elapsed_time = toc;
fprintf('\nTotal Execution Time: %.2f seconds\n', elapsed_time);