tic;  % Start timing
clear; clc;

%% Configuration
roiPool = {'EBA', 'FBA_1', 'FBA_2', 'FFA_1', 'FFA_2', 'hV4', 'MTL', 'OFA', ...
           'V1d', 'V1v', 'V2d', 'V2v', 'V3d', 'V3v'};
laterality_options = {'r', 'l', 'm'};  % 'r' for right, 'l' for left, 'm' for merged

input_excel = 'D:\ML_project\Variance\var_excel\updated_sanitized_allmodels\collapsed_randoms\modified_subject_1_variance_partitioning.xlsx';

% User Input
roi_input = input('Enter ROI name (e.g., rEBA, lFFA_1, mV1d): ', 's');
laterality = roi_input(1); % Extract first character for laterality
roi = roi_input(2:end); % Extract ROI name
if laterality == 'r'
    laterality_label = 'Right';
elseif laterality == 'l'
    laterality_label = 'Left';
else
    laterality_label = 'Merged';
end 

display_only = input('Display only? (1 for Yes, 0 for No): ');
save_option = input('Save image? (1 for Yes, 0 for No): ');

if ~ismember(roi, roiPool) || ~ismember(laterality, laterality_options)
    error('Invalid ROI or laterality. Please enter a valid combination (e.g., rEBA, lFFA_1, mV1d).');
end

roiLabel = sprintf('%s%s', lower(laterality), roi);

% Define colors
brightGreen = [0, 0.8, 0];       % Unique Pose Estimation (green)
strongOrange = [1, 0.6, 0];      % Unique Body Segmentation (orange)
neonYellow = [1, 1, 0];          % Shared Variance (yellow)
gray = [0.6, 0.6, 0.6];          % Shared Pose Random (gray)

%% Read the Excel file
dataTable = readtable(input_excel);
idx = strcmp(dataTable.ROI, roiLabel);
if ~any(idx)
    error('ROI not found in the dataset.');
end

% Extract required data
real_unique_pose = dataTable.unique_pose(idx);
real_unique_seg = dataTable.unique_seg(idx);
real_shared_variance = dataTable.shared_pose_seg(idx);
real_shared_pose_random = dataTable.shared_pose_random(idx);

% Ensure non-positive values are treated as zero for visualization
unique_pose = max(0, real_unique_pose - real_shared_pose_random);
unique_seg = max(0, real_unique_seg);
shared_variance = max(0, real_shared_variance);
shared_pose_random = max(0, real_shared_pose_random);

% Calculate proportional areas
scaling = 10000;
total_variance = unique_pose + unique_seg + shared_variance;
A_pose = (unique_pose / total_variance) * scaling;
A_seg  = (unique_seg / total_variance) * scaling;
A_shared = (shared_variance / total_variance) * scaling;
A_shared_random = (shared_pose_random / total_variance) * scaling;

% Compute radii ensuring correct proportions
r_pose = sqrt(A_pose/pi);
r_seg  = sqrt(A_seg/pi);
r_shared_random = sqrt(A_shared_random/pi);

% Determine spacing for circles
if shared_variance > 0
    d = r_pose + r_seg - sqrt(A_shared/pi);  % Overlapping distance based on shared variance
else
    d = r_pose + r_seg + 0.1 * max(r_pose, r_seg);  % Ensure full separation if no shared variance
end

% Create Venn Diagram
figVenn = figure('Visible', 'on', 'Position', [100, 100, 1200, 800]);
hold on;

% Main Circles
theta = linspace(0, 2*pi, 200);
if unique_pose > 0
    x_left = -d/2 + r_pose * cos(theta);
    y_left = r_pose * sin(theta);
    poly_left = polyshape(x_left, y_left);
    plot(poly_left, 'FaceColor', brightGreen, 'FaceAlpha', 1, 'EdgeColor', 'none');
end

if unique_seg > 0
    x_right = d/2 + r_seg * cos(theta);
    y_right = r_seg * sin(theta);
    poly_right = polyshape(x_right, y_right);
    plot(poly_right, 'FaceColor', strongOrange, 'FaceAlpha', 1, 'EdgeColor', 'none');
end

if shared_variance > 0
    poly_int = intersect(poly_left, poly_right);
    if poly_int.NumRegions > 0
        plot(poly_int, 'FaceColor', neonYellow, 'FaceAlpha', 1, 'EdgeColor', 'none');
    end
end

if shared_pose_random > 0
    x_inner = -d/2 + r_shared_random * cos(theta);
    y_inner = r_shared_random * sin(theta);
    poly_inner = polyshape(x_inner, y_inner);
    plot(poly_inner, 'FaceColor', gray, 'FaceAlpha', 1, 'EdgeColor', 'none');
end

% Labels and Legend
legendEntries = {
    sprintf('Pose Estimation (%.4f)', real_unique_pose),
    sprintf('Body Segmentation (%.4f)', real_unique_seg),
    sprintf('Pose Est. and Body seg. (%.4f)', real_shared_variance),
    sprintf('Random models & Pose Est. (%.4f)', real_shared_pose_random)
};
legendColors = [brightGreen; strongOrange; neonYellow; gray];
legendObjects = [];
for i = 1:length(legendEntries)
    legendObjects(i) = plot(nan, nan, 's', 'MarkerSize', 10, 'MarkerFaceColor', legendColors(i, :), 'MarkerEdgeColor', 'none');
end
legend(legendObjects, legendEntries, 'Location', 'northeast', 'FontSize', 8);

title(sprintf('Variance Partitioning in %s %s of Subject 1 Across Predictors', laterality_label, roi), 'Interpreter', 'none', 'FontSize', 14);
axis equal; axis off;
hold off;

% Save the figure if the user chooses
if save_option
    output_folder = 'D:\ML_project\Variance\var_graphs\custom';
    if ~exist(output_folder, 'dir')
        mkdir(output_folder);
    end
    timestamp = datestr(now, 'yyyy-mm-dd_HH-MM-SS');
    vennFile = fullfile(output_folder, sprintf('%s_%s_Venn_%s.png', laterality_label, roi, timestamp));
    saveas(figVenn, vennFile);
    fprintf('Saved Venn diagram for %s %s to %s\n', laterality_label, roi, vennFile);
end

toc;
