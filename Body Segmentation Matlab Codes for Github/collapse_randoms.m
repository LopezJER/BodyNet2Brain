% Define subject ID
subject = 8; % Change this value as needed

% Define the file path pattern
input_excel = sprintf('D:\\ML_project\\Variance\\var_excel\\updated_sanitized_allmodels\\subject_%d_variance_partitioning.xlsx', subject);

% Check if the file exists
if ~isfile(input_excel)
    error('File not found: %s', input_excel);
end

% Read input Excel file
data = readtable(input_excel);

% Check if required columns exist
required_columns = {'unique_rpose', 'unique_rseg', 'shared_rpose_rseg', ...
                    'shared_pose_rpose', 'shared_pose_rseg', 'shared_pose_rpose_rseg', ...
                    'shared_seg_rpose', 'shared_seg_rseg', 'shared_seg_rpose_rseg', ...
                    'shared_pose_seg_rpose', 'shared_pose_seg_rseg', 'shared_all', 'r2_rpose', 'r2_rseg'};

missing_columns = setdiff(required_columns, data.Properties.VariableNames);
if ~isempty(missing_columns)
    error('Missing required columns: %s', strjoin(missing_columns, ', '));
end

% Perform column transformations
data.r2_random = data.r2_rpose + data.r2_rseg;
data.unique_random = data.unique_rpose + data.unique_rseg + data.shared_rpose_rseg;
data.shared_pose_random = data.shared_pose_rpose + data.shared_pose_rseg + data.shared_pose_rpose_rseg;
data.shared_seg_random = data.shared_seg_rpose + data.shared_seg_rseg + data.shared_seg_rpose_rseg;
data.shared_all_combined = data.shared_pose_seg_rpose + data.shared_pose_seg_rseg + data.shared_all;

% Remove original columns
data(:, required_columns) = [];

% Define output file path
output_excel = sprintf('D:\\ML_project\\Variance\\var_excel\\updated_sanitized_allmodels\\modified_subject_%d_variance_partitioning.xlsx', subject);

% Write output to a new Excel file
writetable(data, output_excel);

disp(['Modified Excel file saved as: ' output_excel]);
