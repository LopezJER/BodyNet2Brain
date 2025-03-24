tic;  % Start timing
clear; clc;

%% Configuration
% subject = 8;
% input_excel = sprintf('D:\\ML_project\\Variance\\var_excel\\updated_sanitized_allmodels\\subject_%d_variance_partitioning.xlsx', subject);
%input_excel = sprintf('D:\\ML_project\\Variance\\var_excel\\updated_sanitized_allmodels\\Aggregated_ROI_Results.xlsx');


roiPool = {'EBA', 'FBA_1', 'FBA_2', 'FFA_1', 'FFA_2', 'hV4', 'MTL', 'OFA', ...
           'V1d', 'V1v', 'V2d', 'V2v', 'V3d', 'V3v'};
laterality_options = {'r', 'l', 'm'};  % 'r' for right, 'l' for left, 'm' for merged

% Define colors for the figures
brightGreen = [0, 0.8, 0];       % Unique Pose Estimation (green)
strongOrange = [1, 0.6, 0];      % Unique Body Segmentation (orange)
neonYellow = [1, 1, 0];          % Shared Variance (for Venn diagram)
blackColor   = [0, 0, 0];        % Full Model (black)
redColor = [1, 0, 0];            % Random Pose Estimation (red)
greyColor = [0.5, 0.5, 0.5];     % Random Body Segmentation (grey)

% Define fixed y-axis ticks and limits for all bar graphs
uniformYTicks = -0.1:0.02:0.2;   % ticks from -0.1 to 0.2 every 0.02

% Define output folder and create a subject-specific subfolder
base_output_folder = 'D:\ML_project\Variance\var_graphs\new_sanitized';
% subject_folder = fullfile(base_output_folder, sprintf('subject_%d', subject));
subject_folder = fullfile(base_output_folder, sprintf('Average'));

if ~exist(subject_folder, 'dir')
    mkdir(subject_folder);
    fprintf('DEBUG: Created subject folder: %s\n', subject_folder);
end

%% Read the Excel file
fprintf('DEBUG: Reading Excel file: %s\n', input_excel);
dataTable = readtable(input_excel);
fprintf('DEBUG: Excel file loaded. Found %d rows.\n', height(dataTable));

%% Loop over ROI pool and laterality options
for r = 1:length(roiPool)
    roi = roiPool{r};
    for lat = 1:length(laterality_options)
        laterality = laterality_options{lat};  % 'r', 'l', or 'm'
        % For file naming (and Excel lookup), use the laterality prefix:
        roiLabel = sprintf('%s%s', lower(laterality), roi);
        % For the graph title, if merged, omit any laterality word.
        if strcmpi(laterality, 'm')
            roiTitle = roi;
        elseif strcmpi(laterality, 'r')
            roiTitle = sprintf('Right %s', roi);
        elseif strcmpi(laterality, 'l')
            roiTitle = sprintf('Left %s', roi);
        end
        
        fprintf('DEBUG: Processing ROI: %s\n', roiTitle);
        
        % Create a subfolder for the current ROI within the subject folder
        roi_subfolder = fullfile(subject_folder, roiLabel);
        if ~exist(roi_subfolder, 'dir')
            mkdir(roi_subfolder);
            fprintf('DEBUG: Created ROI subfolder: %s\n', roi_subfolder);
        end
        
        % Find the row corresponding to this ROI label in the Excel table.
        idx = strcmp(dataTable.ROI, roiLabel);
        if ~any(idx)
            fprintf('DEBUG: No data found for %s. Skipping.\n', roiTitle);
            continue;
        end
        
%         % Extract raw values from the table (using the original column names)
%         unique_pose = dataTable.unique_pose(idx);
%         unique_seg = dataTable.unique_seg(idx);
%         % Here, shared variance is read from column "shared_pose_seg"
%         shared_variance = dataTable.shared_pose_seg(idx);
%         full_R2 = dataTable.Full_R2(idx);
%         unique_rpose = dataTable.unique_rpose(idx);  % new column
%         unique_rseg = dataTable.unique_rseg(idx);      % new column
        % Extract raw values from the table (using the original column names)
        unique_pose = dataTable.Avg_unique_pose(idx);
        unique_seg = dataTable.Avg_unique_seg(idx);
        % Here, shared variance is read from column "shared_pose_seg"
        shared_variance = dataTable.Avg_shared_pose_seg(idx);
        full_R2 = dataTable.Avg_Full_R2(idx);
        unique_rpose = dataTable.Avg_unique_rpose(idx);  % new column
        unique_rseg = dataTable.Avg_unique_rseg(idx);      % new column
        
        % In case multiple rows exist, take the first one.
        unique_pose = unique_pose(1);
        unique_seg = unique_seg(1);
        shared_variance = shared_variance(1);
        full_R2 = full_R2(1);
        unique_rpose = unique_rpose(1);
        unique_rseg = unique_rseg(1);
        
        fprintf('DEBUG: Raw values for %s: Unique Pose = %.4f, Unique Seg = %.4f, Shared Variance = %.4f, Full R^2 = %.4f, Random Pose = %.4f, Random Seg = %.4f\n', ...
            roiTitle, unique_pose, unique_seg, shared_variance, full_R2, unique_rpose, unique_rseg);
        
        % For the Venn diagram, clip negative values to zero for display.
        venn_unique_pose = max(0, unique_pose);
        venn_unique_seg = max(0, unique_seg);
        venn_shared_variance = max(0, shared_variance);
        
        %% Create Venn Diagram (Version 1: With Numbers in Legend)
        % Calculate circle areas for Pose and Segmentation using (value + shared)
        scaling = 100;  % scaling factor for display
        area_pose = venn_unique_pose + venn_shared_variance;
        area_seg  = venn_unique_seg + venn_shared_variance;
        A_pose = area_pose * scaling^2;
        A_seg  = area_seg  * scaling^2;
        A_shared = venn_shared_variance * scaling^2;  % desired intersection area
        
        % Compute radii for each circle.
        r_pose = sqrt(A_pose/pi);   % for Pose Estimation (green)
        r_seg  = sqrt(A_seg/pi);    % for Body Segmentation (orange)
        fprintf('DEBUG: Venn diagram radii for %s: r_pose (green) = %.4f, r_seg (orange) = %.4f\n', roiTitle, r_pose, r_seg);
        
        % Define function for intersection area.
        clamp = @(x) min(max(x, -1), 1);
        intersection_area = @(d, r1, r2) ...
            r1^2 * acos(clamp((d^2 + r1^2 - r2^2)/(2*d*r1))) + ...
            r2^2 * acos(clamp((d^2 + r2^2 - r1^2)/(2*d*r2))) - ...
            0.5 * sqrt(max(0, (-d+r1+r2)*(d+r1-r2)*(d-r1+r2)*(d+r1+r2)));
        
        % Determine center distance d
        lb = max(abs(r_pose - r_seg), eps);
        ub = r_pose + r_seg - eps;
        max_int = intersection_area(lb, r_pose, r_seg);
        if A_shared > max_int
            fprintf('DEBUG: Adjusting A_shared for %s from %.4f to max_int %.4f\n', roiTitle, A_shared, max_int);
            A_shared = max_int;
        end
        
        if venn_shared_variance <= 0
            d = r_pose + r_seg;  % No overlap.
        else
            f = @(d) intersection_area(d, r_pose, r_seg) - A_shared;
            if abs(f(lb)) < 1e-8
                d = lb;
            elseif abs(f(ub)) < 1e-8
                d = ub;
            elseif f(lb)*f(ub) > 0
                warning('DEBUG: No sign change found for %s; defaulting d to lb.', roiTitle);
                d = lb;
            else
                d = fzero(f, [lb, ub]);
            end
        end
        
        fprintf('DEBUG: Computed center distance d for %s: %.4f\n', roiTitle, d);
        
        % Position circles so that the Pose Estimation (green) circle is on the left.
        center_left = [-d/2, 0];
        center_right = [d/2, 0];
        theta = linspace(0, 2*pi, 200); theta(end) = [];
        x_left = center_left(1) + r_pose * cos(theta);
        y_left = center_left(2) + r_pose * sin(theta);
        x_right = center_right(1) + r_seg * cos(theta);
        y_right = center_right(2) + r_seg * sin(theta);
        
        poly_left = polyshape(x_left, y_left);
        poly_right = polyshape(x_right, y_right);
        poly_int = intersect(poly_left, poly_right);
        
        % Create a larger figure for Venn diagram and reposition the legend.
        figVennNum = figure('Visible', 'off', 'Position', [100, 100, 1200, 800]);
        hold on;
        p1 = plot(poly_left, 'FaceColor', brightGreen, 'FaceAlpha', 1, 'EdgeColor', 'none');
        p2 = plot(poly_right, 'FaceColor', strongOrange, 'FaceAlpha', 1, 'EdgeColor', 'none');
        if poly_int.NumRegions > 0
            p3 = plot(poly_int, 'FaceColor', neonYellow, 'FaceAlpha', 1, 'EdgeColor', 'none');
        else
            p3 = patch(NaN, NaN, neonYellow, 'FaceAlpha', 1, 'EdgeColor', 'none');
        end
        delete(findobj(gca, 'Type', 'text'));
        %title(sprintf('Variance Partitioning in %s of Subject %d', roiTitle, subject), 'FontSize', 14);
        title(sprintf('Variance Partitioning in %s', roiTitle), 'FontSize', 14);
        lgd = legend([p1, p2, p3], {sprintf('Unique Pose Estimation (R^2 = %.4f)', unique_pose), ...
                                     sprintf('Unique Body Segmentation (R^2 = %.4f)', unique_seg), ...
                                     sprintf('Shared Variance (R^2 = %.4f)', shared_variance)}, ...
                      'Location', 'northeast', 'FontSize', 8);
        lgd.Position(1) = 0.75;
        axis equal;
        axis off;
        hold off;
        vennFileNum = fullfile(roi_subfolder, sprintf('%s_Venn_Numbers.png', roiLabel));
        saveas(figVennNum, vennFileNum);
        close(figVennNum);
        fprintf('DEBUG: Saved Venn diagram (with numbers) for %s to %s\n', roiTitle, vennFileNum);
        
        % Create clean Venn diagram.
        figVennClean = figure('Visible', 'off', 'Position', [100, 100, 1200, 800]);
        hold on;
        plot(poly_left, 'FaceColor', brightGreen, 'FaceAlpha', 1, 'EdgeColor', 'none');
        plot(poly_right, 'FaceColor', strongOrange, 'FaceAlpha', 1, 'EdgeColor', 'none');
        if poly_int.NumRegions > 0
            plot(poly_int, 'FaceColor', neonYellow, 'FaceAlpha', 1, 'EdgeColor', 'none');
        end
        delete(findobj(gca, 'Type', 'text'));
        %title(sprintf('Variance Partitioning in %s of Subject %d', roiTitle, subject), 'FontSize', 14);
        title(sprintf('Variance Partitioning in %s', roiTitle), 'FontSize', 14);

        lgd = legend({'Unique Pose Estimation', 'Unique Body Segmentation', 'Shared Variance'}, ...
                     'Location', 'northeast', 'FontSize', 8);
        lgd.Position(1) = 0.8;
        axis equal;
        axis off;
        hold off;
        vennFileClean = fullfile(roi_subfolder, sprintf('%s_Venn_Clean.png', roiLabel));
        saveas(figVennClean, vennFileClean);
        close(figVennClean);
        fprintf('DEBUG: Saved Venn diagram (clean) for %s to %s\n', roiTitle, vennFileClean);
        
        %% Create Regular Bar Graph (Unique Contributions Only)
        x_labels = {'Pose Estimation', 'Body Segmentation', 'Random Pose Estimation', 'Random Body Segmentation'};
        x = categorical(x_labels, x_labels, 'Ordinal', true);
        figBar = figure('Visible', 'off', 'Position', [100, 100, 1200, 800]);
        hold on;
        bar(x(1), unique_pose, 'FaceColor', brightGreen);
        bar(x(2), unique_seg, 'FaceColor', strongOrange);
        bar(x(3), unique_rpose, 'FaceColor', redColor);
        bar(x(4), unique_rseg, 'FaceColor', greyColor);
        hold off;
        ylabel('Variance Explained (R^2)', 'FontSize', 12);
        %title(sprintf('Unique Variance Contributions in %s of Subject %d', roiTitle, subject), 'FontSize', 14);
        title(sprintf('Unique Variance Contributions in %s', roiTitle), 'FontSize', 14);

        lgd = legend({'Pose Estimation', 'Body Segmentation', 'Random Pose Estimation', 'Random Body Segmentation'}, ...
                     'Location', 'northeast', 'FontSize', 8);
        lgd.Position(1) = 0.8;
        set(gca, 'YTick', uniformYTicks);
        ylim([-0.1, 0.2]);
        ytickformat('%.4f');
        set(gca, 'FontSize', 7); % x-axis and tick labels at font size 7
        barFile = fullfile(roi_subfolder, sprintf('%s_Bar.png', roiLabel));
        saveas(figBar, barFile);
        close(figBar);
        fprintf('DEBUG: Saved Bar graph (unique only) for %s to %s\n', roiTitle, barFile);
        
        %% Create Alternative Bar Graph (Unique Contributions + Full Model)
        x_labels2 = {'Pose Estimation', 'Body Segmentation', 'Random Pose Estimation', 'Random Body Segmentation', 'Full Model'};
        x2 = categorical(x_labels2, x_labels2, 'Ordinal', true);
        figBar2 = figure('Visible', 'off', 'Position', [100, 100, 1200, 800]);
        hold on;
        bar(x2(1), unique_pose, 'FaceColor', brightGreen);
        bar(x2(2), unique_seg, 'FaceColor', strongOrange);
        bar(x2(3), unique_rpose, 'FaceColor', redColor);
        bar(x2(4), unique_rseg, 'FaceColor', greyColor);
        bar(x2(5), full_R2, 'FaceColor', blackColor);
        hold off;
        ylabel('Variance Explained (R^2)', 'FontSize', 12);
        %title(sprintf('Overall Variance Explained in %s of Subject %d', roiTitle, subject), 'FontSize', 14);
        title(sprintf('Overall Variance Explained in %s', roiTitle), 'FontSize', 14);

        lgd = legend({sprintf('Pose Estimation (R^2 = %.4f)', unique_pose), ...
                      sprintf('Body Segmentation (R^2 = %.4f)', unique_seg), ...
                      sprintf('Random Pose Estimation (R^2 = %.4f)', unique_rpose), ...
                      sprintf('Random Body Segmentation (R^2 = %.4f)', unique_rseg), ...
                      sprintf('Full Model (R^2 = %.4f)', full_R2)}, ...
                      'Location', 'northeast', 'FontSize', 8);
        lgd.Position(1) = 0.75;
        set(gca, 'YTick', uniformYTicks);
        ylim([-0.1, 0.2]);
        ytickformat('%.4f');
        set(gca, 'FontSize', 7);
        barFile2 = fullfile(roi_subfolder, sprintf('%s_BarFull.png', roiLabel));
        saveas(figBar2, barFile2);
        close(figBar2);
        fprintf('DEBUG: Saved Alternative Bar graph for %s to %s\n', roiTitle, barFile2);
        
    end
end

toc;  % Display elapsed time
