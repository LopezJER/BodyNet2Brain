% Create bar and venn graphs from excel - separated lateralities.

%% Variance Partitioning Graphs from Excel with Fixed Y-Ticks
% This script reads an Excel file (with columns: ROI, Unique_Pose, Unique_Seg, Shared_Variance, Full_R2)
% and iterates over a preset pool of ROIs and laterality options (e.g., rEBA, lEBA, etc.).
% For each valid entry, it creates and saves:
%   1. Two Venn diagrams:
%      - Version 1: With a legend that includes R^2 numbers.
%      - Version 2: With a "clean" legend (only descriptive text).
%   2. A grouped bar graph (unique contributions only) that uses a legend with only the names.
%   3. An alternative grouped bar graph (unique contributions + full model) that retains numbers in its legend.
%
% The Venn diagrams clip negative values to zero.
% For all bar graphs, a fixed y-axis range from -0.1 to 0.2 is used with ticks every 0.02.
%
% Figures are saved in a subject-specific subfolder under the output folder.
% Debug messages are printed throughout, and execution time is measured with tic/toc.

tic;  % Start timing
clear; clc;

%% Configuration
subject = 8;
%subjectID = num2str(subject);
input_excel = sprintf('D:\\ML_project\\Variance\\var_excel\\allmodels_sanitized\\subject_%d_variance_partitioning.xlsx', subject);

roiPool = {'EBA', 'FBA_1', 'FBA_2', 'FFA_1', 'FFA_2', 'hV4', 'MTL', 'OFA', ...
           'V1d', 'V1v', 'V2d', 'V2v', 'V3d', 'V3v'};
laterality_options = {'r', 'l', 'm'};  % 'r' for right, 'l' for left, 'm' for merged

% Define colors for the figures
brightGreen = [0, 0.8, 0];       % Unique Pose Estimation (green)
strongOrange = [1, 0.6, 0];      % Unique Body Segmentation (orange)
neonYellow = [1, 1, 0];          % Shared Variance (for Venn diagram)
blackColor   = [0, 0, 0];        % Full Model (black)

% Define fixed y-axis ticks and limits for all bar graphs
uniformYTicks = -0.1:0.02:0.2;   % ticks from -0.1 to 0.2 every 0.02

% Define output folder and create a subject-specific subfolder
base_output_folder = 'D:\ML_project\Variance\var_graphs';
subject_folder = fullfile(base_output_folder, sprintf('subject_%d', subject));
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
        laterality = laterality_options{lat};  % 'r' or 'l'
        % Determine laterality word for titles (e.g., 'Right' or 'Left')
        if strcmpi(laterality, 'r')
            latWord = 'Right';
        elseif strcmpi(laterality, 'l')
            latWord = 'Left';
        else
            latWord = 'Merged';
        end
        
        % Create ROI label for file naming (e.g., 'rEBA')
        roiLabel = sprintf('%s%s', lower(laterality), roi);
        % Create a publication-ready ROI title (e.g., 'Right EBA')
        roiTitle = sprintf('%s %s', latWord, roi);
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
        
        % Extract raw values from the table (for bar graphs)
        unique_pose = dataTable.Unique_Pose(idx);
        unique_seg = dataTable.Unique_Seg(idx);
        shared_variance = dataTable.Shared_Variance(idx);
        full_R2 = dataTable.Full_R2(idx);
        
        % In case multiple rows exist, take the first one.
        unique_pose = unique_pose(1);
        unique_seg = unique_seg(1);
        shared_variance = shared_variance(1);
        full_R2 = full_R2(1);
        
        fprintf('DEBUG: Raw values for %s: Unique Pose = %.4f, Unique Seg = %.4f, Shared = %.4f, Full R^2 = %.4f\n', ...
            roiTitle, unique_pose, unique_seg, shared_variance, full_R2);
        
        % For the Venn diagram, clip negative values to zero.
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
            0.5 * sqrt( max(0, (-d+r1+r2)*(d+r1-r2)*(d-r1+r2)*(d+r1+r2) ) );
        
        % Determine center distance d
        lb = max(abs(r_pose - r_seg), eps);
        ub = r_pose + r_seg - eps;
        max_int = intersection_area(lb, r_pose, r_seg);
        if A_shared > max_int
            fprintf('DEBUG: Adjusting A_shared for %s from %.4f to max_int %.4f\n', roiTitle, A_shared, max_int);
            A_shared = max_int;
        end
        
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
        fprintf('DEBUG: Computed center distance d for %s: %.4f\n', roiTitle, d);
        
        % Position circles so that the Pose Estimation (green) circle is on the left.
        center_left = [-d/2, 0];
        center_right = [d/2, 0];
        theta = linspace(0, 2*pi, 200);
        theta(end) = [];
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
        % Plot Pose (green) on the left and Segmentation (orange) on the right.
        plot(poly_left, 'FaceColor', brightGreen, 'FaceAlpha', 1, 'EdgeColor', 'none');
        plot(poly_right, 'FaceColor', strongOrange, 'FaceAlpha', 1, 'EdgeColor', 'none');
        if poly_int.NumRegions > 0
            plot(poly_int, 'FaceColor', neonYellow, 'FaceAlpha', 1, 'EdgeColor', 'none');
        end
        delete(findobj(gca, 'Type', 'text'));
        title(sprintf('Variance Partitioning in %s of Subject %d', roiTitle, subject), 'FontSize', 14);
        lgd = legend({sprintf('Unique Pose Estimation (R^2 = %.4f)', venn_unique_pose), ...
                      sprintf('Unique Body Segmentation (R^2 = %.4f)', venn_unique_seg), ...
                      sprintf('Shared Variance (R^2 = %.4f)', venn_shared_variance)}, ...
                      'Location', 'northeast', 'FontSize', 8);
        % Move the legend further to the right.
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
        title(sprintf('Variance Partitioning in %s of Subject %d', roiTitle, subject), 'FontSize', 14);
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
        % Define x-axis categories: Pose Estimation (green, left) and Body Segmentation (orange, right)
        x_labels = {'Pose Estimation', 'Body Segmentation'};
        x = categorical(x_labels, x_labels, 'Ordinal', true);
        figBar = figure('Visible', 'off', 'Position', [100, 100, 1200, 800]);
        hold on;
        bar(x(1), unique_pose, 'FaceColor', brightGreen);
        bar(x(2), unique_seg, 'FaceColor', strongOrange);
        hold off;
        ylabel('Variance Explained (R^2)', 'FontSize', 12);
        title(sprintf('Unique Variance Contributions in %s of Subject %d', roiTitle, subject), 'FontSize', 14);
        lgd = legend({'Pose Estimation', 'Body Segmentation'}, 'Location', 'northeast', 'FontSize', 8);
        lgd.Position(1) = 0.8;
        set(gca, 'YTick', uniformYTicks);
        ylim([-0.1, 0.2]);
        ytickformat('%.4f');
        barFile = fullfile(roi_subfolder, sprintf('%s_Bar.png', roiLabel));
        saveas(figBar, barFile);
        close(figBar);
        fprintf('DEBUG: Saved Bar graph (unique only) for %s to %s\n', roiTitle, barFile);
        
        %% Create Alternative Bar Graph (Unique Contributions + Full Model)
        % Order: "Pose Estimation" (green, left), "Body Segmentation" (orange, middle), "Full Model" (black, right)
        x_labels2 = {'Pose Estimation', 'Body Segmentation', 'Full Model'};
        x2 = categorical(x_labels2, x_labels2, 'Ordinal', true);
        figBar2 = figure('Visible', 'off', 'Position', [100, 100, 1200, 800]);
        hold on;
        bar(x2(1), unique_pose, 'FaceColor', brightGreen);
        bar(x2(2), unique_seg, 'FaceColor', strongOrange);
        bar(x2(3), full_R2, 'FaceColor', blackColor);
        hold off;
        ylabel('Variance Explained (R^2)', 'FontSize', 12);
        title(sprintf('Overall Variance Explained in %s of Subject %d', roiTitle, subject), 'FontSize', 14);
        lgd = legend({sprintf('Pose Estimation (R^2 = %.4f)', unique_pose), ...
                      sprintf('Body Segmentation (R^2 = %.4f)', unique_seg), ...
                      sprintf('Full Model (R^2 = %.4f)', full_R2)}, ...
                      'Location', 'northeast', 'FontSize', 8);
        lgd.Position(1) = 0.8;
        set(gca, 'YTick', uniformYTicks);
        ylim([-0.1, 0.2]);
        ytickformat('%.4f');
        barFile2 = fullfile(roi_subfolder, sprintf('%s_BarFull.png', roiLabel));
        saveas(figBar2, barFile2);
        close(figBar2);
        fprintf('DEBUG: Saved Alternative Bar graph for %s to %s\n', roiTitle, barFile2);
        
    end
end

toc;  % Display elapsed time
