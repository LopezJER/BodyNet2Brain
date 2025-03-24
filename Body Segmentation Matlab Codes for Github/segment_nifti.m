function segment_nifti(input_nii_path)
    % Load NIfTI handling libraries
    if ~exist('niftiinfo', 'file') || ~exist('niftiread', 'file') || ~exist('niftiwrite', 'file')
        error('NIfTI handling functions not available. Please ensure you have the required toolbox installed.');
    end

    % Read the input NIfTI file
    nii_info = niftiinfo(input_nii_path);
    nii_data = niftiread(input_nii_path);

    % Check the first letter of the file name for laterality
    [file_path, file_name, ~] = fileparts(input_nii_path);
    laterality = lower(file_name(1)); % Extract first letter

    if ~ismember(laterality, {'r', 'l'})
        error('The file name does not start with a valid laterality indicator (r/l).');
    end

    % Define label names corresponding to intensity values
    labels = {1, 'EBA'; 2, 'FBA_1'; 3, 'FBA_2'; 4, 'MTL'};

    % Output folder for segmented ROIs
    output_folder = fullfile(file_path, [file_name '_segmented']);
    if ~exist(output_folder, 'dir')
        mkdir(output_folder);
    end

    % Loop through each label and create a new NIfTI file for each ROI
    for i = 1:size(labels, 1)
        label_value = labels{i, 1};
        label_name = labels{i, 2};

        % Create a binary mask for the current label
        roi_mask = (nii_data == label_value);

        if any(roi_mask(:))
            % Generate output file name
            output_nii_name = sprintf('%s%s.nii', laterality, label_name);
            output_nii_path = fullfile(output_folder, output_nii_name);

            % Prepare segmented data in the original datatype
            segmented_data = cast(roi_mask, class(nii_data));

            % Update the NIfTI metadata for the new file
            nii_info_out = nii_info; % Copy metadata
            nii_info_out.Filename = output_nii_path;

            % Write the segmented data
            niftiwrite(segmented_data, output_nii_path, nii_info_out);

            fprintf('Saved segmented ROI: %s\n', output_nii_path);
        else
            fprintf('No voxels found for label %d (%s). Skipping.\n', label_value, label_name);
        end
    end

    fprintf('Segmentation completed. Files saved in: %s\n', output_folder);
end
