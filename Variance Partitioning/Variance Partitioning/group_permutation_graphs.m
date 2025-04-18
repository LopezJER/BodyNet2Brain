%% -------------------- GROUP-LEVEL SIGNIFICANCE + SYMBOL MAPS --------------------
clear; clc; close all;

%% --- CONFIGURATION ---
input_file = 'D:/ML_project/Variance/var_excel/sapiens/subject_excels/group_level_2/all_subjects_summary.xlsx';
timestamp = datestr(now, 'ddmm_HHMM');
out_dir = fullfile(fileparts(input_file), ['group_sig_outputs_' timestamp]);
if ~exist(out_dir, 'dir'), mkdir(out_dir); end

raw = readtable(input_file);
stats = {'Unique_Pose','Unique_Seg','Full_R2','Shared'};
pcols = {'p_Unique_Pose','p_Unique_Seg','p_Full_R2','p_Shared'};
rois = unique(raw.ROI);
subjects = unique(raw.Subject);
subject_labels = compose('S%d', subjects);

%% --- FDR CORRECTION & SUMMARY TABLE ---
summary_out = {};
for s = 1:length(stats)
    stat_name = stats{s};
    p_name = pcols{s};
    p_vals = raw.(p_name);

    valid_idx = ~isnan(p_vals);
    [~, ~, adj_p] = fdr_bh(p_vals(valid_idx));
    corrected_p = nan(height(raw), 1);
    corrected_p(valid_idx) = adj_p;
    raw.([p_name '_FDR']) = corrected_p;

    sig_raw = p_vals < 0.05;
    sig_fdr = corrected_p < 0.05;

    raw.([p_name '_Sig']) = sig_raw;
    raw.([p_name '_Sig_FDR']) = sig_fdr;

    for r = 1:length(rois)
        roi_name = rois{r};
        roi_idx = strcmp(raw.ROI, roi_name);
        sub_vals = raw{roi_idx, p_name};
        fdr_vals = raw{roi_idx, [p_name '_FDR']};
        sigs = raw{roi_idx, [p_name '_Sig']};
        sigs_fdr = raw{roi_idx, [p_name '_Sig_FDR']};
        n_present = sum(~isnan(sub_vals));

        new_row = {
            stat_name, roi_name, n_present, ...
            sum(sub_vals < 0.05), sum(sub_vals < 0.01), sum(sub_vals < 0.001), sum(sub_vals < 0.0001), ...
            sum(fdr_vals < 0.05), sum(fdr_vals < 0.01), sum(fdr_vals < 0.001), sum(fdr_vals < 0.0001), ...
            sum(sigs), sum(sigs_fdr), sum(sigs & ~sigs_fdr)
        };
        summary_out(end+1, :) = new_row;
    end
end

headers = {'Stat','ROI','N_Subjects','N_p_005','N_p_001','N_p_0001','N_p_00001', ...
           'N_FDR_p_005','N_FDR_p_001','N_FDR_p_0001','N_FDR_p_00001', ...
           'N_Sig_Before','N_Sig_After','N_Lost_Significance'};
summary_out = cell2table(summary_out, 'VariableNames', headers);
writetable(summary_out, fullfile(out_dir, 'fdr_condensed_summary.xlsx'));
writetable(raw, fullfile(out_dir, 'all_subjects_with_fdr.xlsx'));

%% --- SYMBOL MAPS ---
for s = 1:length(stats)
    stat = stats{s};
    pcol = pcols{s};
    fdrcol = [pcol '_FDR'];

    for mode = 1:3
        switch mode
            case 1
                label = 'raw'; col = pcol;
                title_str = sprintf('Symbol Map - %s (Uncorrected)', replace(stat, 'Full_R2', 'Full Model R^2'));
            case 2
                label = 'fdr'; col = fdrcol;
                title_str = sprintf('Symbol Map - %s (FDR Corrected)', replace(stat, 'Full_R2', 'Full Model R^2'));
            case 3
                label = 'comparative';
                title_str = sprintf('Symbol Map - %s (Comparative View)', replace(stat, 'Full_R2', 'Full Model R^2'));
        end

        % Build symbol matrix
        symbol_matrix = strings(length(rois), length(subjects));
        for r = 1:length(rois)
            for sub = 1:length(subjects)
                idx = strcmp(raw.ROI, rois{r}) & raw.Subject == subjects(sub);
                if any(idx)
                    p_raw = raw{idx, pcol};
                    p_fdr = raw{idx, fdrcol};
                    if isnan(p_raw)
                        symbol_matrix(r, sub) = "X";
                    elseif mode == 3 && get_symbol(p_raw) ~= get_symbol(p_fdr)
                        symbol_matrix(r, sub) = get_symbol(p_fdr) + " (" + get_symbol(p_raw) + ")";
                    else
                        symbol_matrix(r, sub) = get_symbol(p_raw);
                    end
                else
                    symbol_matrix(r, sub) = "X";
                end
            end
        end

        % Create UI figure
        fig = uifigure('Color','w','Position',[100 100 1600 950]);

        % Create UI table with taller rows and no outer white box
        t = uitable(fig, 'Data', symbol_matrix, ...
            'ColumnName', subject_labels, 'RowName', rois, ...
            'FontName', 'Consolas', 'FontSize', 18, ...
            'Position', [20 100 1300 750]);
        t.RowStriping = 'off';
        t.RowHeight = repmat({35}, length(rois), 1); % taller rows

        % Title
        uilabel(fig, 'Text', title_str, ...
            'FontSize', 20, 'FontWeight','bold', 'Position', [20 880 1000 30]);

        % Legend box
        panel = uipanel(fig, 'Title', 'Legend', 'FontSize', 12, ...
            'Position', [1350 330 220 250]);

        uilabel(panel, 'Text', sprintf(...
            ['**** : p < 0.0001\n' ...
             '***  : p < 0.001\n' ...
             '**   : p < 0.01\n' ...
             '*    : p < 0.05\n' ...
             '+    : p < 0.1\n' ...
             '-    : ns\n' ...
             'X    : missing']), ...
             'Position', [10 60 200 170], 'FontSize', 11);

        if mode == 2 || mode == 3
            n_tests = sum(~isnan(raw{:,fdrcol}));
            uilabel(panel, 'Text', sprintf("Benjamini-Hochberg FDR\napplied to %d comparisons", n_tests), ...
                'FontSize', 10, 'Position', [10 10 200 40]);
        end

        % Export
        exportapp(fig, fullfile(out_dir, sprintf('symbolmap_%s_%s.png', label, stat)));
        close(fig);
    end
end

fprintf('\n--- FDR + Symbol Maps Complete!\nSaved to: %s\n', out_dir);

%% --- BH FDR FUNCTION ---
function [h, crit_p, adj_p] = fdr_bh(pvals, q, method, report)
    if nargin < 2 || isempty(q), q = 0.05; end
    if nargin < 3 || isempty(method), method = 'pdep'; end
    if nargin < 4, report = 'no'; end
    p = pvals(:);
    [p_sorted, sort_ids] = sort(p);
    [~, unsort_ids] = sort(sort_ids);
    V = length(p);
    I = (1:V)';
    if strcmpi(method,'pdep')
        cVID = 1;
    elseif strcmpi(method,'dep')
        cVID = sum(1./(1:V));
    else
        error('Method must be "pdep" or "dep".');
    end
    thresh = (I/V)*q/cVID;
    wtd_p = V*p_sorted ./ I;
    adj_p = cummin(wtd_p(end:-1:1));
    adj_p = adj_p(end:-1:1);
    h = p_sorted <= thresh;
    crit_p = max(p_sorted(h));
    if isempty(crit_p), crit_p = 0; end
    adj_p = min(1, adj_p);
    adj_p = adj_p(unsort_ids);
end

%% --- SYMBOL MAPPING FUNCTION ---
function sym = get_symbol(p)
    if isnan(p)
        sym = "X";
    elseif p < 1e-4
        sym = "****";
    elseif p < 0.001
        sym = "***";
    elseif p < 0.01
        sym = "**";
    elseif p < 0.05
        sym = "*";
    elseif p < 0.1
        sym = "+";
    else
        sym = "-";
    end
end
