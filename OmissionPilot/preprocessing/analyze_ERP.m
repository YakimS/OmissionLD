% analyze_ERP.m
% ERP analysis for OmissionPilot - central electrodes
%
% Plots 3x3 grid: rows=ISI type, cols=trial type

clear all;
close all;
clc;

%% ============= PATHS =============
restoredefaultpath
addpath 'D:\matlab_libs\eeglab2025.1.0'
eeglab nogui;

PREPROCESSED_DIR = 'D:\omissionPilot\preprocessed';
OUTPUT_DIR = 'D:\omissionPilot\figures';

%% ============= PARAMETERS =============
% Set default white background and black text for all figures
set(0, 'DefaultFigureColor', 'w');
set(0, 'DefaultAxesColor', 'w');
set(0, 'DefaultAxesXColor', 'k');
set(0, 'DefaultAxesYColor', 'k');
set(0, 'DefaultTextColor', 'k');
set(0, 'DefaultAxesTickLabelInterpreter', 'none');
set(0, 'DefaultLegendInterpreter', 'none');
set(0, 'DefaultColorbarTickLabelInterpreter', 'none');

CENTRAL_CHANNELS = {'E9', 'E186', 'E45', 'E81', 'E132', 'E257'};

% Conditions
CONDITIONS = {'SNGL_STD', 'SNGL_OMI', 'SNGL_DEV', ...
              'DBL100_STD', 'DBL100_OMI', 'DBL100_DEV', ...
              'DBL250_STD', 'DBL250_OMI', 'DBL250_DEV'};

% Colors for trial types
COLORS.standard = [0.2 0.6 0.2];   % green
COLORS.omission = [0.8 0.2 0.2];   % red
COLORS.deviant  = [0.2 0.2 0.8];   % blue

%%
if ~exist(OUTPUT_DIR, 'dir')
    mkdir(OUTPUT_DIR);
end

set_files = dir(fullfile(PREPROCESSED_DIR, '*_preprocessed.set'));

%% Process each file
for f = 1:length(set_files)
    filename = set_files(f).name;
    filepath = PREPROCESSED_DIR;

    EEG = pop_loadset('filename', filename, 'filepath', filepath);
    [~, basename, ~] = fileparts(filename);
    basename = strrep(basename, '_preprocessed', '');

    % Find central channel indices
    chan_idx = [];
    if isfield(EEG.chanlocs, 'labels')
        all_labels = {EEG.chanlocs.labels};
        for c = 1:length(CENTRAL_CHANNELS)
            idx = find(strcmpi(all_labels, CENTRAL_CHANNELS{c}));
            if ~isempty(idx)
                chan_idx = [chan_idx, idx(1)];
            end
        end
    end

    % Extract ERPs per condition (mean and SEM)
    erp_data = struct();
    erp_sem = struct();
    epoch_counts = struct();

    for c = 1:length(CONDITIONS)
        cond = CONDITIONS{c};

        % Find epochs matching this condition
        epoch_idx = [];
        for ep = 1:length(EEG.epoch)
            evt_types = EEG.epoch(ep).eventtype;
            if iscell(evt_types)
                if any(strcmp(evt_types, cond))
                    epoch_idx = [epoch_idx, ep];
                end
            else
                if strcmp(evt_types, cond)
                    epoch_idx = [epoch_idx, ep];
                end
            end
        end

        if ~isempty(epoch_idx)
            % Average across central channels first, then compute mean and SEM across epochs
            data = EEG.data(chan_idx, :, epoch_idx);
            chan_avg = squeeze(mean(data, 1));  % [timepoints x epochs]
            if length(epoch_idx) == 1
                chan_avg = chan_avg';  % ensure correct orientation for single epoch
            end
            erp_data.(cond) = mean(chan_avg, 2)';  % mean across epochs
            erp_sem.(cond) = std(chan_avg, 0, 2)' / sqrt(length(epoch_idx));  % SEM
            epoch_counts.(cond) = length(epoch_idx);
        else
            erp_data.(cond) = nan(1, EEG.pnts);
            erp_sem.(cond) = nan(1, EEG.pnts);
            epoch_counts.(cond) = 0;
        end
    end

    times = EEG.times;

    % Compute global y-limits from all conditions
    all_vals = [];
    for c = 1:length(CONDITIONS)
        cond = CONDITIONS{c};
        if epoch_counts.(cond) > 0
            all_vals = [all_vals, erp_data.(cond) + erp_sem.(cond), erp_data.(cond) - erp_sem.(cond)];
        end
    end
    y_min = min(all_vals);
    y_max = max(all_vals);
    y_margin = (y_max - y_min) * 0;
    Y_LIMITS = [-3,3];%[y_min - y_margin, y_max + y_margin];

    % 3x3 grid - rows=ISI type, cols=trial type
    figure('Position', [100 100 1000 800], 'Color', 'w');

    row_labels = {'Single', 'Double 100ms', 'Double 250ms'};
    col_labels = {'Standard', 'Omission', 'Deviant'};
    row_prefixes = {'SNGL', 'DBL100', 'DBL250'};
    col_suffixes = {'STD', 'OMI', 'DEV'};
    col_colors = {COLORS.standard, COLORS.omission, COLORS.deviant};

    for row = 1:3
        for col = 1:3
            subplot(3, 3, (row-1)*3 + col);
            set(gca, 'Color', 'w');

            cond = [row_prefixes{row} '_' col_suffixes{col}];

            if isfield(erp_data, cond) && epoch_counts.(cond) > 0
                hold on;
                % Plot SEM shading
                upper = erp_data.(cond) + erp_sem.(cond);
                lower = erp_data.(cond) - erp_sem.(cond);
                fill([times, fliplr(times)], [upper, fliplr(lower)], col_colors{col}, ...
                     'FaceAlpha', 0.2, 'EdgeColor', 'none');
                % Plot mean line
                plot(times, erp_data.(cond), 'Color', col_colors{col}, 'LineWidth', 1.5);
                xline(0, 'k--', 'LineWidth', 0.5);
                yline(0, 'k-', 'LineWidth', 0.5);

                % Mark tone onsets for double conditions
                if row == 2  % DBL100
                    xline(100, 'k:', 'LineWidth', 0.5);
                elseif row == 3  % DBL250
                    xline(250, 'k:', 'LineWidth', 0.5);
                end
                hold off;

                title(sprintf('n=%d', epoch_counts.(cond)));
            else
                text(0.5, 0.5, 'No data', 'HorizontalAlignment', 'center', 'Units', 'normalized');
            end

            xlim([times(1) times(end)]);
            ylim(Y_LIMITS);

            if row == 3
                xlabel('Time (ms)');
            end
            if col == 1
                ylabel(sprintf('%s\n(µV)', row_labels{row}));
            end
            if row == 1
                title(sprintf('%s (n=%d)', col_labels{col}, epoch_counts.(cond)), 'FontWeight', 'bold');
            end
        end
    end

    sgtitle(sprintf('%s - ERP by Condition', basename), 'Interpreter', 'none', 'FontSize', 12);

    saveas(gcf, fullfile(OUTPUT_DIR, [basename '_ERP_grid.png']));

    close all;
end
