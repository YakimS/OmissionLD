% analyze_omission_response.m
% Omission response analysis comparing SNGL_STD vs DBL250_OMI vs DBL100_OMI
% Individual subject plots and grand average across subjects

clear; close all; clc;

%% ============= DATA TYPE SELECTION =============
% Choose which data to analyze: 'gel', 'saline', or 'both'
DATA_TYPE = 'both';  % Options: 'gel', 'saline', 'both'

addpath 'D:\matlab_libs\eeglab2025.1.0'
eeglab nogui;

PREPROCESSED_BASE = 'D:\omissionPilot\preprocessed';
OUTPUT_BASE = 'D:\omissionPilot\figures';

% Configure data types to process
if strcmp(DATA_TYPE, 'both')
    data_types = {'gel', 'saline'};
elseif strcmp(DATA_TYPE, 'gel')
    data_types = {'gel'};
elseif strcmp(DATA_TYPE, 'saline')
    data_types = {'saline'};
else
    error('Invalid DATA_TYPE: %s. Must be ''gel'', ''saline'', or ''both''', DATA_TYPE);
end

% Electrode sets per data type
ELECTRODES.gel.central = {'E106', 'E7', 'E80', 'E55', 'E31'};
ELECTRODES.gel.frontal = {'E4', 'E5', 'E10', 'E11', 'E12', 'E16', 'E18', 'E19'};
ELECTRODES.saline.central = {'E9', 'E186', 'E45', 'E81', 'E132', 'E257'};
ELECTRODES.saline.frontal = {'E6', 'E7', 'E14', 'E15', 'E16', 'E22', 'E23'};
ELECTRODE_SETS = {'central', 'frontal'};

% Set default white background and black text for all figures
set(0, 'DefaultFigureColor', 'w');
set(0, 'DefaultAxesColor', 'w');
set(0, 'DefaultAxesXColor', 'k');
set(0, 'DefaultAxesYColor', 'k');
set(0, 'DefaultTextColor', 'k');
set(0, 'DefaultAxesTickLabelInterpreter', 'none');
set(0, 'DefaultLegendInterpreter', 'none');
set(0, 'DefaultColorbarTickLabelInterpreter', 'none');

% Colors
COLOR_STD = [0.2 0.6 0.2];       % green - SNGL_STD
COLOR_OMI250 = [0.8 0.2 0.2];   % red - DBL250_OMI
COLOR_OMI100 = [0.2 0.2 0.8];   % blue - DBL100_OMI

%% Process each data type
for dt = 1:length(data_types)
    current_type = data_types{dt};
    PREPROCESSED_DIR = fullfile(PREPROCESSED_BASE, current_type);
    OUTPUT_DIR = fullfile(OUTPUT_BASE, current_type);

    if ~exist(OUTPUT_DIR, 'dir')
        mkdir(OUTPUT_DIR);
    end

    set_files = dir(fullfile(PREPROCESSED_DIR, '*_preprocessed.set'));
    n_files = length(set_files);

    if n_files == 0
        continue;
    end

%% Process each electrode set
for es = 1:length(ELECTRODE_SETS)
    current_elec_set = ELECTRODE_SETS{es};
    ROI_CHANNELS = ELECTRODES.(current_type).(current_elec_set);

    % Storage for grand average
    all_erp_std = [];
    all_erp_omi250 = [];
    all_erp_omi100 = [];
    times = [];

%% Individual subject plots
figure('Position', [100 100 1200 300*n_files], 'Color', 'w');

for f = 1:n_files
    EEG = pop_loadset('filename', set_files(f).name, 'filepath', PREPROCESSED_DIR);
    [~, basename] = fileparts(set_files(f).name);
    basename = strrep(basename, '_preprocessed', '');

    % Store times (same for all subjects)
    if isempty(times)
        times = EEG.times;
    end

    % Find channels
    chan_idx = find(ismember({EEG.chanlocs.labels}, ROI_CHANNELS));

    % Get epoch indices for each condition
    std_idx = find(cellfun(@(x) any(strcmp(x, 'SNGL_STD')), {EEG.epoch.eventtype}));
    omi250_idx = find(cellfun(@(x) any(strcmp(x, 'DBL250_OMI')), {EEG.epoch.eventtype}));
    omi100_idx = find(cellfun(@(x) any(strcmp(x, 'DBL100_OMI')), {EEG.epoch.eventtype}));

    % Compute ERPs (mean across central channels, then across epochs)
    if ~isempty(std_idx)
        erp_std = mean(mean(EEG.data(chan_idx, :, std_idx), 1), 3);
    else
        erp_std = nan(1, EEG.pnts);
    end
    if ~isempty(omi250_idx)
        erp_omi250 = mean(mean(EEG.data(chan_idx, :, omi250_idx), 1), 3);
    else
        erp_omi250 = nan(1, EEG.pnts);
    end
    if ~isempty(omi100_idx)
        erp_omi100 = mean(mean(EEG.data(chan_idx, :, omi100_idx), 1), 3);
    else
        erp_omi100 = nan(1, EEG.pnts);
    end

    % Store for grand average
    all_erp_std = [all_erp_std; erp_std];
    all_erp_omi250 = [all_erp_omi250; erp_omi250];
    all_erp_omi100 = [all_erp_omi100; erp_omi100];

    % Plot
    subplot(n_files, 1, f);
    set(gca, 'Color', 'w');
    hold on;
    h1 = plot(EEG.times, erp_std, 'Color', COLOR_STD, 'LineWidth', 1.5);
    h2 = plot(EEG.times, erp_omi250, 'Color', COLOR_OMI250, 'LineWidth', 1.5);
    h3 = plot(EEG.times, erp_omi100, 'Color', COLOR_OMI100, 'LineWidth', 1.5);
    xline(0, 'k--');
    yline(0, 'k-');
    % Mark expected tone2 onset times
    xline(100, 'Color', COLOR_OMI100, 'LineStyle', ':', 'LineWidth', 1);
    xline(250, 'Color', COLOR_OMI250, 'LineStyle', ':', 'LineWidth', 1);
    xlabel('Time (ms)');
    ylabel('Amplitude (µV)');
    title(sprintf('%s - STD(n=%d) OMI250(n=%d) OMI100(n=%d)', ...
        basename, length(std_idx), length(omi250_idx), length(omi100_idx)), 'Interpreter', 'none');
    legend([h1, h2, h3], {'SNGL\_STD', 'DBL250\_OMI', 'DBL100\_OMI'}, 'Location', 'best');
    hold off;
end

sgtitle(sprintf('Omission Response - Individual Subjects (%s)', current_elec_set), 'FontSize', 14, 'FontWeight', 'bold');

saveas(gcf, fullfile(OUTPUT_DIR, ['omission_response_individual_' current_elec_set '.png']));

%% Grand average across subjects
figure('Position', [100 100 800 600], 'Color', 'w');

% Compute grand averages
grand_std = mean(all_erp_std, 1, 'omitnan');
grand_omi250 = mean(all_erp_omi250, 1, 'omitnan');
grand_omi100 = mean(all_erp_omi100, 1, 'omitnan');

% SEM for shading
n_valid_std = sum(~isnan(all_erp_std(:,1)));
n_valid_omi250 = sum(~isnan(all_erp_omi250(:,1)));
n_valid_omi100 = sum(~isnan(all_erp_omi100(:,1)));

sem_std = std(all_erp_std, 0, 1, 'omitnan') / sqrt(n_valid_std);
sem_omi250 = std(all_erp_omi250, 0, 1, 'omitnan') / sqrt(n_valid_omi250);
sem_omi100 = std(all_erp_omi100, 0, 1, 'omitnan') / sqrt(n_valid_omi100);

set(gca, 'Color', 'w');
hold on;

% SNGL_STD with SEM shading
fill([times, fliplr(times)], [grand_std + sem_std, fliplr(grand_std - sem_std)], ...
    COLOR_STD, 'FaceAlpha', 0.2, 'EdgeColor', 'none');
h1 = plot(times, grand_std, 'Color', COLOR_STD, 'LineWidth', 2);

% DBL250_OMI with SEM shading
fill([times, fliplr(times)], [grand_omi250 + sem_omi250, fliplr(grand_omi250 - sem_omi250)], ...
    COLOR_OMI250, 'FaceAlpha', 0.2, 'EdgeColor', 'none');
h2 = plot(times, grand_omi250, 'Color', COLOR_OMI250, 'LineWidth', 2);

% DBL100_OMI with SEM shading
fill([times, fliplr(times)], [grand_omi100 + sem_omi100, fliplr(grand_omi100 - sem_omi100)], ...
    COLOR_OMI100, 'FaceAlpha', 0.2, 'EdgeColor', 'none');
h3 = plot(times, grand_omi100, 'Color', COLOR_OMI100, 'LineWidth', 2);

xline(0, 'k--');
yline(0, 'k-');
% Mark expected tone2 onset times
xline(100, 'Color', COLOR_OMI100, 'LineStyle', ':', 'LineWidth', 1.5);
xline(250, 'Color', COLOR_OMI250, 'LineStyle', ':', 'LineWidth', 1.5);

xlabel('Time (ms)');
ylabel('Amplitude (µV)');
title(sprintf('Grand Average Omission Response (N=%d subjects)', n_files));
legend([h1, h2, h3], {'SNGL\_STD', 'DBL250\_OMI', 'DBL100\_OMI'}, 'Location', 'best');
hold off;

sgtitle(sprintf('Omission Response - Grand Average (%s)', current_elec_set), 'FontSize', 14, 'FontWeight', 'bold');

saveas(gcf, fullfile(OUTPUT_DIR, ['omission_response_grand_average_' current_elec_set '.png']));

end  % end electrode_sets loop

end  % end data_types loop
