% quick_ERP_check.m
% Simple ERP verification for SNGL_STD vs SNGL_DEV
% Includes grand average across subjects with difference wave

clear; close all; clc;

addpath 'D:\matlab_libs\eeglab2025.1.0'
eeglab nogui;

PREPROCESSED_DIR = 'D:\omissionPilot\preprocessed';
ROI_CHANNELS = {'E15', 'E14', 'E22', 'E6', 'E7', 'E16','E23'};


set_files = dir(fullfile(PREPROCESSED_DIR, '*_preprocessed.set'));
n_files = length(set_files);

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
COLOR_STD = [0.2 0.6 0.2];   % green
COLOR_DEV = [0.2 0.2 0.8];   % blue
COLOR_DIFF = [0.5 0.0 0.5];  % purple for difference wave

% Storage for grand average
all_erp_std = [];
all_erp_dev = [];
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

    % Find central channels
    chan_idx = find(ismember({EEG.chanlocs.labels}, ROI_CHANNELS));

    % Get epoch indices for each condition
    std_idx = find(cellfun(@(x) any(strcmp(x, 'SNGL_STD')), {EEG.epoch.eventtype}));
    dev_idx = find(cellfun(@(x) any(strcmp(x, 'SNGL_DEV')), {EEG.epoch.eventtype}));

    % Compute ERPs (mean across central channels first)
    % data_std: [timepoints x epochs]
    data_std = squeeze(mean(EEG.data(chan_idx, :, std_idx), 1));
    data_dev = squeeze(mean(EEG.data(chan_idx, :, dev_idx), 1));

    % Mean ERP across epochs
    erp_std = mean(data_std, 2)';
    erp_dev = mean(data_dev, 2)';

    % Within-subject SEM (across epochs)
    sem_std_subj = std(data_std, 0, 2)' / sqrt(length(std_idx));
    sem_dev_subj = std(data_dev, 0, 2)' / sqrt(length(dev_idx));

    % Store for grand average
    all_erp_std = [all_erp_std; erp_std];
    all_erp_dev = [all_erp_dev; erp_dev];

    % Plot
    subplot(n_files, 1, f);
    set(gca, 'Color', 'w');
    hold on;

    % Standard with SEM shading
    fill([EEG.times, fliplr(EEG.times)], [erp_std + sem_std_subj, fliplr(erp_std - sem_std_subj)], ...
        COLOR_STD, 'FaceAlpha', 0.2, 'EdgeColor', 'none');
    h1 = plot(EEG.times, erp_std, 'Color', COLOR_STD, 'LineWidth', 1.5);

    % Deviant with SEM shading
    fill([EEG.times, fliplr(EEG.times)], [erp_dev + sem_dev_subj, fliplr(erp_dev - sem_dev_subj)], ...
        COLOR_DEV, 'FaceAlpha', 0.2, 'EdgeColor', 'none');
    h2 = plot(EEG.times, erp_dev, 'Color', COLOR_DEV, 'LineWidth', 1.5);

    xline(0, 'k--');
    yline(0, 'k-');
    xlabel('Time (ms)');
    ylabel('Amplitude (µV)');
    title(sprintf('%s - STD(n=%d) DEV(n=%d)', basename, length(std_idx), length(dev_idx)), 'Interpreter', 'none');
    legend([h1, h2], {'Standard', 'Deviant'}, 'Location', 'best');
    hold off;
end

sgtitle('Individual Subject ERPs', 'FontSize', 14, 'FontWeight', 'bold');

%% Grand average across subjects
figure('Position', [100 100 1200 500], 'Color', 'w');

% Compute grand averages
grand_std = mean(all_erp_std, 1);
grand_dev = mean(all_erp_dev, 1);
grand_diff = grand_std - grand_dev;

% SEM for shading
sem_std = std(all_erp_std, 0, 1) / sqrt(n_files);
sem_dev = std(all_erp_dev, 0, 1) / sqrt(n_files);
sem_diff = std(all_erp_std - all_erp_dev, 0, 1) / sqrt(n_files);

% Plot 1: Grand average ERPs
subplot(1, 2, 1);
set(gca, 'Color', 'w');
hold on;

% Standard with SEM shading
fill([times, fliplr(times)], [grand_std + sem_std, fliplr(grand_std - sem_std)], ...
    COLOR_STD, 'FaceAlpha', 0.2, 'EdgeColor', 'none');
h1 = plot(times, grand_std, 'Color', COLOR_STD, 'LineWidth', 2);

% Deviant with SEM shading
fill([times, fliplr(times)], [grand_dev + sem_dev, fliplr(grand_dev - sem_dev)], ...
    COLOR_DEV, 'FaceAlpha', 0.2, 'EdgeColor', 'none');
h2 = plot(times, grand_dev, 'Color', COLOR_DEV, 'LineWidth', 2);

xline(0, 'k--');
yline(0, 'k-');
xlabel('Time (ms)');
ylabel('Amplitude (µV)');
title(sprintf('Grand Average (N=%d subjects)', n_files));
legend([h1, h2], {'Standard', 'Deviant'}, 'Location', 'best');
hold off;

% Plot 2: Grand average difference wave
subplot(1, 2, 2);
set(gca, 'Color', 'w');
hold on;
fill([times, fliplr(times)], [grand_diff + sem_diff, fliplr(grand_diff - sem_diff)], ...
    COLOR_DIFF, 'FaceAlpha', 0.2, 'EdgeColor', 'none');
plot(times, grand_diff, 'Color', COLOR_DIFF, 'LineWidth', 2);
xline(0, 'k--');
yline(0, 'k-');
xlabel('Time (ms)');
ylabel('Amplitude (µV)');
title(sprintf('Grand Average Difference (STD - DEV), N=%d', n_files));
hold off;

sgtitle('Grand Average ERPs Across Subjects', 'FontSize', 14, 'FontWeight', 'bold');
