% analyze_omission_response.m
% Omission response analysis comparing SNGL_STD vs DBL250_OMI vs DBL100_OMI
% Individual subject plots and grand average across subjects

clear; close all; clc;

addpath 'D:\matlab_libs\eeglab2025.1.0'
eeglab nogui;

PREPROCESSED_DIR = 'D:\omissionPilot\preprocessed';
OUTPUT_DIR = 'D:\omissionPilot\figures';
% central {'E9', 'E186', 'E45', 'E81', 'E132', 'E257'}
% frontal {'E9', 'E186', 'E45', 'E81', 'E132', 'E257'}
ROI_CHANNELS = {'E15', 'E14', 'E22', 'E6', 'E7', 'E16','E23'};

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

if ~exist(OUTPUT_DIR, 'dir')
    mkdir(OUTPUT_DIR);
end

set_files = dir(fullfile(PREPROCESSED_DIR, '*_preprocessed.set'));
n_files = length(set_files);

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

    % Find central channels
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

sgtitle('Omission Response - Individual Subjects', 'FontSize', 14, 'FontWeight', 'bold');

saveas(gcf, fullfile(OUTPUT_DIR, 'omission_response_individual.png'));

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

sgtitle('Omission Response - Grand Average', 'FontSize', 14, 'FontWeight', 'bold');

saveas(gcf, fullfile(OUTPUT_DIR, 'omission_response_grand_average.png'));

fprintf('Figures saved to %s\n', OUTPUT_DIR);
