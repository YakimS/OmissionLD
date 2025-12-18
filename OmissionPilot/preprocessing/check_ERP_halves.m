% check_ERP_halves.m
% Compare first half vs last half of SNGL_STD trials
% To check for fatigue/drift effects across the experiment

clear; close all; clc;

addpath 'D:\matlab_libs\eeglab2025.1.0'
eeglab nogui;

PREPROCESSED_DIR = 'D:\omissionPilot\preprocessed';
OUTPUT_DIR = 'D:\omissionPilot\figures';
CENTRAL_CHANNELS = {'E9', 'E186', 'E45', 'E81', 'E132', 'E257'};

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
COLOR_FIRST = [0.2 0.6 0.2];   % green - first half
COLOR_LAST = [0.8 0.2 0.2];    % red - last half

if ~exist(OUTPUT_DIR, 'dir')
    mkdir(OUTPUT_DIR);
end

set_files = dir(fullfile(PREPROCESSED_DIR, '*_preprocessed.set'));
n_files = length(set_files);

% Storage for grand average
all_erp_first = [];
all_erp_last = [];
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
    chan_idx = find(ismember({EEG.chanlocs.labels}, CENTRAL_CHANNELS));

    % Get epoch indices for SNGL_STD
    std_idx = find(cellfun(@(x) any(strcmp(x, 'SNGL_STD')), {EEG.epoch.eventtype}));
    n_trials = length(std_idx);
    half_point = floor(n_trials / 2);

    % Split into first and last half
    first_half_idx = std_idx(1:half_point);
    last_half_idx = std_idx(half_point+1:end);

    % Compute ERPs (mean across central channels first)
    data_first = squeeze(mean(EEG.data(chan_idx, :, first_half_idx), 1));
    data_last = squeeze(mean(EEG.data(chan_idx, :, last_half_idx), 1));

    % Mean ERP across epochs
    erp_first = mean(data_first, 2)';
    erp_last = mean(data_last, 2)';

    % Within-subject SEM (across epochs)
    sem_first = std(data_first, 0, 2)' / sqrt(length(first_half_idx));
    sem_last = std(data_last, 0, 2)' / sqrt(length(last_half_idx));

    % Store for grand average
    all_erp_first = [all_erp_first; erp_first];
    all_erp_last = [all_erp_last; erp_last];

    % Plot
    subplot(n_files, 1, f);
    set(gca, 'Color', 'w');
    hold on;

    % First half with SEM shading
    fill([EEG.times, fliplr(EEG.times)], [erp_first + sem_first, fliplr(erp_first - sem_first)], ...
        COLOR_FIRST, 'FaceAlpha', 0.2, 'EdgeColor', 'none');
    h1 = plot(EEG.times, erp_first, 'Color', COLOR_FIRST, 'LineWidth', 1.5);

    % Last half with SEM shading
    fill([EEG.times, fliplr(EEG.times)], [erp_last + sem_last, fliplr(erp_last - sem_last)], ...
        COLOR_LAST, 'FaceAlpha', 0.2, 'EdgeColor', 'none');
    h2 = plot(EEG.times, erp_last, 'Color', COLOR_LAST, 'LineWidth', 1.5);

    xline(0, 'k--');
    yline(0, 'k-');
    xlabel('Time (ms)');
    ylabel('Amplitude (µV)');
    title(sprintf('%s - First half (n=%d) vs Last half (n=%d)', basename, length(first_half_idx), length(last_half_idx)), 'Interpreter', 'none');
    legend([h1, h2], {'First half', 'Last half'}, 'Location', 'best');
    hold off;
end

sgtitle('SNGL\_STD: First Half vs Last Half', 'FontSize', 14, 'FontWeight', 'bold');

saveas(gcf, fullfile(OUTPUT_DIR, 'ERP_halves_individual.png'));

%% Grand average across subjects
figure('Position', [100 100 800 500], 'Color', 'w');

% Compute grand averages
grand_first = mean(all_erp_first, 1);
grand_last = mean(all_erp_last, 1);

% SEM for shading (across subjects)
sem_first_grand = std(all_erp_first, 0, 1) / sqrt(n_files);
sem_last_grand = std(all_erp_last, 0, 1) / sqrt(n_files);

set(gca, 'Color', 'w');
hold on;

% First half with SEM shading
fill([times, fliplr(times)], [grand_first + sem_first_grand, fliplr(grand_first - sem_first_grand)], ...
    COLOR_FIRST, 'FaceAlpha', 0.2, 'EdgeColor', 'none');
h1 = plot(times, grand_first, 'Color', COLOR_FIRST, 'LineWidth', 2);

% Last half with SEM shading
fill([times, fliplr(times)], [grand_last + sem_last_grand, fliplr(grand_last - sem_last_grand)], ...
    COLOR_LAST, 'FaceAlpha', 0.2, 'EdgeColor', 'none');
h2 = plot(times, grand_last, 'Color', COLOR_LAST, 'LineWidth', 2);

xline(0, 'k--');
yline(0, 'k-');
xlabel('Time (ms)');
ylabel('Amplitude (µV)');
title(sprintf('Grand Average SNGL\\_STD: First vs Last Half (N=%d subjects)', n_files));
legend([h1, h2], {'First half', 'Last half'}, 'Location', 'best');
hold off;

saveas(gcf, fullfile(OUTPUT_DIR, 'ERP_halves_grand_average.png'));

fprintf('Figures saved to %s\n', OUTPUT_DIR);
