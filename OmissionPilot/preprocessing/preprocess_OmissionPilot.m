% preprocess_OmissionPilot.m
% Preprocessing pipeline for OmissionPilot EEG data using EEGLAB
%
% Processes only events where CNAM contains a specified string (e.g., ITI1200ms)
% Two-stage preprocessing with intermediate saving:
%   Stage 1: Channel selection, event relabeling, bandpass + notch filter
%   Stage 2: Bad channel detection/interpolation, re-reference, epoch, baseline, artifact rejection
%
% Conditions (epoched from trial onset, 1.2s epochs):
%   Single: SNGL_STD, SNGL_OMI, SNGL_DEV
%   Double ISI=100ms: DBL100_STD, DBL100_OMI, DBL100_DEV
%   Double ISI=250ms: DBL250_STD, DBL250_OMI, DBL250_DEV

clear all;
close all;
clc;

%% ============= PATHS =============
restoredefaultpath
addpath 'D:\matlab_libs\eeglab2025.1.0'
eeglab %nogui;

if ~exist('pop_mffimport', 'file')
    plugin_askinstall('mffmatlabio', 'pop_mffimport', true);
end
if ~exist('zapline_plus', 'file')
    plugin_askinstall('zapline-plus', 'zapline_plus', true);
end
if ~exist('pop_iclabel', 'file')
    plugin_askinstall('ICLabel', 'pop_iclabel', true);
end

RAW_DATA_DIR = 'D:\omissionPilot\raw_data';
OUTPUT_DIR = 'D:\omissionPilot\preprocessed';

%% ============= PARAMETERS =============
% Event filter - only process events where CNAM contains this string
CNAM_FILTER = 'ITI1200ms';

% Processing control flags
FORCE_REPROCESS = true;  % Set to true to reprocess even if output exists

% Filtering parameters
HIGHPASS_FREQ = 0.1;
LOWPASS_FREQ = 50;
NOTCH_FREQ = 50;  % Hz, line noise

% Epoching parameters
EPOCH_START = -0.1;
EPOCH_END = 1.2;
BASELINE_START = -100;  % ms
BASELINE_END = 0;       % ms
ARTIFACT_THRESHOLD = 150;  % µV

% Channel rejection thresholds (leaner alternative to clean_artifacts + ICA)
VARIANCE_THRESHOLD = 4;     % z-score threshold for variance-based rejection
CORRELATION_THRESHOLD = 0.4; % minimum average correlation with neighbors
KURTOSIS_THRESHOLD = 5;     % z-score threshold for kurtosis-based rejection

% ICA parameters - conservative thresholds (>90% confidence for artifact removal)
USE_ICA = true;             % Set to false to skip ICA
ICA_THRESHOLD = 0.90;       % Only remove components with >90% probability of being artifact

% Condition event codes (created during event relabeling)
CONDITION_EVENTS = {'SNGL_STD', 'SNGL_OMI', 'SNGL_DEV', ...
                    'DBL100_STD', 'DBL100_OMI', 'DBL100_DEV', ...
                    'DBL250_STD', 'DBL250_OMI', 'DBL250_DEV'};


%%
if ~exist(OUTPUT_DIR, 'dir')
    mkdir(OUTPUT_DIR);
end

mff_files = dir(fullfile(RAW_DATA_DIR, '*.mff'));

%% Process each file
for f = 1:length(mff_files)
    filename = mff_files(f).name;
    filepath = fullfile(RAW_DATA_DIR, filename);
    [~, basename, ~] = fileparts(filename);

    % Define output filenames
    stage1_filename = [basename '_stage1.set'];
    stage1_path = fullfile(OUTPUT_DIR, stage1_filename);
    output_filename = [basename '_preprocessed.set'];
    output_path = fullfile(OUTPUT_DIR, output_filename);
    summary_filename = [basename '_preprocessing_summary.txt'];
    summary_path = fullfile(OUTPUT_DIR, summary_filename);

    % Initialize summary structure
    summary = struct();
    summary.filename = filename;
    summary.cnam_filter = CNAM_FILTER;
    summary.processing_date = datestr(now, 'yyyy-mm-dd HH:MM:SS');

    % Check if final output already exists
    if exist(output_path, 'file') && ~FORCE_REPROCESS
        fprintf('Skipping %s - output already exists\n', basename);
        continue;
    end

    fprintf('\n========== Processing %s ==========\n', basename);

    % ============= STAGE 1: Channel selection, event filtering/relabeling, filtering =============
    if exist(stage1_path, 'file') && ~FORCE_REPROCESS
        fprintf('Loading Stage 1 output from %s\n', stage1_filename);
        EEG = pop_loadset('filename', stage1_filename, 'filepath', OUTPUT_DIR);
        % Load stage1 summary if exists
        stage1_summary_file = fullfile(OUTPUT_DIR, [basename '_stage1_summary.mat']);
        if exist(stage1_summary_file, 'file')
            load(stage1_summary_file, 'stage1_summary');
            summary.original_events = stage1_summary.original_events;
            summary.events_after_cnam_filter = stage1_summary.events_after_cnam_filter;
            summary.original_channels = stage1_summary.original_channels;
            summary.channels_kept = stage1_summary.channels_kept;
        end
    else
        fprintf('--- Stage 1: Loading, channel selection, event filtering, filtering ---\n');

        EEG = pop_mffimport(filepath, {'code'});
        EEG = eeg_checkset(EEG);
        EEG.setname = basename;

        summary.original_events = length(EEG.event);
        summary.original_channels = EEG.nbchan;

        % Filter events by CNAM field
        EEG = filter_events_by_cnam(EEG, CNAM_FILTER);
        EEG = eeg_checkset(EEG);
        summary.events_after_cnam_filter = length(EEG.event);
        fprintf('Events after CNAM filter (%s): %d (from %d)\n', CNAM_FILTER, summary.events_after_cnam_filter, summary.original_events);

        % Select only scalp channels based on Z coordinate (Z → inferior (−) ↔ superior / vertex (+))
        labels = {EEG.chanlocs.labels};
        Z = [EEG.chanlocs.Z];
        SCALP_Z_THRESHOLD = -0.1;
        scalp_idx = Z > SCALP_Z_THRESHOLD; 
        KEEP_CHANNELS = labels(scalp_idx); % Exclude very low  (neck/face channels)
        EEG = pop_select(EEG, 'channel', find(scalp_idx));
        EEG = eeg_checkset(EEG);
        summary.channels_kept = EEG.nbchan;
        summary.scalp_channels = KEEP_CHANNELS;
        fprintf('Selected %d scalp channels (Z > %.2f)\n', EEG.nbchan, SCALP_Z_THRESHOLD);

        % Relabel events based on condition type and trial type
        [EEG, relabel_counts] = relabel_events(EEG, basename);
        EEG = eeg_checkset(EEG);
        summary.relabeled_events = relabel_counts;

        % Bandpass filtering
        EEG = pop_eegfiltnew(EEG, 'locutoff', HIGHPASS_FREQ);
        EEG = pop_eegfiltnew(EEG, 'hicutoff', LOWPASS_FREQ);
        EEG = eeg_checkset(EEG);

        % Notch filter for line noise removal (50Hz)
        EEG = pop_eegfiltnew(EEG, 'locutoff', NOTCH_FREQ-5, 'hicutoff', NOTCH_FREQ+5, 'revfilt', 1);
        EEG = eeg_checkset(EEG);

        % Save Stage 1 output
        fprintf('Saving Stage 1 output to %s\n', stage1_filename);
        EEG = pop_saveset(EEG, 'filename', stage1_filename, 'filepath', OUTPUT_DIR);

        % Save stage1 summary
        stage1_summary = struct();
        stage1_summary.original_events = summary.original_events;
        stage1_summary.events_after_cnam_filter = summary.events_after_cnam_filter;
        stage1_summary.original_channels = summary.original_channels;
        stage1_summary.channels_kept = summary.channels_kept;
        save(fullfile(OUTPUT_DIR, [basename '_stage1_summary.mat']), 'stage1_summary');
    end

    %% ============= STAGE 2: Bad channel detection, interpolation, re-reference, epoch, artifact rejection =============
    fprintf('--- Stage 2: Bad channel detection, interpolation, epoching, artifact rejection ---\n');

    % Lean channel rejection (variance, correlation, kurtosis-based)
    [bad_chans, bad_chan_labels, bad_chan_details] = detect_bad_channels(EEG, VARIANCE_THRESHOLD, CORRELATION_THRESHOLD, KURTOSIS_THRESHOLD);
    nRemovedChans = length(bad_chans);
    summary.bad_channels_detected = nRemovedChans;
    summary.bad_channel_labels = bad_chan_labels;
    summary.bad_channel_details = bad_chan_details;
    fprintf('Detected %d bad channels\n', nRemovedChans);

    % Interpolate bad channels
    if nRemovedChans > 0
        EEG = pop_interp(EEG, bad_chans, 'spherical');
        EEG = eeg_checkset(EEG);
        fprintf('Interpolated %d channels\n', nRemovedChans);
    end
    summary.channels_interpolated = nRemovedChans;

    % Re-reference to average
    EEG = pop_reref(EEG, []);
    EEG = eeg_checkset(EEG);

    % Epoch data
    n_events_before_epoch = length(EEG.event);
    EEG = pop_epoch(EEG, CONDITION_EVENTS, [EPOCH_START EPOCH_END]);
    EEG = eeg_checkset(EEG);
    summary.epochs_created = EEG.trials;
    fprintf('Created %d epochs\n', EEG.trials);

    % Baseline correction
    EEG = pop_rmbase(EEG, [BASELINE_START BASELINE_END]);
    EEG = eeg_checkset(EEG);

    % ICA for artifact removal (eye movements and muscle artifacts)
    if USE_ICA
        fprintf('--- Running ICA ---\n');

        % Run ICA (using runica algorithm)
        EEG = pop_runica(EEG, 'icatype', 'runica', 'extended', 1);
        EEG = eeg_checkset(EEG);

        % Classify components using ICLabel
        EEG = pop_iclabel(EEG, 'default');

        % Get ICLabel classifications
        % Categories: Brain, Muscle, Eye, Heart, Line Noise, Channel Noise, Other
        % Index:      1      2       3    4      5           6              7
        ic_classes = EEG.etc.ic_classification.ICLabel.classifications;

        % Find components with >90% probability of being Muscle (2) or Eye (3)
        muscle_probs = ic_classes(:, 2);
        eye_probs = ic_classes(:, 3);

        ics_to_remove = find(muscle_probs > ICA_THRESHOLD | eye_probs > ICA_THRESHOLD);

        summary.ica_components_total = size(EEG.icaweights, 1);
        summary.ica_muscle_removed = sum(muscle_probs > ICA_THRESHOLD);
        summary.ica_eye_removed = sum(eye_probs > ICA_THRESHOLD);
        summary.ica_components_removed = length(ics_to_remove);

        fprintf('ICA: %d components total, removing %d (muscle: %d, eye: %d) with >%.0f%% confidence\n', ...
            summary.ica_components_total, length(ics_to_remove), ...
            summary.ica_muscle_removed, summary.ica_eye_removed, ICA_THRESHOLD*100);

        % Remove artifact components
        if ~isempty(ics_to_remove)
            EEG = pop_subcomp(EEG, ics_to_remove, 0);
            EEG = eeg_checkset(EEG);
        end
    else
        summary.ica_components_total = 0;
        summary.ica_components_removed = 0;
        summary.ica_muscle_removed = 0;
        summary.ica_eye_removed = 0;
    end

    n_epochs_before_rejection = EEG.trials;
    original_trial_indices = 1:n_epochs_before_rejection;
    EEG = pop_eegthresh(EEG, 1, 1:EEG.nbchan, -ARTIFACT_THRESHOLD, ARTIFACT_THRESHOLD, ...
                       EPOCH_START, EPOCH_END, 0, 0);  % last 0 = don't remove yet

    % Get rejection mask
    if isfield(EEG, 'reject') && isfield(EEG.reject, 'rejthresh')
        rejected_mask = EEG.reject.rejthresh;
    else
        rejected_mask = zeros(1, n_epochs_before_rejection);
    end
    rejected_indices = find(rejected_mask);
    retained_indices = find(~rejected_mask);

    % Save rejection plot
    plot_rejection_filename = [basename '_rejection_plot.png'];
    plot_rejection_path = fullfile(OUTPUT_DIR, plot_rejection_filename);
    save_rejection_plot(rejected_mask, basename, plot_rejection_path);
    fprintf('Saved rejection plot to %s\n', plot_rejection_filename);

    % Now actually remove rejected epochs
    EEG = pop_rejepoch(EEG, rejected_mask, 0);
    EEG = eeg_checkset(EEG);

    summary.epochs_rejected = length(rejected_indices);
    summary.epochs_retained = EEG.trials;
    summary.rejection_rate = 100 * summary.epochs_rejected / n_epochs_before_rejection;
    summary.rejected_trial_indices = rejected_indices;
    fprintf('Rejected %d epochs (%.1f%%), retained %d epochs\n', ...
            summary.epochs_rejected, summary.rejection_rate, summary.epochs_retained);

    % Count epochs per condition
    epoch_counts = struct();
    for c = 1:length(CONDITION_EVENTS)
        cond = CONDITION_EVENTS{c};
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
        epoch_counts.(cond) = length(epoch_idx);
    end
    summary.epochs_per_condition = epoch_counts;

    % Save final output
    fprintf('Saving final output to %s\n', output_filename);
    EEG = pop_saveset(EEG, 'filename', output_filename, 'filepath', OUTPUT_DIR);

    % Write summary text file
    write_preprocessing_summary(summary_path, summary, CONDITION_EVENTS);
    fprintf('Wrote preprocessing summary to %s\n', summary_filename);
end

%% Helper function: Filter events to keep only those with CNAM containing filter string
function EEG = filter_events_by_cnam(EEG, cnam_filter)
    % Keep only events where mffkey_CNAM contains the filter string
    keep_idx = [];
    for e = 1:length(EEG.event)
        evt = EEG.event(e);
        if isfield(evt, 'mffkey_CNAM')
            cnam = evt.mffkey_CNAM;
            if ischar(cnam) && contains(cnam, cnam_filter)
                keep_idx = [keep_idx, e];
            end
        end
    end
    EEG.event = EEG.event(keep_idx);
end

%% Helper function: Relabel events based on condition type and trial type
function [EEG, relabel_counts] = relabel_events(EEG, basename)
    relabel_counts = struct();

    for e = 1:length(EEG.event)
        evt = EEG.event(e);
        evtype = evt.type;

        % Skip non-trial events
        if ~any(strcmp(evtype, {'STND', 'TON1', 'DEVT', 'OMIS'}))
            continue;
        end

        % Get condition info from event keys
        ctyp = '';  % condition type: single/double
        isi_val = 0;
        ttyp = '';  % trial type: standard/omission/deviant

        if isfield(evt, 'mffkey_CTYP')
            ctyp = evt.mffkey_CTYP;
            if isnumeric(ctyp), ctyp = ''; end
        end
        if isfield(evt, 'mffkey_ISI')
            isi_raw = evt.mffkey_ISI;
            if ischar(isi_raw)
                isi_val = round(str2double(isi_raw) * 1000);
            else
                isi_val = round(isi_raw * 1000);
            end
        end
        if isfield(evt, 'mffkey_TTYP')
            ttyp = evt.mffkey_TTYP;
            if isnumeric(ttyp), ttyp = ''; end
        end

        % Build new event label
        new_label = '';

        if strcmp(ctyp, 'single')
            if strcmp(evtype, 'STND') && strcmp(ttyp, 'standard')
                new_label = 'SNGL_STD';
            elseif strcmp(evtype, 'OMIS')
                new_label = 'SNGL_OMI';
            elseif strcmp(evtype, 'DEVT')
                new_label = 'SNGL_DEV';
            end

        elseif strcmp(ctyp, 'double')
            if strcmp(evtype, 'TON1')
                if isi_val == 100
                    prefix = 'DBL100';
                elseif isi_val == 250
                    prefix = 'DBL250';
                else
                    prefix = sprintf('DBL%d', isi_val);
                end

                if strcmp(ttyp, 'standard')
                    new_label = [prefix '_STD'];
                elseif strcmp(ttyp, 'omission')
                    new_label = [prefix '_OMI'];
                elseif strcmp(ttyp, 'deviant')
                    new_label = [prefix '_DEV'];
                end
            end
        end

        if ~isempty(new_label)
            EEG.event(e).type = new_label;
            if isfield(relabel_counts, new_label)
                relabel_counts.(new_label) = relabel_counts.(new_label) + 1;
            else
                relabel_counts.(new_label) = 1;
            end
        end
    end

    % Print summary
    fprintf('Event relabeling summary for %s:\n', basename);
    labels = fieldnames(relabel_counts);
    for i = 1:length(labels)
        fprintf('  %s: %d\n', labels{i}, relabel_counts.(labels{i}));
    end
end

%% Helper function: Detect bad channels using statistical measures
function [bad_chans, bad_chan_labels, bad_chan_details] = detect_bad_channels(EEG, var_thresh, corr_thresh, kurt_thresh)
    % Fast channel rejection based on variance, correlation, and kurtosis
    % Returns indices of bad channels, their labels, and details about why each was rejected

    data = EEG.data;
    nchan = size(data, 1);

    bad_variance = [];
    bad_correlation = [];
    bad_kurtosis = [];

    % 1. Variance-based detection (channels with abnormal variance)
    chan_var = var(data, 0, 2);
    z_var = (chan_var - median(chan_var)) / mad(chan_var, 1);
    bad_variance = find(abs(z_var) > var_thresh);

    % 2. Correlation-based detection (channels poorly correlated with neighbors)
    if ~isempty(EEG.chanlocs) && isfield(EEG.chanlocs, 'X')
        % Get channel positions
        X = [EEG.chanlocs.X]';
        Y = [EEG.chanlocs.Y]';
        Z = [EEG.chanlocs.Z]';

        % Only proceed if we have valid coordinates
        if ~any(isnan(X)) && ~any(isnan(Y)) && ~any(isnan(Z))
            pos = [X, Y, Z];

            % Compute pairwise distances
            dist_mat = squareform(pdist(pos));

            % For each channel, find average correlation with 4 nearest neighbors
            avg_corr = zeros(nchan, 1);
            for ch = 1:nchan
                [~, sorted_idx] = sort(dist_mat(ch, :));
                neighbors = sorted_idx(2:min(5, nchan));  % 4 nearest neighbors

                corrs = zeros(length(neighbors), 1);
                for n = 1:length(neighbors)
                    r = corrcoef(data(ch, :), data(neighbors(n), :));
                    corrs(n) = r(1, 2);
                end
                avg_corr(ch) = mean(corrs);
            end

            bad_correlation = find(avg_corr < corr_thresh);
        end
    end

    % 3. Kurtosis-based detection (channels with abnormal distribution)
    chan_kurt = kurtosis(data, 1, 2);
    z_kurt = (chan_kurt - median(chan_kurt)) / mad(chan_kurt, 1);
    bad_kurtosis = find(abs(z_kurt) > kurt_thresh);

    % Combine all bad channels
    bad_chans = unique([bad_variance; bad_correlation; bad_kurtosis]);

    % Get labels
    if ~isempty(bad_chans)
        bad_chan_labels = {EEG.chanlocs(bad_chans).labels};
    else
        bad_chan_labels = {};
    end

    % Build details struct
    bad_chan_details = struct();
    if ~isempty(bad_variance)
        bad_chan_details.high_variance = {EEG.chanlocs(bad_variance).labels};
    else
        bad_chan_details.high_variance = {};
    end
    if ~isempty(bad_correlation)
        bad_chan_details.low_correlation = {EEG.chanlocs(bad_correlation).labels};
    else
        bad_chan_details.low_correlation = {};
    end
    if ~isempty(bad_kurtosis)
        bad_chan_details.high_kurtosis = {EEG.chanlocs(bad_kurtosis).labels};
    else
        bad_chan_details.high_kurtosis = {};
    end

    % Print details
    if ~isempty(bad_variance)
        fprintf('  High variance channels: %s\n', strjoin({EEG.chanlocs(bad_variance).labels}, ', '));
    end
    if ~isempty(bad_correlation)
        fprintf('  Low correlation channels: %s\n', strjoin({EEG.chanlocs(bad_correlation).labels}, ', '));
    end
    if ~isempty(bad_kurtosis)
        fprintf('  High kurtosis channels: %s\n', strjoin({EEG.chanlocs(bad_kurtosis).labels}, ', '));
    end
end

%% Helper function: Write preprocessing summary to text file
function write_preprocessing_summary(filepath, summary, condition_events)
    fid = fopen(filepath, 'w');
    if fid == -1
        warning('Could not create summary file: %s', filepath);
        return;
    end

    fprintf(fid, '==============================================\n');
    fprintf(fid, 'PREPROCESSING SUMMARY\n');
    fprintf(fid, '==============================================\n\n');

    fprintf(fid, 'File: %s\n', summary.filename);
    fprintf(fid, 'Processing date: %s\n', summary.processing_date);
    fprintf(fid, 'CNAM filter: %s\n\n', summary.cnam_filter);

    fprintf(fid, '--- EVENT FILTERING ---\n');
    if isfield(summary, 'original_events')
        fprintf(fid, 'Original events: %d\n', summary.original_events);
    end
    if isfield(summary, 'events_after_cnam_filter')
        fprintf(fid, 'Events after CNAM filter: %d\n', summary.events_after_cnam_filter);
    end
    fprintf(fid, '\n');

    fprintf(fid, '--- CHANNEL PROCESSING ---\n');
    if isfield(summary, 'original_channels')
        fprintf(fid, 'Original channels: %d\n', summary.original_channels);
    end
    if isfield(summary, 'channels_kept')
        fprintf(fid, 'Channels kept (after selection): %d\n', summary.channels_kept);
    end
    if isfield(summary, 'bad_channels_detected')
        fprintf(fid, 'Bad channels detected: %d\n', summary.bad_channels_detected);
    end
    if isfield(summary, 'bad_channel_labels') && ~isempty(summary.bad_channel_labels)
        fprintf(fid, 'Bad channel labels: %s\n', strjoin(summary.bad_channel_labels, ', '));
    end
    if isfield(summary, 'bad_channel_details')
        if ~isempty(summary.bad_channel_details.high_variance)
            fprintf(fid, '  - High variance: %s\n', strjoin(summary.bad_channel_details.high_variance, ', '));
        end
        if ~isempty(summary.bad_channel_details.low_correlation)
            fprintf(fid, '  - Low correlation: %s\n', strjoin(summary.bad_channel_details.low_correlation, ', '));
        end
        if ~isempty(summary.bad_channel_details.high_kurtosis)
            fprintf(fid, '  - High kurtosis: %s\n', strjoin(summary.bad_channel_details.high_kurtosis, ', '));
        end
    end
    if isfield(summary, 'channels_interpolated')
        fprintf(fid, 'Channels interpolated: %d\n', summary.channels_interpolated);
    end
    fprintf(fid, '\n');

    fprintf(fid, '--- EPOCH PROCESSING ---\n');
    if isfield(summary, 'epochs_created')
        fprintf(fid, 'Epochs created: %d\n', summary.epochs_created);
    end
    if isfield(summary, 'epochs_rejected')
        fprintf(fid, 'Epochs rejected: %d\n', summary.epochs_rejected);
    end
    if isfield(summary, 'epochs_retained')
        fprintf(fid, 'Epochs retained: %d\n', summary.epochs_retained);
    end
    if isfield(summary, 'rejection_rate')
        fprintf(fid, 'Rejection rate: %.1f%%\n', summary.rejection_rate);
    end
    fprintf(fid, '\n');

    fprintf(fid, '--- EPOCHS PER CONDITION ---\n');
    if isfield(summary, 'epochs_per_condition')
        for c = 1:length(condition_events)
            cond = condition_events{c};
            if isfield(summary.epochs_per_condition, cond)
                fprintf(fid, '%s: %d\n', cond, summary.epochs_per_condition.(cond));
            else
                fprintf(fid, '%s: 0\n', cond);
            end
        end
    end
    fprintf(fid, '\n');

    fprintf(fid, '==============================================\n');

    fclose(fid);
end

%% Helper function: Save rejection plot showing which trials were rejected
function save_rejection_plot(rejected_mask, basename, filepath)
    n_trials = length(rejected_mask);
    rejected_indices = find(rejected_mask);
    n_rejected = length(rejected_indices);

    fig = figure('Position', [100 100 1000 400], 'Color', 'w', 'Visible', 'off');

    % Set black text color for this figure
    set(fig, 'DefaultAxesXColor', 'k', 'DefaultAxesYColor', 'k', 'DefaultTextColor', 'k');

    % Subplot 1: Trial rejection across experiment
    subplot(1, 2, 1);
    set(gca, 'Color', 'w', 'XColor', 'k', 'YColor', 'k');
    hold on;

    % Plot all trials as green dots (retained)
    retained_indices = find(~rejected_mask);
    scatter(retained_indices, ones(size(retained_indices)), 20, [0.2 0.6 0.2], 'filled', 'MarkerFaceAlpha', 0.5);

    % Plot rejected trials as red dots
    if ~isempty(rejected_indices)
        scatter(rejected_indices, ones(size(rejected_indices)), 40, [0.8 0.2 0.2], 'filled');
    end

    xlabel('Trial number (chronological order)');
    ylabel('');
    yticks([]);
    xlim([0 n_trials+1]);
    title(sprintf('%s: Trial Rejections (%d/%d = %.1f%%)', ...
        basename, n_rejected, n_trials, 100*n_rejected/n_trials), 'Interpreter', 'none');
    legend({'Retained', 'Rejected'}, 'Location', 'best');
    hold off;

    % Subplot 2: Histogram of rejections across experiment (binned)
    subplot(1, 2, 2);
    set(gca, 'Color', 'w', 'XColor', 'k', 'YColor', 'k');

    n_bins = 10;
    bin_edges = linspace(0, n_trials, n_bins+1);
    bin_centers = (bin_edges(1:end-1) + bin_edges(2:end)) / 2;

    % Count rejections per bin
    rejections_per_bin = histcounts(rejected_indices, bin_edges);
    trials_per_bin = histcounts(1:n_trials, bin_edges);
    rejection_rate_per_bin = 100 * rejections_per_bin ./ trials_per_bin;

    bar(bin_centers, rejection_rate_per_bin, 'FaceColor', [0.8 0.2 0.2], 'EdgeColor', 'none');
    xlabel('Trial number (binned)');
    ylabel('Rejection rate (%)');
    xlim([0 n_trials]);
    ylim([0 max(100, max(rejection_rate_per_bin)*1.1)]);
    title('Rejection Rate Across Experiment');

    % Add horizontal line for overall rejection rate
    yline(100*n_rejected/n_trials, 'k--', 'LineWidth', 1.5);

    sgtitle(sprintf('%s - Artifact Rejection Summary', basename), 'Interpreter', 'none', 'FontSize', 12);

    saveas(fig, filepath);
    close(fig);
end
