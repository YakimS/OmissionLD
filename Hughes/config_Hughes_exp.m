% config_Hughes_exp.m
% Generates experiment configuration and metadata for Hughes et al. (2001) omission oddball experiment
% "Responses of Human Auditory Association Cortex to the Omission of an Expected Acoustic Event"
%
% This script creates all experiment parameters, trial sequences, and stimuli,
% then saves everything to a .mat file for use by run_Hughes_exp.m
%
% Usage:
%   1. Modify parameters below
%   2. Run this script to generate and preview the configuration
%   3. Configuration is saved to OUTPUT_FILE
%   4. Use the saved .mat file with run_Hughes_exp.m

clear all;
close all;
clc;
%% ============= OUTPUT FILE NAME ============='
OUTPUT_FILE = 'Config-test_timingTest.mat';   %Set to false for automatic datetime naming: exp_config_YYYYMMDD_HHMMSS.mat

%% ============= EXPERIMENT PARAMETERS =============
% Stimulus parameters
exp_meta.STIM_DURATION = 0.1;      % Stimulus duration in seconds
exp_meta.STIM_FREQ = 1000;            % Tone frequency in Hz
exp_meta.SAMPLE_RATE = 44100;        % Audio sampling rate in Hz
exp_meta.ISI_MIN = 1;            % Minimum ISI in seconds
exp_meta.ISI_MAX = 1;            % Maximum ISI in seconds
exp_meta.ISI = (exp_meta.ISI_MIN + exp_meta.ISI_MAX) / 2;  % Mean ISI
exp_meta.WITHIN_PAIR_ISI = [0.5];    % ISI values within tone pairs in seconds (array)

% Trial parameters
exp_meta.N_TRIALS = 100;             % Total trials per procedure
exp_meta.OMISSION_FREQ = 0.001;        % Omission rate (proportion)
exp_meta.N_STANDARD_BEGIN = 10;      % Number of standard trials at beginning
exp_meta.MIN_STANDARD_BETWEEN = 2;   % Minimum standard trials between omissions

% Display parameters
exp_meta.SHOW_TRIAL_OUTPUT = true;   % Show trial-by-trial output
exp_meta.WAIT_BETWEEN_BLOCKS = 5;    % Wait time between blocks in seconds

%% Generate output filename
if OUTPUT_FILE == false
    timestamp = datestr(now, 'yyyymmdd_HHMMSS');
    OUTPUT_FILE = sprintf('exp_config_%s.mat', timestamp);
end

%% Calculate derived parameters
exp_meta.n_omissions = round(exp_meta.N_TRIALS * exp_meta.OMISSION_FREQ);
exp_meta.n_standards = exp_meta.N_TRIALS - exp_meta.n_omissions;
isi_values = linspace(exp_meta.ISI_MIN, exp_meta.ISI_MAX, exp_meta.N_TRIALS);

fprintf('\n=== Generating Audio Stimuli ===\n');
% Create sine wave beep
t = 0:1/exp_meta.SAMPLE_RATE:exp_meta.STIM_DURATION;
beep_sound = sin(2*pi*exp_meta.STIM_FREQ*t);

% Add onset/offset ramps to avoid clicks
ramp_duration = 0.005;  % Ramp duration in seconds
ramp_samples = round(ramp_duration * exp_meta.SAMPLE_RATE);
ramp = linspace(0, 1, ramp_samples);
beep_sound(1:ramp_samples) = beep_sound(1:ramp_samples) .* ramp;
beep_sound(end-ramp_samples+1:end) = beep_sound(end-ramp_samples+1:end) .* fliplr(ramp);

% Store beep sound in metadata
exp_meta.beep_sound = beep_sound;
fprintf('Audio stimulus generated: %d Hz, %.3f sec, %d Hz sampling\n', ...
        exp_meta.STIM_FREQ, exp_meta.STIM_DURATION, exp_meta.SAMPLE_RATE);

%% Create conditions and randomize order
n_isi_conditions = length(exp_meta.WITHIN_PAIR_ISI);
exp_meta.n_total_blocks = n_isi_conditions * 2;  
conditions = [];
for isi_idx = 1:n_isi_conditions
    conditions = [conditions; 1, exp_meta.WITHIN_PAIR_ISI(isi_idx), isi_idx];  % Procedure 1
    conditions = [conditions; 2, exp_meta.WITHIN_PAIR_ISI(isi_idx), isi_idx];  % Procedure 2
end
rand_order = randperm(exp_meta.n_total_blocks);
exp_meta.conditions = conditions(rand_order, :);

%% Create trial sequences for all conditions
try
    % Create unique trial sequence for each block
    trial_sequences = cell(exp_meta.n_total_blocks, 1);
    isi_arrays = cell(exp_meta.n_total_blocks, 1);

    for block = 1:exp_meta.n_total_blocks
        trial_sequences{block} = create_omission_sequence(exp_meta.N_TRIALS, exp_meta.n_omissions, ...
                                                           exp_meta.N_STANDARD_BEGIN, exp_meta.MIN_STANDARD_BETWEEN);
        isi_arrays{block} = isi_values(randperm(exp_meta.N_TRIALS));  % Randomized ISI for this block
    end
    exp_meta.trial_sequences = trial_sequences;
    exp_meta.isi_arrays = isi_arrays;

    % Calculate and display duration estimates
    fprintf('\n=== Duration Estimates ===\n');
    total_experiment_duration = 0;
    for block = 1:exp_meta.n_total_blocks
        proc_type = exp_meta.conditions(block, 1);
        isi_val = exp_meta.conditions(block, 2);
        isi_array_block = isi_arrays{block};
        proc_name = {'Single Tones', 'Tone Pairs'};

        block_duration = sum(exp_meta.STIM_DURATION + isi_array_block);

        total_experiment_duration = total_experiment_duration + block_duration;

        fprintf('Block %d - Procedure %d (%s, ISI=%.3fs): %.1f sec (%.2f min)\n', ...
                block, proc_type, proc_name{proc_type}, isi_val, ...
                block_duration, block_duration/60);
    end

    fprintf('TOTAL EXPERIMENT DURATION: %.1f sec (%.2f min)\n', total_experiment_duration, total_experiment_duration/60);
    fprintf('=========================\n');
    exp_meta.total_duration = total_experiment_duration;

catch ME
    error('Failed to create trial sequences: %s', ME.message);
end

%% Save configuration
exp_meta.created_date = datestr(now);
exp_meta.output_file = OUTPUT_FILE;
save(OUTPUT_FILE);

%% Helper function to create omission sequence with constraints
function trial_sequence = create_omission_sequence(n_trials, n_omissions, n_standard_begin, min_standard_between)
    % Initialize all trials as standards
    trial_sequence = ones(1, n_trials);

    % Determine valid positions for omissions (after initial standards)
    valid_start = n_standard_begin + 1;

    % Place omissions with minimum spacing constraint
    omission_positions = [];
    available_positions = valid_start:n_trials;

    for i = 1:n_omissions
        if isempty(available_positions)
            error('Cannot place %d omissions with constraints. Try fewer omissions or more trials.', n_omissions);
        end

        % Randomly select from available positions
        idx = randi(length(available_positions));
        pos = available_positions(idx);
        omission_positions(end+1) = pos;

        % Remove this position and surrounding positions from available list
        exclude_range = max(1, pos-min_standard_between):min(n_trials, pos+min_standard_between);
        available_positions = setdiff(available_positions, exclude_range);
    end

    % Set omissions in sequence (0 = omission, 1 = standard)
    trial_sequence(omission_positions) = 0;
end
