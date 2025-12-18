% config_OmissionPilot.m
% Generates experiment configuration for OmissionPilot experiment
%
% Conditions:
%   - Each condition has multiple blocks
%   - Deviants: omission and deviant tone (randomly placed)
%
% Usage:
%   1. Modify parameters below
%   2. Run this script to generate configuration
%   3. Use saved .mat file with run_OmissionPilot.m

clear all;
close all;
clc;

%% ============= OUTPUT FILE NAME =============
OUTPUT_FILE = 'Config-DoubleOmission_iti-?.mat';  % Set to false for auto datetime naming

%% ============= EXPERIMENT PARAMETERS =============
% Stimulus parameters
exp_meta.STIM_DURATION = 0.05;       % Stimulus duration in seconds
exp_meta.SAMPLE_RATE = 44100;       % Audio sampling rate in Hz

% Tone frequencies (will be shuffled per block as standard/deviant)
exp_meta.FREQ_1 = 1000;             % First tone frequency (Hz) 
exp_meta.FREQ_2 = 1200;             % Second tone frequency (Hz) 

% Timing parameters
exp_meta.ITI_VALUES = [1.2 ];   % Inter-trial interval values (seconds)
exp_meta.ISI_VALUES = [0.1, 0.25];   % Inter-stimulus intervals for double conditions (seconds)
exp_meta.ITI_JITTER_MIN = -0.05;  % Minimum jitter in seconds (-50ms)
exp_meta.ITI_JITTER_MAX = 0.05;   % Maximum jitter in seconds (+50ms)
exp_meta.ITI_JITTER_MEAN = (exp_meta.ITI_JITTER_MIN + exp_meta.ITI_JITTER_MAX) / 2;  % Mean jitter (0)

% Trial parameters
exp_meta.N_TRIALS = 600;            % Total trials per condition (double(isi x iti) + single(iti))
exp_meta.BLOCKS_PER_CONDITION = 6;  % Number of blocks per condition
exp_meta.OMISSION_RATE = 0.1;       % Proportion of omission trials
exp_meta.DEVIANT_RATE = 0.1;        % Proportion of deviant tone trials
exp_meta.N_STANDARD_BEGIN =10;      % Standard trials at beginning of each block
exp_meta.MIN_STANDARD_BETWEEN = 3;  % Minimum standards between ANY deviants

% Display/timing parameters
exp_meta.SHOW_TRIAL_OUTPUT = true;
exp_meta.WAIT_BETWEEN_BLOCKS = 5;   % Seconds between blocks

%% Generate output filename
if OUTPUT_FILE == false
    timestamp = datestr(now, 'yyyymmdd_HHMMSS');
    OUTPUT_FILE = sprintf('Config-OmissionPilot_%s.mat', timestamp);
end

%% Calculate derived parameters
trials_per_block = floor(exp_meta.N_TRIALS / exp_meta.BLOCKS_PER_CONDITION);
n_omissions_per_block = round(trials_per_block * exp_meta.OMISSION_RATE);
n_deviants_per_block = round(trials_per_block * exp_meta.DEVIANT_RATE);
total_deviants_per_block = n_omissions_per_block + n_deviants_per_block;

% Store in exp_meta
exp_meta.trials_per_block = trials_per_block;
exp_meta.n_omissions_per_block = n_omissions_per_block;
exp_meta.n_deviants_per_block = n_deviants_per_block;

%% Validate constraints
fprintf('\n=== Validating Constraints ===\n');
available_positions = trials_per_block - exp_meta.N_STANDARD_BEGIN;
% Maximum deviants that can fit: floor((available + min_between) / (min_between + 1))
max_deviants_possible = floor((available_positions + exp_meta.MIN_STANDARD_BETWEEN) / (exp_meta.MIN_STANDARD_BETWEEN + 1));

fprintf('Trials per block: %d\n', trials_per_block);
fprintf('Omissions per block: %d\n', n_omissions_per_block);
fprintf('Deviant tones per block: %d\n', n_deviants_per_block);
fprintf('Total deviants per block: %d\n', total_deviants_per_block);
fprintf('Available positions (after %d standards): %d\n', exp_meta.N_STANDARD_BEGIN, available_positions);
fprintf('Max deviants possible (with spacing %d): %d\n', exp_meta.MIN_STANDARD_BETWEEN, max_deviants_possible);

if total_deviants_per_block > max_deviants_possible
    error(['Cannot place %d deviants with constraints.\n' ...
           'Max possible with spacing %d: %d deviants.\n' ...
           'Try: fewer deviants, more trials, smaller N_STANDARD_BEGIN, or smaller MIN_STANDARD_BETWEEN.'], ...
           total_deviants_per_block, exp_meta.MIN_STANDARD_BETWEEN, max_deviants_possible);
end
fprintf('Constraint check PASSED.\n');

%% Generate audio stimuli
fprintf('\n=== Generating Audio Stimuli ===\n');

% Create sine wave for frequency 1
t = 0:1/exp_meta.SAMPLE_RATE:exp_meta.STIM_DURATION;
beep_freq1 = sin(2*pi*exp_meta.FREQ_1*t);
beep_freq2 = sin(2*pi*exp_meta.FREQ_2*t);

% Add onset/offset ramps to avoid clicks
ramp_duration = 0.005;
ramp_samples = round(ramp_duration * exp_meta.SAMPLE_RATE);
ramp = linspace(0, 1, ramp_samples);

beep_freq1(1:ramp_samples) = beep_freq1(1:ramp_samples) .* ramp;
beep_freq1(end-ramp_samples+1:end) = beep_freq1(end-ramp_samples+1:end) .* fliplr(ramp);
beep_freq2(1:ramp_samples) = beep_freq2(1:ramp_samples) .* ramp;
beep_freq2(end-ramp_samples+1:end) = beep_freq2(end-ramp_samples+1:end) .* fliplr(ramp);

exp_meta.beep_freq1 = beep_freq1;
exp_meta.beep_freq2 = beep_freq2;

fprintf('Generated tones: %d Hz and %d Hz, %.3f sec, %d Hz sampling\n', ...
        exp_meta.FREQ_1, exp_meta.FREQ_2, exp_meta.STIM_DURATION, exp_meta.SAMPLE_RATE);

%% Build condition list
% Conditions: for each ITI, 1 single + 2 double (one per ISI)
fprintf('\n=== Building Conditions ===\n');

conditions = {};
cond_idx = 1;
for iti_idx = 1:length(exp_meta.ITI_VALUES)
    iti_val = exp_meta.ITI_VALUES(iti_idx);

    % Single condition
    conditions{cond_idx}.type = 'single';
    conditions{cond_idx}.iti = iti_val;
    conditions{cond_idx}.isi = NaN;
    conditions{cond_idx}.name = sprintf('Single_ITI%.0fms', iti_val*1000);
    fprintf('Condition %d: %s\n', cond_idx, conditions{cond_idx}.name);
    cond_idx = cond_idx + 1;

    % Double conditions (one per ISI value)
    for isi_idx = 1:length(exp_meta.ISI_VALUES)
        isi_val = exp_meta.ISI_VALUES(isi_idx);
        conditions{cond_idx}.type = 'double';
        conditions{cond_idx}.iti = iti_val;
        conditions{cond_idx}.isi = isi_val;
        conditions{cond_idx}.name = sprintf('Double_ITI%.0fms_ISI%.0fms', iti_val*1000, isi_val*1000);
        fprintf('Condition %d: %s\n', cond_idx, conditions{cond_idx}.name);
        cond_idx = cond_idx + 1;
    end
end
exp_meta.conditions = conditions;
exp_meta.n_conditions = length(conditions);

%% Build block list and shuffle
fprintf('\n=== Building Block List ===\n');

n_total_blocks = exp_meta.n_conditions * exp_meta.BLOCKS_PER_CONDITION;
exp_meta.n_total_blocks = n_total_blocks;

% Create block info: which condition and which block number within condition
block_info = [];
for cond = 1:exp_meta.n_conditions
    for blk = 1:exp_meta.BLOCKS_PER_CONDITION
        block_info = [block_info; cond, blk];
    end
end

% Shuffle block order
shuffle_order = randperm(n_total_blocks);
block_info = block_info(shuffle_order, :);
exp_meta.block_info = block_info;

fprintf('Total blocks: %d (shuffled)\n', n_total_blocks);

%% Assign standard/deviant frequencies per block (shuffled)
fprintf('\n=== Assigning Frequencies Per Block ===\n');

block_freq_assignment = zeros(n_total_blocks, 2);  % [standard_freq, deviant_freq]
for b = 1:n_total_blocks
    if rand() < 0.5
        block_freq_assignment(b, :) = [exp_meta.FREQ_1, exp_meta.FREQ_2];
    else
        block_freq_assignment(b, :) = [exp_meta.FREQ_2, exp_meta.FREQ_1];
    end
    fprintf('Block %d: Standard=%dHz, Deviant=%dHz\n', b, ...
            block_freq_assignment(b,1), block_freq_assignment(b,2));
end
exp_meta.block_freq_assignment = block_freq_assignment;

%% Generate events table
fprintf('\n=== Generating Events Table ===\n');

% Preallocate cell arrays for table
all_events = {};
event_row = 1;

for block = 1:n_total_blocks
    cond_idx = block_info(block, 1);
    block_within_cond = block_info(block, 2);
    cond = conditions{cond_idx};

    standard_freq = block_freq_assignment(block, 1);
    deviant_freq = block_freq_assignment(block, 2);

    % Generate trial types for this block
    trial_types = create_deviant_sequence(trials_per_block, n_omissions_per_block, ...
                                          n_deviants_per_block, exp_meta.N_STANDARD_BEGIN, ...
                                          exp_meta.MIN_STANDARD_BETWEEN);

    % Generate ITI jitter for each trial (evenly spaced, then shuffled for deterministic mean)
    iti_jitters = linspace(exp_meta.ITI_JITTER_MIN, exp_meta.ITI_JITTER_MAX, trials_per_block);
    iti_jitters = iti_jitters(randperm(trials_per_block));

    for trial = 1:trials_per_block
        trial_type = trial_types{trial};  % 'standard', 'omission', or 'deviant'
        iti_jitter = iti_jitters(trial);

        if strcmp(cond.type, 'single')
            % Single condition: one event per trial
            % Total trial time = ITI (onset-to-onset)
            base_iti = cond.iti + iti_jitter;

            switch trial_type
                case 'standard'
                    event_type = 'tone';
                    freq = standard_freq;
                    trigger = 'STND';
                    display_str = sprintf('STANDARD (%dHz)', freq);
                    wait_after = base_iti - exp_meta.STIM_DURATION;  % Account for tone duration
                case 'omission'
                    event_type = 'omission';
                    freq = NaN;
                    trigger = 'OMIS';
                    display_str = 'OMISSION';
                    wait_after = base_iti;  % No tone, so full ITI
                case 'deviant'
                    event_type = 'tone';
                    freq = deviant_freq;
                    trigger = 'DEVT';
                    display_str = sprintf('DEVIANT (%dHz)', freq);
                    wait_after = base_iti - exp_meta.STIM_DURATION;  % Account for tone duration
            end

            all_events{event_row, 1} = block;
            all_events{event_row, 2} = trial;
            all_events{event_row, 3} = 1;  % event_in_trial
            all_events{event_row, 4} = cond.type;
            all_events{event_row, 5} = cond.iti;
            all_events{event_row, 6} = cond.isi;
            all_events{event_row, 7} = trial_type;
            all_events{event_row, 8} = event_type;
            all_events{event_row, 9} = freq;
            all_events{event_row, 10} = wait_after;
            all_events{event_row, 11} = trigger;
            all_events{event_row, 12} = display_str;
            all_events{event_row, 13} = cond.name;
            event_row = event_row + 1;

        else
            % Double condition: two events per trial
            % Total trial time = ITI (onset-to-onset)
            % Tone1 at t=0, Tone2 at t=ISI, next trial at t=ITI
            base_iti = cond.iti + iti_jitter;

            % Event 1: First tone (always standard)
            % Wait after = ISI - STIM_DURATION (so tone2 onset is at ISI from tone1 onset)
            all_events{event_row, 1} = block;
            all_events{event_row, 2} = trial;
            all_events{event_row, 3} = 1;  % event_in_trial
            all_events{event_row, 4} = cond.type;
            all_events{event_row, 5} = cond.iti;
            all_events{event_row, 6} = cond.isi;
            all_events{event_row, 7} = trial_type;
            all_events{event_row, 8} = 'tone';
            all_events{event_row, 9} = standard_freq;
            all_events{event_row, 10} = cond.isi - exp_meta.STIM_DURATION;  % ISI minus tone1 duration
            all_events{event_row, 11} = 'TON1';
            all_events{event_row, 12} = sprintf('TONE1 (%dHz)', standard_freq);
            all_events{event_row, 13} = cond.name;
            event_row = event_row + 1;

            % Event 2: Second tone (standard, omission, or deviant)
            % Remaining time = ITI - ISI
            remaining_time = base_iti - cond.isi;

            switch trial_type
                case 'standard'
                    event_type = 'tone';
                    freq = standard_freq;
                    trigger = 'TON2';
                    display_str = sprintf('TONE2 (%dHz)', freq);
                    wait_after = remaining_time - exp_meta.STIM_DURATION;  % Account for tone duration
                case 'omission'
                    event_type = 'omission';
                    freq = NaN;
                    trigger = 'OMIS';
                    display_str = 'OMISSION';
                    wait_after = remaining_time;  % No tone, full remaining time
                case 'deviant'
                    event_type = 'tone';
                    freq = deviant_freq;
                    trigger = 'DEVT';
                    display_str = sprintf('DEVIANT (%dHz)', freq);
                    wait_after = remaining_time - exp_meta.STIM_DURATION;  % Account for tone duration
            end

            all_events{event_row, 1} = block;
            all_events{event_row, 2} = trial;
            all_events{event_row, 3} = 2;  % event_in_trial
            all_events{event_row, 4} = cond.type;
            all_events{event_row, 5} = cond.iti;
            all_events{event_row, 6} = cond.isi;
            all_events{event_row, 7} = trial_type;
            all_events{event_row, 8} = event_type;
            all_events{event_row, 9} = freq;
            all_events{event_row, 10} = wait_after;
            all_events{event_row, 11} = trigger;
            all_events{event_row, 12} = display_str;
            all_events{event_row, 13} = cond.name;
            event_row = event_row + 1;
        end
    end
end

% Convert to table
events_table = cell2table(all_events, 'VariableNames', ...
    {'block', 'trial', 'event_in_trial', 'condition_type', 'iti', 'isi', ...
     'trial_type', 'event_type', 'freq', 'wait_after', 'trigger', 'display_str', 'condition_name'});

exp_meta.events_table = events_table;
fprintf('Generated %d events across %d blocks\n', height(events_table), n_total_blocks);

%% Calculate duration estimates
fprintf('\n=== Duration Estimates ===\n');

total_duration = 0;
for block = 1:n_total_blocks
    block_events = events_table(events_table.block == block, :);
    block_duration = sum(block_events.wait_after) + ...
                     sum(strcmp(block_events.event_type, 'tone')) * exp_meta.STIM_DURATION;

    cond_idx = block_info(block, 1);
    cond_name = conditions{cond_idx}.name;

    fprintf('Block %d (%s): %.1f sec (%.2f min)\n', block, cond_name, block_duration, block_duration/60);
    total_duration = total_duration + block_duration;
end

% Add wait between blocks
total_duration = total_duration + (n_total_blocks - 1) * exp_meta.WAIT_BETWEEN_BLOCKS;

fprintf('\nTOTAL EXPERIMENT DURATION: %.1f sec (%.2f min)\n', total_duration, total_duration/60);
exp_meta.total_duration = total_duration;

%% Display summary statistics
fprintf('\n=== Events Summary ===\n');
fprintf('Standards: %d\n', sum(strcmp(events_table.trial_type, 'standard')));
fprintf('Omissions: %d\n', sum(strcmp(events_table.trial_type, 'omission')));
fprintf('Deviants: %d\n', sum(strcmp(events_table.trial_type, 'deviant')));

%% Create parameters summary string
exp_meta.params_string = sprintf(['STIM_DURATION=%.3f, SAMPLE_RATE=%d, ' ...
    'FREQ_1=%d, FREQ_2=%d, ' ...
    'ITI_VALUES=[%s], ISI_VALUES=[%s], ' ...
    'ITI_JITTER_MIN=%.3f, ITI_JITTER_MAX=%.3f, ' ...
    'N_TRIALS=%d, BLOCKS_PER_CONDITION=%d, ' ...
    'OMISSION_RATE=%.2f, DEVIANT_RATE=%.2f, ' ...
    'N_STANDARD_BEGIN=%d, MIN_STANDARD_BETWEEN=%d, ' ...
    'WAIT_BETWEEN_BLOCKS=%d, n_conditions=%d, n_total_blocks=%d, ' ...
    'trials_per_block=%d, n_omissions_per_block=%d, n_deviants_per_block=%d'], ...
    exp_meta.STIM_DURATION, exp_meta.SAMPLE_RATE, ...
    exp_meta.FREQ_1, exp_meta.FREQ_2, ...
    num2str(exp_meta.ITI_VALUES), num2str(exp_meta.ISI_VALUES), ...
    exp_meta.ITI_JITTER_MIN, exp_meta.ITI_JITTER_MAX, ...
    exp_meta.N_TRIALS, exp_meta.BLOCKS_PER_CONDITION, ...
    exp_meta.OMISSION_RATE, exp_meta.DEVIANT_RATE, ...
    exp_meta.N_STANDARD_BEGIN, exp_meta.MIN_STANDARD_BETWEEN, ...
    exp_meta.WAIT_BETWEEN_BLOCKS, exp_meta.n_conditions, exp_meta.n_total_blocks, ...
    exp_meta.trials_per_block, exp_meta.n_omissions_per_block, exp_meta.n_deviants_per_block);

%% Save configuration
exp_meta.created_date = datestr(now);
exp_meta.output_file = OUTPUT_FILE;
save(OUTPUT_FILE, 'exp_meta');
fprintf('\n=== Configuration saved to %s ===\n', OUTPUT_FILE);

%% Helper function to create deviant sequence with constraints
function trial_types = create_deviant_sequence(n_trials, n_omissions, n_deviants, n_standard_begin, min_between)
    % Deterministic algorithm:
    % 1. Start with n_standard_begin standards
    % 2. Place deviants with exactly min_between standards between each
    % 3. Randomly distribute remaining standards in the gaps
    % 4. Shuffle gaps to avoid consecutive gaps of same size

    total_deviants = n_omissions + n_deviants;

    % Calculate available positions and required minimum
    n_available = n_trials - n_standard_begin;
    min_standards_needed = (total_deviants - 1) * min_between;  % standards between deviants

    if total_deviants + min_standards_needed > n_available
        error('Cannot place all deviants with constraints. Need %d positions, have %d.', ...
              total_deviants + min_standards_needed, n_available);
    end

    % Extra standards to distribute randomly
    extra_standards = n_available - total_deviants - min_standards_needed;

    % Create gaps array: (total_deviants + 1) gaps
    % Gap 1: before first deviant (can be 0 or more)
    % Gaps 2 to total_deviants: between deviants (minimum min_between each)
    % Gap total_deviants+1: after last deviant (can be 0 or more)
    n_gaps = total_deviants + 1;
    gaps = zeros(1, n_gaps);

    % Set minimum spacing for gaps between deviants (indices 2 to n_gaps-1)
    gaps(2:end-1) = min_between;

    % Randomly distribute extra standards across all gaps
    for i = 1:extra_standards
        gap_idx = randi(n_gaps);
        gaps(gap_idx) = gaps(gap_idx) + 1;
    end

    % Shuffle the middle gaps (between deviants) to avoid consecutive same sizes
    middle_gaps = gaps(2:end-1);
    middle_gaps = shuffle_no_consecutive_repeats(middle_gaps);
    gaps(2:end-1) = middle_gaps;

    % Build trial sequence
    trial_types = {};

    % Add initial standards (n_standard_begin)
    for i = 1:n_standard_begin
        trial_types{end+1} = 'standard';
    end

    % Add gap before first deviant
    for i = 1:gaps(1)
        trial_types{end+1} = 'standard';
    end

    % Create deviant type assignments (shuffled)
    deviant_types = [repmat({'omission'}, 1, n_omissions), repmat({'deviant'}, 1, n_deviants)];
    deviant_types = deviant_types(randperm(length(deviant_types)));

    % Add deviants with gaps between them
    for d = 1:total_deviants
        trial_types{end+1} = deviant_types{d};

        % Add gap after this deviant (if not the last one, or if there's a trailing gap)
        if d < total_deviants
            for i = 1:gaps(d+1)
                trial_types{end+1} = 'standard';
            end
        else
            % Trailing gap after last deviant
            for i = 1:gaps(end)
                trial_types{end+1} = 'standard';
            end
        end
    end
end

%% Helper function to shuffle array avoiding consecutive repeats
function arr = shuffle_no_consecutive_repeats(arr)
    if length(arr) <= 1
        return;
    end

    % Try to rearrange to avoid consecutive repeats
    max_attempts = 1000;

    for attempt = 1:max_attempts
        % Shuffle the array
        arr = arr(randperm(length(arr)));

        % Check for consecutive repeats
        has_repeat = false;
        for i = 1:length(arr)-1
            if arr(i) == arr(i+1)
                has_repeat = true;
                break;
            end
        end

        if ~has_repeat
            return;
        end
    end

    % If random shuffle failed, use greedy algorithm
    arr = greedy_no_consecutive(arr);
end

%% Greedy algorithm to arrange without consecutive repeats
function result = greedy_no_consecutive(arr)
    % Count occurrences of each value
    unique_vals = unique(arr);
    counts = zeros(size(unique_vals));
    for i = 1:length(unique_vals)
        counts(i) = sum(arr == unique_vals(i));
    end

    result = zeros(size(arr));
    last_val = NaN;

    for pos = 1:length(arr)
        % Find valid candidates (not equal to last value, with count > 0)
        valid_mask = counts > 0;
        if ~isnan(last_val)
            last_idx = find(unique_vals == last_val, 1);
            if ~isempty(last_idx)
                valid_mask(last_idx) = false;
            end
        end

        if ~any(valid_mask)
            % No valid option - must repeat (unavoidable)
            valid_mask = counts > 0;
        end

        % Among valid options, pick the one with highest count (greedy)
        valid_counts = counts;
        valid_counts(~valid_mask) = -1;
        [~, best_idx] = max(valid_counts);

        result(pos) = unique_vals(best_idx);
        counts(best_idx) = counts(best_idx) - 1;
        last_val = result(pos);
    end
end
