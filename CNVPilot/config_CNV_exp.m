% config_CNV_exp.m
% Generates experiment configuration for the auditory CNV (Contingent Negative
% Variation) experiment (Phase 1).
%
% Paradigm (after Gauthier et al. 1985; see CNV_paradigm_summary_and_sleep_mapping.docx):
%   Each trial is a couple: a warning stimulus S1 -> (S1-S2 interval, ISI) -> S2.
%   S1 is ACOUSTICALLY IDENTICAL across both couples, so the couple type is unknown
%   to the subject until S2 arrives. Two couples:
%     'same' : S1 -> S1 (same tone repeated)     -> TARGET, respond ASAP with eyes
%     'diff' : S1 -> S2 (a different tone)        -> no response
%   The anticipatory CNV builds during the S1-S2 interval (common to both couples);
%   the couples diverge only at/after S2.
%
% Response is by EYE MOVEMENT only, analysed offline from EOG. Nothing is collected
% online. The subject watches a silent nature movie (played externally) throughout.
%
% Usage:
%   1. Modify parameters below.
%   2. Run this script to generate the configuration .mat file.
%   3. Validate it with analyze_CNV_config('<file>.mat').
%   4. Run the experiment with run_CNV_exp.m.

clear all;
close all;
clc;

%% ============= OUTPUT FILE NAME =============
OUTPUT_FILE = 'Config-CNVPilot_100trialsX3ISIs.mat';  % Set to false for auto datetime naming

%% ============= EXPERIMENT PARAMETERS =============
% --- Stimulus parameters ---
exp_meta.STIM_DURATION = 0.2;        % Tone duration in seconds (S1 and S2), <=200 ms
exp_meta.SAMPLE_RATE   = 44100;      % Audio sampling rate in Hz

exp_meta.S1_FREQ    = 2000;          % S1 (warning) tone frequency (Hz)
exp_meta.S2_FREQ    = 500;           % S2 (different) tone frequency (Hz)
exp_meta.SWAP_FREQS = false;         % true -> swap S1_FREQ and S2_FREQ (per-subject CB)

% --- Timing parameters ---
exp_meta.ISI_VALUES      = [0.8 1 1.2];  % S1-S2 intervals (s); held constant within a block
exp_meta.N_BLOCKS_PER_ISI = 2;         % Blocks per ISI value

% Rest interval AFTER S2 (offset) before next S1. Jittered to keep S1 onset
% temporally unpredictable (protects the pre-S1 baseline from anticipatory drift).
exp_meta.REST_MEAN   = 6.0;          % Mean rest interval (s)
exp_meta.REST_JITTER = 3.0;          % Half-range of rest jitter (s) -> rest in [REST_MEAN-REST_JITTER, REST_MEAN+REST_JITTER] s

% --- Trial parameters ---
% Two couples: 'same' (target) and 'diff'. Choose the number of test trials PER ISI value
% (block type) and the percentage that are 'diff'; the rest are 'same'/target. Each ISI's
% trials are split as evenly as possible across its N_BLOCKS_PER_ISI blocks. No two 'diff'
% trials ever occur back to back (see create_couple_sequence); 'same' trials may repeat.
exp_meta.N_TRIALS     = 100;          % test trials PER ISI value (block type)
exp_meta.PERCENT_DIFF = 20;          % percent of test trials that are 'diff' (rest = 'same')
exp_meta.N_WARMUP     = 5;           % Warmup 'same'/target trials at each block start (excluded)

% --- Display / estimate parameters ---
exp_meta.SHOW_TRIAL_OUTPUT = true;
exp_meta.NOMINAL_BREAK_SEC = 20;     % Nominal self-paced break length, for duration ESTIMATE only

%% Apply per-subject frequency counterbalancing (once, for the whole session)
if exp_meta.SWAP_FREQS
    [exp_meta.S1_FREQ, exp_meta.S2_FREQ] = deal(exp_meta.S2_FREQ, exp_meta.S1_FREQ);
end

%% Generate output filename
if OUTPUT_FILE == false
    timestamp = datestr(now, 'yyyymmdd_HHMMSS');
    OUTPUT_FILE = sprintf('Config-CNVPilot_%s.mat', timestamp);
end

%% Derived counts + block layout + per-ISI trial split
couples = {'same', 'diff'};
n_isi          = numel(exp_meta.ISI_VALUES);
n_total_blocks = exp_meta.N_BLOCKS_PER_ISI * n_isi;

% Block ISI order (shuffled so no two same-ISI blocks are adjacent).
block_isi_pool = repelem(exp_meta.ISI_VALUES, exp_meta.N_BLOCKS_PER_ISI);  % e.g. [0.8 0.8 1 1 1.2 1.2]
block_isi = block_isi_pool;
for attempt = 1:1000
    order = block_isi_pool(randperm(numel(block_isi_pool)));
    if all(diff(order) ~= 0)   % no two same-ISI blocks adjacent
        block_isi = order; break;
    end
    block_isi = order;         % keep last attempt if guard never satisfied
end

% N_TRIALS is PER ISI value (block type): each ISI gets N_TRIALS test trials, split
% PERCENT_DIFF 'diff' / rest 'same', and distributed across that ISI's blocks.
block_total    = zeros(1, n_total_blocks);
diff_per_block = zeros(1, n_total_blocks);
same_per_block = zeros(1, n_total_blocks);
for iv = 1:n_isi
    idx = find(block_isi == exp_meta.ISI_VALUES(iv));       % this ISI's blocks
    tot = distribute_even(exp_meta.N_TRIALS, numel(idx));
    nd  = distribute_even(round(exp_meta.N_TRIALS * exp_meta.PERCENT_DIFF / 100), numel(idx));
    block_total(idx)    = tot;
    diff_per_block(idx) = nd;
    same_per_block(idx) = tot - nd;
end

% Feasibility of the "no two consecutive diff" rule: need enough 'same' to separate them.
if any(diff_per_block > same_per_block + 1)
    error('CNV:config', ['PERCENT_DIFF too high for N_TRIALS/blocks: a block would need ' ...
        'diff > same+1, so two diff trials could not be kept apart. Lower PERCENT_DIFF.']);
end

exp_meta.couples        = couples;
exp_meta.block_isi      = block_isi;
exp_meta.block_total    = block_total;
exp_meta.diff_per_block = diff_per_block;
exp_meta.same_per_block = same_per_block;
exp_meta.n_diff_total   = sum(diff_per_block);
exp_meta.n_same_total   = sum(same_per_block);
exp_meta.n_total_blocks = n_total_blocks;

%% Generate audio stimuli (two ramped sine tones: S1_FREQ and S2_FREQ)
fprintf('\n=== Generating Audio Stimuli ===\n');

t = 0:1/exp_meta.SAMPLE_RATE:exp_meta.STIM_DURATION;
beep_S1 = sin(2*pi*exp_meta.S1_FREQ*t);
beep_S2 = sin(2*pi*exp_meta.S2_FREQ*t);

% Onset/offset ramps to avoid clicks
ramp_duration = 0.005;
ramp_samples  = round(ramp_duration * exp_meta.SAMPLE_RATE);
ramp = linspace(0, 1, ramp_samples);

beep_S1(1:ramp_samples)          = beep_S1(1:ramp_samples) .* ramp;
beep_S1(end-ramp_samples+1:end)  = beep_S1(end-ramp_samples+1:end) .* fliplr(ramp);
beep_S2(1:ramp_samples)          = beep_S2(1:ramp_samples) .* ramp;
beep_S2(end-ramp_samples+1:end)  = beep_S2(end-ramp_samples+1:end) .* fliplr(ramp);

exp_meta.beep_S1 = beep_S1;
exp_meta.beep_S2 = beep_S2;

fprintf('S1 = %d Hz (warning), S2 = %d Hz (different), %.3f s, %d Hz sampling\n', ...
        exp_meta.S1_FREQ, exp_meta.S2_FREQ, exp_meta.STIM_DURATION, exp_meta.SAMPLE_RATE);

%% Block layout
fprintf('\n=== Blocks ===\n');
for b = 1:n_total_blocks
    fprintf('Block %d: ISI = %.0f ms | %d test (%d same, %d diff)\n', b, block_isi(b)*1000, ...
            block_total(b), same_per_block(b), diff_per_block(b));
end

%% Generate events table
fprintf('\n=== Generating Events Table ===\n');

col_names = {'block','trial','event_in_trial','couple','trial_kind','is_target', ...
             'isi','rest','event_type','freq','s1_freq','s2_freq','wait_after', ...
             'trigger','condition_name','display_str','ctyp'};
all_events = cell(0, numel(col_names));
event_row = 1;

for block = 1:n_total_blocks
    isi_val  = block_isi(block);
    cnam     = sprintf('CNV_ISI%.0fms', isi_val*1000);

    % Couple sequence for this block: N_WARMUP 'same' + interleaved test trials
    % (same/diff) with NO two 'diff' back to back.
    [couple_seq, kind_seq] = create_couple_sequence(same_per_block(block), ...
                                diff_per_block(block), exp_meta.N_WARMUP);
    n_trials_b = numel(couple_seq);

    % Per-trial rest jitter: evenly spaced then shuffled (deterministic mean)
    rest_jit  = linspace(-exp_meta.REST_JITTER, exp_meta.REST_JITTER, n_trials_b);
    rest_jit  = rest_jit(randperm(n_trials_b));
    rest_vals = exp_meta.REST_MEAN + rest_jit;

    for tr = 1:n_trials_b
        couple    = couple_seq{tr};
        kind      = kind_seq{tr};
        rest      = rest_vals(tr);
        is_target = double(strcmp(couple, 'same'));

        % ---- Event 1: S1 (warning tone, identical for both couples) ----
        all_events(event_row, :) = { ...
            block, tr, 1, couple, kind, is_target, ...
            isi_val, rest, 'tone', exp_meta.S1_FREQ, ...
            exp_meta.S1_FREQ, s2_freq_of(couple, exp_meta), ...
            isi_val - exp_meta.STIM_DURATION, ...          % wait_after: S2 onset lands at ISI
            'WARN', cnam, sprintf('S1(%dHz)', exp_meta.S1_FREQ), 'cnv'};
        event_row = event_row + 1;

        % ---- Event 2: S2 (same tone or different tone) ----
        if strcmp(couple, 'same')
            s2_freq = exp_meta.S1_FREQ;  s2_trig = 'S2SM';
            s2_disp = sprintf('S2same(%dHz)', s2_freq);
        else   % 'diff'
            s2_freq = exp_meta.S2_FREQ;  s2_trig = 'S2DF';
            s2_disp = sprintf('S2diff(%dHz)', s2_freq);
        end
        s2_wait = rest;                                    % remaining time to next S1

        all_events(event_row, :) = { ...
            block, tr, 2, couple, kind, is_target, ...
            isi_val, rest, 'tone', s2_freq, ...
            exp_meta.S1_FREQ, s2_freq_of(couple, exp_meta), ...
            s2_wait, s2_trig, cnam, s2_disp, 'cnv'};
        event_row = event_row + 1;
    end
end

events_table = cell2table(all_events, 'VariableNames', col_names);
exp_meta.events_table = events_table;
fprintf('Generated %d events across %d blocks\n', height(events_table), n_total_blocks);

%% Duration estimate (stimulus time; self-paced breaks estimated separately)
fprintf('\n=== Duration Estimate ===\n');
stim_duration = 0;
for block = 1:n_total_blocks
    be = events_table(events_table.block == block, :);
    block_dur = sum(be.wait_after) + sum(strcmp(be.event_type, 'tone')) * exp_meta.STIM_DURATION;
    fprintf('Block %d (ISI %.0f ms): %.1f s (%.2f min)\n', block, block_isi(block)*1000, ...
            block_dur, block_dur/60);
    stim_duration = stim_duration + block_dur;
end
nominal_breaks = (n_total_blocks - 1) * exp_meta.NOMINAL_BREAK_SEC;
fixed_waits    = 2 + 2;   % pre-start and pre-end WaitSecs(2) in the runner
total_duration = stim_duration + nominal_breaks + fixed_waits;

fprintf('\nStimulus time (excl. breaks): %.1f s (%.2f min)\n', stim_duration, stim_duration/60);
fprintf('Nominal breaks (%d x %ds):    %.1f s\n', n_total_blocks-1, exp_meta.NOMINAL_BREAK_SEC, nominal_breaks);
fprintf('ESTIMATED TOTAL:              %.1f s (%.2f min)\n', total_duration, total_duration/60);

exp_meta.stim_duration_est = stim_duration;
exp_meta.total_duration    = total_duration;

%% Summary statistics
fprintf('\n=== Events Summary ===\n');
is_test = strcmp(events_table.trial_kind, 'test') & events_table.event_in_trial == 1;
is_warm = strcmp(events_table.trial_kind, 'warmup') & events_table.event_in_trial == 1;
n_same = sum(is_test & strcmp(events_table.couple,'same'));
n_diff = sum(is_test & strcmp(events_table.couple,'diff'));
fprintf('Test trials:   %d total = %d/ISI x %d ISIs (same=%d, diff=%d -> %.1f%% diff)\n', ...
        sum(is_test), exp_meta.N_TRIALS, n_isi, n_same, n_diff, 100*n_diff/max(1,sum(is_test)));
fprintf('Warmup trials: %d (all target/same)\n', sum(is_warm));
for iv = 1:n_isi
    isi_val = exp_meta.ISI_VALUES(iv);
    n_isi_trials = sum(is_test & events_table.isi == isi_val);
    fprintf('  ISI %.0f ms: %d test trials\n', isi_val*1000, n_isi_trials);
end

%% Parameters summary string (stored on the config, echoed into NetStation)
exp_meta.params_string = sprintf(['STIM_DURATION=%.3f, SAMPLE_RATE=%d, ' ...
    'S1_FREQ=%d, S2_FREQ=%d, SWAP_FREQS=%d, ' ...
    'ISI_VALUES=[%s], N_BLOCKS_PER_ISI=%d, ' ...
    'REST_MEAN=%.2f, REST_JITTER=%.2f, ' ...
    'N_TRIALS=%d, PERCENT_DIFF=%d, N_WARMUP=%d, ' ...
    'n_total_blocks=%d, n_same_total=%d, n_diff_total=%d'], ...
    exp_meta.STIM_DURATION, exp_meta.SAMPLE_RATE, ...
    exp_meta.S1_FREQ, exp_meta.S2_FREQ, exp_meta.SWAP_FREQS, ...
    num2str(exp_meta.ISI_VALUES), exp_meta.N_BLOCKS_PER_ISI, ...
    exp_meta.REST_MEAN, exp_meta.REST_JITTER, ...
    exp_meta.N_TRIALS, exp_meta.PERCENT_DIFF, exp_meta.N_WARMUP, ...
    exp_meta.n_total_blocks, exp_meta.n_same_total, exp_meta.n_diff_total);

%% Save configuration
exp_meta.created_date = datestr(now);
exp_meta.output_file  = OUTPUT_FILE;
save(OUTPUT_FILE, 'exp_meta');
fprintf('\n=== Configuration saved to %s ===\n', OUTPUT_FILE);
fprintf('Validate with: analyze_CNV_config(''%s'')\n', OUTPUT_FILE);

%% ===================== Helper functions =====================

function f = s2_freq_of(couple, exp_meta)
    % Nominal S2 frequency for a couple. Stored as the S2FQ metadata field on BOTH
    % events of a trial, so the couple is fully described from either row.
    if strcmp(couple, 'same')
        f = exp_meta.S1_FREQ;
    else   % 'diff'
        f = exp_meta.S2_FREQ;
    end
end

function v = distribute_even(total, k)
    % Split `total` into k integers that sum to `total`, as evenly as possible.
    % The remainder is spread over a random subset of blocks so no block is
    % systematically larger.
    base = floor(total / k);
    r    = total - base * k;
    v    = base * ones(1, k);
    if r > 0
        extra = randperm(k, r);   % r distinct blocks get one extra trial
        v(extra) = v(extra) + 1;
    end
end

function [couple_seq, kind_seq] = create_couple_sequence(n_same, n_diff, n_warmup)
    % Build one block's couple order (two couples: 'same'/target, 'diff'):
    %   - n_warmup leading 'same'/target trials (tagged 'warmup', excluded from analysis)
    %   - then n_same 'same' + n_diff 'diff' interleaved so that NO TWO 'diff' trials are
    %     adjacent (diff run length == 1). 'same' may repeat (it is the majority/target).
    %
    % Gap method: the n_same 'same' trials create n_same+1 gaps; dropping each 'diff' into
    % a DISTINCT gap guarantees every 'diff' is separated by at least one 'same'.
    if n_diff > n_same + 1
        error('CNV:seq', 'Cannot keep diff apart: n_diff=%d > n_same+1=%d.', n_diff, n_same + 1);
    end

    chosen = randperm(n_same + 1, n_diff);   % distinct gap indices for the diff trials
    is_diff_gap = false(1, n_same + 1);
    is_diff_gap(chosen) = true;

    test = cell(1, n_same + n_diff);
    idx  = 0;
    for g = 1:(n_same + 1)
        if is_diff_gap(g)
            idx = idx + 1;  test{idx} = 'diff';
        end
        if g <= n_same
            idx = idx + 1;  test{idx} = 'same';
        end
    end

    couple_seq = [repmat({'same'}, 1, n_warmup), test];
    kind_seq   = [repmat({'warmup'}, 1, n_warmup), repmat({'test'}, 1, numel(test))];
end
