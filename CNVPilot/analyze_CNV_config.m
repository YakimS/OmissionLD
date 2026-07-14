function analyze_CNV_config(config_file, show_plots)
% analyze_CNV_config  Validate + summarise a CNV experiment configuration.
%
%   analyze_CNV_config()                      % uses 'Config-CNVPilot.mat'
%   analyze_CNV_config('Config-CNVPilot.mat')
%   analyze_CNV_config(FILE, false)           % skip figures (headless)
%
% Read-only QA tool. Loads the config, prints trial counts and durations, and runs
% a battery of ASSERTIONS that catch timing/counting/metadata bugs before a subject
% is in the chair. Errors (nonzero exit) if any check fails.
%
% Two couples: 'same' (target) and 'diff'. No two 'diff' trials may be adjacent.

    if nargin < 1 || isempty(config_file); config_file = 'Config-CNVPilot.mat'; end
    if nargin < 2 || isempty(show_plots);  show_plots  = true;  end

    loaded  = load(config_file);
    exp     = loaded.exp_meta;
    events  = exp.events_table;
    STIM    = exp.STIM_DURATION;
    tol     = 1e-9;

    fprintf('\n=========================================================\n');
    fprintf(' CNV CONFIG ANALYZER  -  %s\n', config_file);
    fprintf('=========================================================\n');
    fprintf('S1=%d Hz  S2=%d Hz  STIM=%.0f ms  SR=%d Hz  SWAP_FREQS=%d\n', ...
            exp.S1_FREQ, exp.S2_FREQ, STIM*1000, exp.SAMPLE_RATE, exp.SWAP_FREQS);
    fprintf('ISIs=[%s] s  blocks=%d (%d per ISI)  rest=%.1f+/-%.1f s\n', ...
            num2str(exp.ISI_VALUES), exp.n_total_blocks, exp.N_BLOCKS_PER_ISI, ...
            exp.REST_MEAN, exp.REST_JITTER);
    fprintf('N_TRIALS=%d/ISI  PERCENT_DIFF=%d%%  N_WARMUP=%d\n', ...
            exp.N_TRIALS, exp.PERCENT_DIFF, exp.N_WARMUP);

    fails = {};

    %% ---- 1. Counts ----
    fprintf('\n--- Trial counts ---\n');
    s1 = events(events.event_in_trial == 1, :);                 % one row per trial
    test = s1(strcmp(s1.trial_kind, 'test'), :);
    warm = s1(strcmp(s1.trial_kind, 'warmup'), :);

    n_test = height(test);
    n_same = sum(strcmp(test.couple, 'same'));
    n_diff = sum(strcmp(test.couple, 'diff'));
    fprintf('Total test trials: %d (same=%d, diff=%d -> %.1f%% diff) | warmup: %d\n', ...
            n_test, n_same, n_diff, 100*n_diff/max(1,n_test), height(warm));

    fails = expect(fails, n_test == exp.N_TRIALS * numel(exp.ISI_VALUES), 'total test == N_TRIALS x nISI');
    fails = expect(fails, n_same == exp.n_same_total, 'total same == n_same_total');
    fails = expect(fails, n_diff == exp.n_diff_total, 'total diff == n_diff_total');
    fails = expect(fails, height(warm) == exp.N_WARMUP * exp.n_total_blocks, 'warmup total');
    fails = expect(fails, all(warm.is_target == 1), 'all warmups are target');

    % per-block same/diff match the stored split
    for b = 1:exp.n_total_blocks
        tb = test(test.block == b, :);
        ns = sum(strcmp(tb.couple, 'same'));  nd = sum(strcmp(tb.couple, 'diff'));
        fprintf('  block %d: %d same, %d diff\n', b, ns, nd);
        fails = expect(fails, ns == exp.same_per_block(b), sprintf('block %d same count', b));
        fails = expect(fails, nd == exp.diff_per_block(b), sprintf('block %d diff count', b));
    end

    % per-ISI test totals: each ISI should get exactly N_TRIALS test trials
    for iv = 1:numel(exp.ISI_VALUES)
        isi = exp.ISI_VALUES(iv);
        n_this = sum(test.isi == isi);
        fprintf('  ISI %4.0f ms: %d test trials\n', isi*1000, n_this);
        fails = expect(fails, n_this == exp.N_TRIALS, sprintf('ISI %.0f test == N_TRIALS', isi*1000));
    end

    %% ---- 2. Timing reconstruction (the key check) ----
    fprintf('\n--- Timing reconstruction ---\n');
    for b = 1:exp.n_total_blocks
        be = events(events.block == b, :);
        n  = height(be);
        onsets = zeros(n, 1);
        for k = 1:n-1
            dur = STIM * strcmp(be.event_type{k}, 'tone');
            onsets(k+1) = onsets(k) + dur + be.wait_after(k);
        end
        n_trials = n / 2;
        for tr = 1:n_trials
            r1 = 2*tr - 1;  r2 = 2*tr;
            isi = be.isi(r1);  rest = be.rest(r1);
            fails = expect_quiet(fails, abs((onsets(r2) - onsets(r1)) - isi) < tol, ...
                           sprintf('block %d trial %d: S2-S1==ISI', b, tr));
            if tr < n_trials
                got = onsets(r2+1) - onsets(r1);
                fails = expect_quiet(fails, abs(got - (isi + STIM + rest)) < tol, ...
                           sprintf('block %d trial %d: nextS1-S1==ISI+STIM+REST', b, tr));
            end
        end
    end
    fprintf('Reconstructed onsets for %d blocks (S2 lands at ISI; ITI = ISI+STIM+REST).\n', ...
            exp.n_total_blocks);

    %% ---- 3. Sequence constraint: no two 'diff' back to back ----
    fprintf('\n--- Sequence constraint (no two consecutive diff) ---\n');
    for b = 1:exp.n_total_blocks
        be = events(events.block == b & events.event_in_trial == 1, :);
        test_b = be(strcmp(be.trial_kind, 'test'), :);
        seq = test_b.couple;
        mrd = max_run_of(seq, 'diff');
        mrs = max_run_of(seq, 'same');
        fails = expect(fails, mrd <= 1, sprintf('block %d: no consecutive diff', b));
        fprintf('  block %d: max diff-run = %d, max same-run = %d, first test = %s\n', ...
                b, mrd, mrs, seq{1});
    end

    %% ---- 4. Metadata integrity ----
    fprintf('\n--- Metadata integrity ---\n');
    valid_trig = {'WARN', 'S2SM', 'S2DF'};
    fails = expect(fails, all(ismember(events.trigger, valid_trig)), 'all triggers valid');
    fails = expect(fails, all(strcmp(events.event_type, 'tone')), 'all events are tones');
    fails = expect(fails, all(events.s1_freq == exp.S1_FREQ), 's1_freq constant');

    % per-couple S2 frequency
    okS2 = true;
    for r = 1:height(events)
        if strcmp(events.couple{r}, 'same'); exp_s2 = exp.S1_FREQ; else; exp_s2 = exp.S2_FREQ; end
        if events.s2_freq(r) ~= exp_s2; okS2 = false; end
    end
    fails = expect(fails, okS2, 's2_freq matches couple');

    % No NaN in numeric metadata-bound columns
    numcols = {'block','trial','event_in_trial','is_target','isi','rest', ...
               'freq','s1_freq','s2_freq','wait_after'};
    anynan = false;
    for cc = 1:numel(numcols)
        if any(isnan(events.(numcols{cc}))); anynan = true; end
    end
    fails = expect(fails, ~anynan, 'no NaN in numeric columns');

    % CNAM milliseconds match ISI seconds
    okCnam = true;
    for r = 1:height(events)
        ms = sscanf(events.condition_name{r}, 'CNV_ISI%dms');
        if isempty(ms) || ms ~= round(events.isi(r)*1000); okCnam = false; end
    end
    fails = expect(fails, okCnam, 'CNAM ms matches ISI');

    %% ---- 5. REST distribution ----
    fprintf('\n--- Rest interval distribution ---\n');
    rests = s1.rest;
    fprintf('mean=%.3f  min=%.3f  max=%.3f  (target %.1f in [%.1f, %.1f])\n', ...
            mean(rests), min(rests), max(rests), exp.REST_MEAN, ...
            exp.REST_MEAN-exp.REST_JITTER, exp.REST_MEAN+exp.REST_JITTER);
    fails = expect(fails, abs(mean(rests) - exp.REST_MEAN) < 1e-6, 'rest mean == REST_MEAN');
    fails = expect(fails, min(rests) >= exp.REST_MEAN - exp.REST_JITTER - tol, 'rest >= floor');
    fails = expect(fails, max(rests) <= exp.REST_MEAN + exp.REST_JITTER + tol, 'rest <= ceil');

    %% ---- 6. Duration recompute ----
    fprintf('\n--- Duration ---\n');
    stim_dur = 0;
    for b = 1:exp.n_total_blocks
        be = events(events.block == b, :);
        bd = sum(be.wait_after) + sum(strcmp(be.event_type, 'tone')) * STIM;
        fprintf('  block %d (ISI %.0f ms): %.1f s (%.2f min)\n', b, be.isi(1)*1000, bd, bd/60);
        stim_dur = stim_dur + bd;
    end
    fprintf('Stimulus time: %.1f s (%.2f min).  Estimated total incl. breaks: %.2f min\n', ...
            stim_dur, stim_dur/60, exp.total_duration/60);
    if isfield(exp, 'stim_duration_est')
        fails = expect(fails, abs(stim_dur - exp.stim_duration_est) < 1e-6, ...
                       'stim duration matches stored estimate');
    end

    %% ---- Optional plot: couple sequence per block ----
    if show_plots
        try
            code_map = containers.Map({'same','diff'}, {1,2});
            maxlen = max(arrayfun(@(b) sum(events.block==b & events.event_in_trial==1), ...
                                  1:exp.n_total_blocks));
            M = nan(exp.n_total_blocks, maxlen);
            for b = 1:exp.n_total_blocks
                be = events(events.block==b & events.event_in_trial==1, :);
                for tr = 1:height(be); M(b, tr) = code_map(be.couple{tr}); end
            end
            figure('Name','CNV couple sequence','Color','w');
            imagesc(M); colormap([0.2 0.6 0.2; 0.7 0.1 0.1]);
            clim([1 2]); xlabel('trial (in block)'); ylabel('block');
            title('Couple sequence  (green=same/target, red=diff)');
            cb = colorbar('Ticks',[1.25 1.75],'TickLabels',{'same','diff'});
            cb.Label.String = 'couple';
        catch plotErr
            fprintf('(plot skipped: %s)\n', plotErr.message);
        end
    end

    %% ---- Verdict ----
    fprintf('\n=========================================================\n');
    if isempty(fails)
        fprintf(' ALL CHECKS PASSED.\n');
        fprintf('=========================================================\n');
    else
        fprintf(' %d CHECK(S) FAILED:\n', numel(fails));
        for k = 1:numel(fails); fprintf('   - %s\n', fails{k}); end
        fprintf('=========================================================\n');
        error('analyze_CNV_config:failed', '%d configuration check(s) failed.', numel(fails));
    end
end

%% ---- local helpers ----
function fails = expect(fails, cond, label)
    if cond
        fprintf('  [PASS] %s\n', label);
    else
        fprintf('  [FAIL] %s\n', label);
        fails{end+1} = label;
    end
end

function fails = expect_quiet(fails, cond, label)
    % Like expect() but only prints on failure (used for the hundreds of timing checks).
    if ~cond
        fprintf('  [FAIL] %s\n', label);
        fails{end+1} = label;
    end
end

function m = max_run_of(cellseq, val)
    % Longest run of consecutive entries equal to `val`.
    m = 0; cur = 0;
    for i = 1:numel(cellseq)
        if strcmp(cellseq{i}, val)
            cur = cur + 1; if cur > m; m = cur; end
        else
            cur = 0;
        end
    end
end
