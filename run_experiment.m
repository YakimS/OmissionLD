function results = run_experiment(config)
% RUN_EXPERIMENT Execute the oddball paradigm auditory experiment
%
% Usage:
%   config = experiment_config();
%   results = run_experiment(config);
%
% Or simply:
%   results = run_experiment(experiment_config());
%
% The experiment consists of:
%   1. Initial message: "Experiment starts now"
%   2. 5-second pause
%   3. Three blocks (random order), one for each standard tone
%   4. Each block has baseline trials followed by mixed trials
%
% Output:
%   results - Structure containing all trial data and timestamps
%
% Event markers are recorded for:
%   - Experiment start
%   - Block start
%   - Each tone presentation
%   - Each omission
%   - Block end
%   - Experiment end

    %% ===== INITIALIZATION =====

    fprintf('\n');
    fprintf('========================================\n');
    fprintf('ODDBALL PARADIGM EXPERIMENT\n');
    fprintf('========================================\n\n');

    % Create results directory if it doesn't exist
    if ~exist(config.results_dir, 'dir')
        mkdir(config.results_dir);
    end

    % Initialize results structure
    results = struct();
    results.config = config;
    results.start_time = datetime('now');
    results.blocks = [];
    results.all_events = [];

    %% ===== LOAD TONE FILES =====

    fprintf('Loading tone files...\n');
    tones = cell(1, 3);
    sample_rates = zeros(1, 3);

    for i = 1:3
        tone_path = fullfile(config.tone_dir, config.tone_files{i});
        if ~exist(tone_path, 'file')
            error('Tone file not found: %s\nPlease run generate_tones() first.', tone_path);
        end
        [tones{i}, sample_rates(i)] = audioread(tone_path);
        fprintf('  Loaded: %s\n', config.tone_files{i});
    end

    % Verify all sample rates are the same
    if length(unique(sample_rates)) > 1
        error('All tone files must have the same sample rate');
    end
    sample_rate = sample_rates(1);

    fprintf('Tones loaded successfully!\n\n');

    %% ===== PREPARE AUDIO PLAYER =====

    % Initialize audio player
    % We'll use audioplayer for precise timing control

    %% ===== EXPERIMENT START =====

    fprintf('========================================\n');
    fprintf('EXPERIMENT STARTS NOW\n');
    fprintf('========================================\n\n');

    % Record experiment start event
    event_start = create_event('experiment_start', 0, 0, '', 0);
    results.all_events = [results.all_events; event_start];

    % Initial pause
    fprintf('Initial pause: %d seconds...\n\n', config.initial_pause_duration);
    pause(config.initial_pause_duration);

    %% ===== DETERMINE BLOCK ORDER =====

    % Randomize the order of blocks (which tone is standard in each block)
    block_order = randperm(3);  % Random permutation of [1, 2, 3]

    fprintf('Block order (by standard tone):\n');
    for i = 1:3
        fprintf('  Block %d: %s\n', i, config.tone_labels{block_order(i)});
    end
    fprintf('\n');

    results.block_order = block_order;

    %% ===== RUN BLOCKS =====

    for block_num = 1:3
        % Determine which tone is standard for this block
        standard_tone_idx = block_order(block_num);

        % Determine which tones are deviants (the other two)
        deviant_tone_indices = setdiff(1:3, standard_tone_idx);

        fprintf('========================================\n');
        fprintf('BLOCK %d / 3\n', block_num);
        fprintf('Standard: %s\n', config.tone_labels{standard_tone_idx});
        fprintf('Deviants: %s, %s\n', ...
                config.tone_labels{deviant_tone_indices(1)}, ...
                config.tone_labels{deviant_tone_indices(2)});
        fprintf('========================================\n\n');

        % Run the block
        block_results = run_block(block_num, standard_tone_idx, ...
                                   deviant_tone_indices, tones, ...
                                   sample_rate, config);

        % Store block results
        results.blocks = [results.blocks; block_results];
        results.all_events = [results.all_events; block_results.events];

        % Inter-block interval (except after last block)
        if block_num < 3
            fprintf('\nInter-block interval: %d seconds...\n\n', config.inter_block_interval);
            pause(config.inter_block_interval);
        end
    end

    %% ===== EXPERIMENT END =====

    fprintf('\n========================================\n');
    fprintf('EXPERIMENT COMPLETED\n');
    fprintf('========================================\n\n');

    % Record experiment end event
    event_end = create_event('experiment_end', 0, 0, '', toc);
    results.all_events = [results.all_events; event_end];

    results.end_time = datetime('now');
    results.duration_seconds = seconds(results.end_time - results.start_time);

    fprintf('Total duration: %.2f minutes\n', results.duration_seconds / 60);

    %% ===== SAVE RESULTS =====

    results_path = fullfile(config.results_dir, config.results_filename);
    save(results_path, 'results');
    fprintf('Results saved to: %s\n\n', results_path);
end


%% ===== HELPER FUNCTION: RUN BLOCK =====

function block_results = run_block(block_num, standard_idx, deviant_indices, tones, sample_rate, config)
    % RUN_BLOCK Execute a single block of the experiment

    block_results = struct();
    block_results.block_number = block_num;
    block_results.standard_tone_idx = standard_idx;
    block_results.deviant_tone_indices = deviant_indices;
    block_results.events = [];

    % Start block timer
    tic;

    %% Generate trial sequence

    % Baseline trials (only standard tone)
    baseline_sequence = repmat(standard_idx, config.baseline_trials, 1);

    % Trials after baseline
    trials_after_baseline = config.trials_per_block - config.baseline_trials;

    % Calculate number of each trial type
    n_omissions = round(config.omission_percentage / 100 * trials_after_baseline);
    n_deviants = round(config.deviant_percentage / 100 * trials_after_baseline);
    n_deviants_per_type = round(n_deviants / 2);  % Split equally between two deviant types
    n_standards = trials_after_baseline - n_omissions - (n_deviants_per_type * 2);

    % Create trial sequence for non-baseline trials
    trial_sequence = [];

    % Add standard trials (using standard_idx)
    trial_sequence = [trial_sequence; repmat(standard_idx, n_standards, 1)];

    % Add deviant trials
    trial_sequence = [trial_sequence; repmat(deviant_indices(1), n_deviants_per_type, 1)];
    trial_sequence = [trial_sequence; repmat(deviant_indices(2), n_deviants_per_type, 1)];

    % Add omission trials (represented by 0)
    trial_sequence = [trial_sequence; zeros(n_omissions, 1)];

    % Randomize the non-baseline sequence
    trial_sequence = trial_sequence(randperm(length(trial_sequence)));

    % Combine baseline and randomized sequence
    full_sequence = [baseline_sequence; trial_sequence];

    %% Execute trials

    fprintf('Running %d trials...\n', length(full_sequence));

    for trial_num = 1:length(full_sequence)
        trial_start_time = toc;

        tone_idx = full_sequence(trial_num);

        % Determine trial type
        if tone_idx == 0
            trial_type = 'omission';
            tone_label = 'Omission';
        elseif tone_idx == standard_idx
            trial_type = 'standard';
            tone_label = config.tone_labels{tone_idx};
        else
            trial_type = 'deviant';
            tone_label = config.tone_labels{tone_idx};
        end

        % Record event
        event = create_event(trial_type, block_num, trial_num, tone_label, trial_start_time);
        block_results.events = [block_results.events; event];

        % Play tone (or omit)
        if tone_idx > 0
            % Play the tone
            player = audioplayer(tones{tone_idx}, sample_rate);
            play(player);
            % Don't wait for playback to finish - continue to maintain SOA timing
        end

        % Progress indicator
        if config.verbose && mod(trial_num, config.progress_interval) == 0
            fprintf('  Trial %d/%d completed\n', trial_num, length(full_sequence));
        end

        % Wait for SOA (inter-stimulus interval)
        % Calculate time to wait to maintain precise SOA
        trial_elapsed = toc - trial_start_time;
        soa_seconds = config.stimulus_onset_asynchrony_ms / 1000;
        wait_time = soa_seconds - trial_elapsed;

        if wait_time > 0
            pause(wait_time);
        end
    end

    block_duration = toc;
    block_results.duration_seconds = block_duration;

    fprintf('Block %d completed in %.2f seconds\n', block_num, block_duration);

    % Summary statistics
    n_standard_played = sum(strcmp({block_results.events.trial_type}, 'standard'));
    n_deviant_played = sum(strcmp({block_results.events.trial_type}, 'deviant'));
    n_omission_played = sum(strcmp({block_results.events.trial_type}, 'omission'));

    fprintf('  Standard: %d, Deviant: %d, Omission: %d\n', ...
            n_standard_played, n_deviant_played, n_omission_played);
end


%% ===== HELPER FUNCTION: CREATE EVENT =====

function event = create_event(trial_type, block_num, trial_num, tone_label, timestamp)
    % CREATE_EVENT Create an event marker structure

    event = struct();
    event.trial_type = trial_type;
    event.block_number = block_num;
    event.trial_number = trial_num;
    event.tone_label = tone_label;
    event.timestamp = timestamp;
    event.datetime = datetime('now');
end
