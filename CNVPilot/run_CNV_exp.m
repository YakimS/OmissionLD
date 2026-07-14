% run_CNV_exp.m
% Runs the auditory CNV experiment from a pre-generated configuration (Phase 2).
%
% Each trial: S1 (warning tone) -> S1-S2 interval (ISI) -> S2 (same tone or different
% tone). NO response is collected online - the subject responds to the target couple
% ('same') with an EYE MOVEMENT, analysed offline from EOG. A silent nature movie
% plays externally throughout.
%
% Usage:
%   1. First run config_CNV_exp.m to generate the configuration .mat file.
%   2. Set CONFIG_FILE below, start the external silent movie, then run this script.
%
% Controls (lab mode): press SPACE (any key) to advance a self-paced break;
%                      press ESC during an inter-trial rest to abort gracefully.

clear all;
close all;
clc;

isLabComputer = false;
if isLabComputer; Priority(1); end

SUBJECT_ID  = 'S00';                       % For logging only
CONFIG_FILE = 'Config-CNVPilot-demo.mat';

% Audio device (lab). Falls back to the system default if not found.
AUDIO_DEVICE_NAME = 'Headphones (Realtek(R) Audio)';
AUDIO_API_NAME    = 'MME';

%% Load experiment configuration
loaded   = load(CONFIG_FILE);
exp_meta = loaded.exp_meta;
events   = exp_meta.events_table;
STIM     = exp_meta.STIM_DURATION;

% NetStation configuration
USE_NETSTATION  = isLabComputer & true;
NETSTATION_HOST = '10.10.10.42';
NETSTATION_PORT = 55513;
NETSTATION_NTP  = '10.10.10.51';

% Keyboard (lab only)
escKey = [];
if isLabComputer
    KbName('UnifyKeyNames');
    escKey = KbName('ESCAPE');
end

fprintf('\n===========================================\n');
fprintf('   AUDITORY CNV EXPERIMENT\n');
fprintf('   Subject: %s | Config: %s\n', SUBJECT_ID, CONFIG_FILE);
fprintf('   S1=%d Hz, S2=%d Hz | %d blocks | ~%.1f min\n', ...
        exp_meta.S1_FREQ, exp_meta.S2_FREQ, exp_meta.n_total_blocks, exp_meta.total_duration/60);
fprintf('===========================================\n');

pamaster     = [];
pahandleplay = [];

try
    %% Connect to NetStation
    if isLabComputer && USE_NETSTATION
        fprintf('Connecting to NetStation at %s:%d...\n', NETSTATION_HOST, NETSTATION_PORT);
        [nsstatus, nserror] = NetStation('Connect', NETSTATION_HOST, NETSTATION_PORT);
        if nsstatus ~= 0
            error('CNV:netstation', 'Could not connect to NetStation: %s', nserror);
        end
        % Use BOTH clock-sync methods (as in Areej's lab script): the legacy ECI
        % 'Synchronize' first, then NTP last so the more accurate NTP offset is the
        % one active for event timestamps.
        NetStation('Synchronize', 10);                      % legacy ECI sync
        WaitSecs(0.5);
        NetStation('GetNTPSynchronize', NETSTATION_NTP);    % NTP sync (primary, more accurate)
        WaitSecs(1);
        fprintf('NetStation connected and synchronized (Synchronize + NTP).\n');
    end

    %% Sound setup
    if isLabComputer
        InitializePsychSound(1);

        OutputDevice = [];
        devices = PsychPortAudio('GetDevices');
        for ii = 1:length(devices)
            if strcmp(devices(ii).DeviceName, AUDIO_DEVICE_NAME) && ...
               strcmp(devices(ii).HostAudioAPIName, AUDIO_API_NAME)
                OutputDevice = devices(ii).DeviceIndex;
            end
        end
        if isempty(OutputDevice)
            warning('Audio device "%s" (%s) not found; using system default.', ...
                    AUDIO_DEVICE_NAME, AUDIO_API_NAME);
        end

        pamaster = PsychPortAudio('Open', OutputDevice, 1+8, 3, exp_meta.SAMPLE_RATE, 1);
        PsychPortAudio('Start', pamaster, 0, 0, 1);
        pahandleplay = PsychPortAudio('OpenSlave', pamaster);
        PsychPortAudio('Volume', pahandleplay, 1);
    end

    %% Map each frequency to its waveform
    beep_sounds = containers.Map('KeyType', 'double', 'ValueType', 'any');
    beep_sounds(exp_meta.S1_FREQ) = exp_meta.beep_S1;
    beep_sounds(exp_meta.S2_FREQ) = exp_meta.beep_S2;

    %% Start recording + experiment-start trigger
    if USE_NETSTATION
        NetStation('StartRecording');
        NetStation('Event', 'STRT', GetSecs, 0.001, 'subj', SUBJECT_ID);
        WaitSecs(2);
    end

    %% Main loop over events (two rows per trial: S1 then S2)
    current_block = 0;
    trial_tS2     = NaN;   % absolute deadline for the current trial's S2 (lab scheduling)

    for i = 1:height(events)
        event = events(i, :);

        % ---- Block transition ----
        if event.block ~= current_block
            if current_block > 0
                fprintf('\nBlock %d complete.\n', current_block);
                selfpaced_break(current_block, exp_meta.n_total_blocks, isLabComputer);
            end
            current_block = event.block;
            cnam = event.condition_name{1};

            fprintf('\n===========================================\n');
            fprintf('Block %d/%d: %s\n', current_block, exp_meta.n_total_blocks, cnam);
            fprintf('===========================================\n');

            if USE_NETSTATION
                NetStation('GetNTPSynchronize', NETSTATION_NTP);
                NetStation('Event', 'BGIN', GetSecs, 0.001, ...
                           'BNUM', current_block, 'CNAM', cnam, ...
                           'ISI', event.isi, 'PRMS', exp_meta.params_string);
            end
        end

        if exp_meta.SHOW_TRIAL_OUTPUT && event.event_in_trial == 1
            fprintf('Trial %2d [%-4s %s]: ', event.trial, event.couple{1}, event.trial_kind{1});
        end

        % ---- Event 1: S1 (warning tone, identical across couples) ----
        if event.event_in_trial == 1
            sound_data   = beep_sounds(event.freq);
            startTime_S1 = play_beep(pahandleplay, sound_data, isLabComputer, exp_meta.SAMPLE_RATE, 0);
            trial_tS2    = startTime_S1 + event.isi;   % absolute S2 onset deadline (lab)

            if exp_meta.SHOW_TRIAL_OUTPUT; fprintf('%s ', event.display_str{1}); end
            send_event(USE_NETSTATION, 'WARN', startTime_S1, event, STIM);

            if ~isLabComputer
                do_wait(event.wait_after, isLabComputer);   % fallback relative pacing
            end
            continue;
        end

        % ---- Event 2: S2 (same tone or different tone) ----
        sound_data = beep_sounds(event.freq);
        if isLabComputer
            startTime_S2 = play_beep(pahandleplay, sound_data, isLabComputer, ...
                                     exp_meta.SAMPLE_RATE, trial_tS2);   % scheduled at ISI
        else
            startTime_S2 = play_beep(pahandleplay, sound_data, isLabComputer, ...
                                     exp_meta.SAMPLE_RATE, 0);
        end
        if exp_meta.SHOW_TRIAL_OUTPUT; fprintf('%s ', event.display_str{1}); end
        send_event(USE_NETSTATION, event.trigger{1}, startTime_S2, event, STIM);

        % ---- Inter-trial rest (bring us to the next S1 onset) ----
        is_last_overall  = (i == height(events));
        next_is_newblock = ~is_last_overall && (events.block(i+1) ~= current_block);

        if isLabComputer
            if ~is_last_overall && ~next_is_newblock
                t_nextS1 = trial_tS2 + STIM + event.rest;   % next S1 = S2 onset + STIM + REST
                if wait_until_esc(t_nextS1, escKey)
                    error('CNV:abort', 'ESC pressed - aborting.');
                end
            end
        else
            do_wait(event.wait_after, isLabComputer);        % fallback relative pacing
        end

        if exp_meta.SHOW_TRIAL_OUTPUT; fprintf('\n'); end
    end

    %% Experiment-end trigger
    fprintf('\nBlock %d complete.\n', current_block);
    if USE_NETSTATION
        WaitSecs(2);
        NetStation('Event', 'SEND', GetSecs, 0.001, 'subj', SUBJECT_ID);
        NetStation('StopRecording');
    end
    fprintf('\n===========================================\n');
    fprintf('EXPERIMENT COMPLETE\n');
    fprintf('===========================================\n');

catch ME
    fprintf(2, '\n*** Experiment stopped: %s ***\n', ME.message);
    if USE_NETSTATION
        try
            NetStation('Event', 'SEND', GetSecs, 0.001, 'abrt', 1);
            NetStation('StopRecording');
        catch
        end
    end
end

%% Cleanup (always)
if USE_NETSTATION
    try NetStation('Disconnect'); catch; end   % release the session so the next run's Connect is clean
end
if isLabComputer
    % close BOTH slave and master to free the audio device
    try PsychPortAudio('Stop', pamaster, 0); catch; end
    try PsychPortAudio('Close'); catch; end    % no handle -> closes master + all slaves
end

%% ===================== Helper functions =====================

function send_event(use_ns, code, when, event, stim)
    % Send one NetStation event with the full CNV metadata block. No-op unless
    % NetStation is enabled. Metadata keys mirror the events-table columns.
    if ~use_ns; return; end
    NetStation('Event', code, when, max(stim, 0.001), ...
        'BNUM', event.block, ...
        'TNUM', event.trial, ...
        'EVNT', event.event_in_trial, ...
        'COUP', event.couple{1}, ...
        'TTYP', event.trial_kind{1}, ...
        'TRGT', event.is_target, ...
        'ISI',  event.isi, ...
        'REST', event.rest, ...
        'S1FQ', event.s1_freq, ...
        'S2FQ', event.s2_freq, ...
        'CNAM', event.condition_name{1}, ...
        'CTYP', event.ctyp{1}, ...
        'WAIT', event.wait_after, ...
        'DISP', event.display_str{1});
end

function startTime = play_beep(pahandle, beep_sound, isLabComputer, sampleRate, when)
    % Play a tone. 'when' is an absolute PsychPortAudio start time (0 = now).
    % Returns the actual onset timestamp. Blocks until playback ends.
    if nargin < 5; when = 0; end
    if isLabComputer
        PsychPortAudio('FillBuffer', pahandle, beep_sound);
        startTime = PsychPortAudio('Start', pahandle, 1, when, 1);
        PsychPortAudio('Stop', pahandle, 1);
    else
        startTime = posixtime(datetime('now'));
        sound(beep_sound, sampleRate);
        pause(numel(beep_sound) / sampleRate);   % block for the tone (sound() is async) so
    end                                          % dry-run pacing matches the lab timeline
end

function abort = wait_until_esc(t_end, escKey)
    % Wait (lab) until the absolute time t_end, polling ESC every ~5 ms.
    % Returns true if ESC was pressed.
    abort = false;
    while GetSecs < t_end - 0.010
        [~, ~, keyCode] = KbCheck;
        if keyCode(escKey); abort = true; return; end
        WaitSecs(0.005);
    end
    WaitSecs('UntilTime', t_end);   % precise final sliver
end

function do_wait(duration, isLabComputer)
    if duration <= 0; return; end
    if isLabComputer
        WaitSecs(duration);
    else
        pause(duration);
    end
end

function selfpaced_break(block_num, n_blocks, isLabComputer)
    fprintf('\n--- BREAK --- Block %d of %d done. ', block_num, n_blocks);
    if isLabComputer
        fprintf('Press any key (SPACE) to start the next block...\n');
        KbReleaseWait;
        KbStrokeWait;
    else
        input('Press Enter to start the next block...', 's');
    end
end
